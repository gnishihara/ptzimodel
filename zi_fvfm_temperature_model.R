# Run the Caulerpa PAM -- Temperature Model.
# Greg Nishihara
# 2020 Sep 17

library(tidyverse)
library(readxl)
library(rstan)
library(brms)
library(tidybayes)
library(bayesplot)

# source("theme_iris.R") # I don't have this script so no theme for me.

# My system is setup for a Japanese locale so any plots with 
# a time class will be in Japanese. Change it to American English.
Sys.setlocale("LC_TIME", "en_US.UTF-8") 

# Load the data set -------------------------------------------------------------
# Name of the relative path of the folder 
# This will be different for each system.
FOLDER = "~/Dropbox/MS_Iris_Greg_Terada/2_Terada_Caulerpa_lentillifera_study/Dataset_for_Figures/Fig_1_PAM_Temperature/"
# Get file and path names.
FNAMES = dir(FOLDER, pattern = "xlsx", full.names = TRUE) 

# The names of the sheets I want to read from the excel workbook.
sheets0 = c("Data_for_analysis_24h", 
            "Data_for_analysis_48h",
            "Data_for_analysis_72h") 
sheets50 = c("Data_for_analysis_50µmol_24h",
             "Data_for_analysis_50µmol_48h",
             "Data_for_analysis_50µmol_72h")

# Build the tibble with the file names and the sheet names.
alldata = tibble(fnames =  rep(FNAMES, each = 3), 
                 sheets = c(sheets0, sheets50))

# Read all the files at once.
# The order of the arguments to the read_xlsx() function is important.
alldata = alldata %>% mutate(data = map2(fnames, sheets, read_xlsx)) 

# Clean up the data
alldata =  alldata %>% 
  mutate(Hour = str_extract(sheets, "24h|48h|72h")) %>% 
  mutate(Light = str_extract(fnames, "[0-9]{1,2}µmol")) %>% 
  select(Light, Hour, data) %>% 
  unnest(data) %>% 
  rename(Temperature = Temp, FvFm = Yield)

# Summary statistics
alldata_summary = alldata %>% group_by(Temperature, Light, Hour) %>% 
  summarise(mean = mean(FvFm),
            n = length(FvFm),
            sd = sd(FvFm),
            se = sd(FvFm) / sqrt(n-1),
            l95 = mean - 1.96*sd,
            u95 = mean + 1.96*sd)

# Take a look at the data.
# So the problem is the zeros occurring below 16 °C and above 32 °C.
temperature_breaks = alldata %>% pull(Temperature) %>% unique()
alldata %>% ggplot() + 
  geom_point(aes(x = Temperature, y = FvFm, color = Hour),  
             position = position_jitter(0.1))  +
  scale_x_continuous(breaks = temperature_breaks) +
  facet_grid(rows = vars(Light))


# Bayesian inference with BRMS starts here -------------------------------------
### Creating model as a string and pass to stanvar()
stan_function = "
  real fvfmmodel (real ps, real ha, real eta, real ktopt, real temperature) {
    real inverse_kelvin = 1.0 / (temperature + 273.15);
    real gas_constant = 8.314/1000.0;
    real hd = ha * eta;
    real x = (1.0 / ktopt - inverse_kelvin);
    real numerator = ha / gas_constant * x;
    real denominator = hd / gas_constant * x;
    return ps * hd * exp(numerator) / (hd - ha * (1.0 - exp(denominator)));
  }
"
STANVARS = stanvar(scode = stan_function, block = "functions")

# Build the brms formula. 
# IMPORTANT: the zero-inflation rate portion of the model must be a 
# linear model. I am not sure if the writing the zi model as 
# zi ~ Temperature + Temperature^2 was working properly, so I created a new
# variable Temperature2 which is the square of the Temperature.

brmsmodel = brmsformula(FvFm ~ fvfmmodel(PS, HA, ET, KT, Temperature),
                        PS ~ 1,
                        HA ~ 1,
                        ET ~ 1,
                        KT ~ 1,
                        zi ~ Temperature + Temperature2,
                        nl = TRUE,
                        family = brms::zero_inflated_beta(link = "identity",
                                                          link_phi = "log",
                                                          link_zi= "logit"))

# Add the squared Temperature to the dataset.
alldata = alldata %>% mutate(Temperature2 = Temperature^2)

# Build the prior here ---------------------------------------------------------
# priors = read_csv("../Saccharina_japonica_study/notes/priors_FvFm.csv")
# The priors are based on the parameter estimates from past papers.
priors = read_csv("priors_FvFm.csv")
priors = priors %>% group_by(parameter) %>% 
  summarise(mean = mean(value), sd = sd(value))

priors = priors %>%
  mutate(mean = ifelse(parameter == "topt", mean + 273.15, mean)) %>%
  mutate(prior = str_glue("normal({mean}, {sd})")) %>%
  select(parameter, prior)

get_prior(brmsmodel, data = alldata) # Run this to see the defaults.

# Build the prior distribution here.
# The prior distribution for the zero-inflation part of the model are the
# default priors.
# Note how there is a lower bound (lb) for the ET parameter and the
# KT parameter.

PRIORS = 
  set_prior(priors %>% filter(parameter == "eta") %>% pull(prior),
            class = "b", lb = 1, nlpar = "ET") +
  set_prior(priors %>% filter(parameter == "topt") %>% pull(prior),
            class = "b", lb = 273.15, nlpar = "KT") +
  set_prior(priors %>% filter(parameter == "ymax") %>% pull(prior),
            class = "b", lb = 0, nlpar = "PS") +
  set_prior(priors %>% filter(parameter == "ha") %>% pull(prior),
            class = "b", lb = 0, nlpar = "HA") +
  set_prior("normal(0, 5)", class = "b", dpar = "zi") +
  set_prior("normal(0, 5)", class = "Intercept", dpar = "zi") +
  set_prior("gamma(0.01,0.01)", class = "phi") 

# I wanted to look a the stan code that was generated by brms.
make_stancode(brmsmodel, data = alldata, 
              stanvars = STANVARS, prior = PRIORS, 
              save = "zimodel.stan")

# Sample from the posterior distribution ---------------------------------------
# Parameters fo the sampler
ITER = 10
ADAPT_DELTA = 0.80
MAX_TREEDEPTH = 10
CHAINS = 4
CORES = CHAINS
CONTROL = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH)

# file.remove("zimodel.rds")
# Preliminary fit to create the brm object.
# Save the object to the file zimodel.rds using the file argument.
if(!file.exists("zimodel.rds")) {
  stime = lubridate::now()
  brmout = brm(brmsmodel,
               data = alldata %>% filter(str_detect(Hour, "4")),
               stanvars = STANVARS, prior = PRIORS,
               iter = ITER, chains = CHAINS, cores = CORES, 
               seed = 2020, control = CONTROL, file = "zimodel")
  etime = lubridate::now()
  etime - stime # How long did it take?
} else {
  brmout = readRDS("zimodel.rds")
}


# This is a function to sample from the posterior distribution after the model 
# object was compiled with brm(). The function is passed to the map() function
# to run the sampler on each subset of the data.
posterior_samples = function(dataset, ADAPT_DELTA, MAX_TREEDEPTH, ITER) {
  CONTROL = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH)
  update(brmout, newdata = dataset,
         stanvars = STANVARS, prior = PRIORS,
         iter = ITER, chains = CHAINS, cores = CORES, seed = 2020,
         control = CONTROL)
}
stime = lubridate::now() # This is here to calculate sampling time.
outdata = alldata %>% group_nest(Hour, Light) %>% 
  mutate(brmout = map(data, posterior_samples,
                      ADAPT_DELTA = 0.99999,
                      MAX_TREEDEPTH = 15,
                      ITER = 5000))
etime = lubridate::now()
etime - stime # How long did it take?

save(outdata, file = "Caulerpa_brms_samples_20200917.RData")
