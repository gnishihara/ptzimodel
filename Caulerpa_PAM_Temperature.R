## Caulerpa lentillifera PAM Temperature Analysis
## Greg Nishihara
## 2020 SEP 07
## 


library(tidyverse)
library(readxl)
library(rstan)
library(brms)
library(tidybayes)
library(bayesplot)
library(ggpubr)
library(knitr)
library(kableExtra)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gtable)
library(xtable)
library(grid)
library(gridExtra)


# source("theme_iris.R") # I dont have this script so no theme for me.

# I need to set the time code for my system.
Sys.setlocale("LC_TIME", "en_US.UTF-8") # アメリカ英語に設定


# Load the dataset.

FOLDER = "~/Dropbox/MS_Iris_Greg_Terada/2_Terada_Caulerpa_lentillifera_study/Dataset_for_Figures/Fig_1_PAM_Temperature/"
FNAMES = dir(FOLDER, pattern = "xlsx", full.names = TRUE)
sheets0 = c("Data_for_analysis_24h","Data_for_analysis_48h","Data_for_analysis_72h") 
sheets50 = c("Data_for_analysis_50µmol_24h","Data_for_analysis_50µmol_48h","Data_for_analysis_50µmol_72h")
Caulerpa_0_24 = read_xlsx(path = FNAMES[1], sheet = sheets0[1])  %>% rename(Temperature = Temp, FvFm = Yield)
Caulerpa_0_48 = read_xlsx(path = FNAMES[1], sheet = sheets0[2])  %>% rename(Temperature = Temp, FvFm = Yield)
Caulerpa_0_72 = read_xlsx(path = FNAMES[1], sheet = sheets0[3])  %>% rename(Temperature = Temp, FvFm = Yield)
Caulerpa_50_24= read_xlsx(path = FNAMES[2], sheet = sheets50[1]) %>% rename(Temperature = Temp, FvFm = Yield)
Caulerpa_50_48= read_xlsx(path = FNAMES[2], sheet = sheets50[2]) %>% rename(Temperature = Temp, FvFm = Yield)
Caulerpa_50_72= read_xlsx(path = FNAMES[2], sheet = sheets50[3]) %>% rename(Temperature = Temp, FvFm = Yield)


# Caulerpa_0_24=read.csv("../Clentillifera/Fig_1_PAM_Temperature/Caulerpa_0PAR_8-40C_24h.csv")
# Caulerpa_0_48=read.csv("../Clentillifera/Fig_1_PAM_Temperature/Caulerpa_0PAR_8-40C_48h.csv")
# Caulerpa_0_72=read.csv("../Clentillifera/Fig_1_PAM_Temperature/Caulerpa_0PAR_8-40C_72h.csv")
# Caulerpa_50_24=read.csv("../Clentillifera/Fig_1_PAM_Temperature/Caulerpa_50PAR_8-40C_24h.csv")
# Caulerpa_50_48=read.csv("../Clentillifera/Fig_1_PAM_Temperature/Caulerpa_50PAR_8-40C_48h.csv")
# Caulerpa_50_72=read.csv("../Clentillifera/Fig_1_PAM_Temperature/Caulerpa_50PAR_8-40C_72h.csv")

Caulerpa_0_24_mean = Caulerpa_0_24 %>% 
  group_by(Temperature) %>% 
  summarise(mean = mean(FvFm),
            sd = sd(FvFm),
            se = sd(FvFm) / sqrt(n()),
            l95 = mean - 1.96*sd,
            u95 = mean + 1.96*sd)
Caulerpa_0_24_mean

Caulerpa_0_48_mean = Caulerpa_0_48 %>% 
  group_by(Temperature) %>% 
  summarise(mean = mean(FvFm),
            sd = sd(FvFm),
            se = sd(FvFm) / sqrt(n()),
            l95 = mean - 1.96*sd,
            u95 = mean + 1.96*sd)
Caulerpa_0_48_mean

Caulerpa_0_72_mean = Caulerpa_0_72 %>% 
  group_by(Temperature) %>% 
  summarise(mean = mean(FvFm),
            sd = sd(FvFm),
            se = sd(FvFm) / sqrt(n()),
            l95 = mean - 1.96*sd,
            u95 = mean + 1.96*sd)
Caulerpa_0_72_mean

Caulerpa_50_24_mean = Caulerpa_50_24 %>% 
  group_by(Temperature) %>% 
  summarise(mean = mean(FvFm),
            sd = sd(FvFm),
            se = sd(FvFm) / sqrt(n()),
            l95 = mean - 1.96*sd,
            u95 = mean + 1.96*sd)
Caulerpa_50_24_mean

Caulerpa_50_48_mean = Caulerpa_50_48 %>% 
  group_by(Temperature) %>% 
  summarise(mean = mean(FvFm),
            sd = sd(FvFm),
            se = sd(FvFm) / sqrt(n()),
            l95 = mean - 1.96*sd,
            u95 = mean + 1.96*sd)
Caulerpa_50_48_mean

Caulerpa_50_72_mean = Caulerpa_50_72 %>% 
  group_by(Temperature) %>% 
  summarise(mean = mean(FvFm),
            sd = sd(FvFm),
            se = sd(FvFm) / sqrt(n()),
            l95 = mean - 1.96*sd,
            u95 = mean + 1.96*sd)
Caulerpa_50_72_mean

Cau_0_24 = tibble(Caulerpa_0_24) %>% 
  select(FvFm, Temperature)
Cau_0_48 = tibble(Caulerpa_0_48) %>% 
  select(FvFm, Temperature)
Cau_0_72 = tibble(Caulerpa_0_72) %>% 
  select(FvFm, Temperature)
Cau_50_24 = tibble(Caulerpa_50_24) %>% 
  select(FvFm, Temperature)
Cau_50_48 = tibble(Caulerpa_50_48) %>% 
  select(FvFm, Temperature)
Cau_50_72 = tibble(Caulerpa_50_72) %>% 
  select(FvFm, Temperature)


alldata = tibble(light = c(0,0,0,50,50,50),
                 hour = c(24,48,72, 24,48,72),
                 data = list(Cau_0_24, Cau_0_48, Cau_0_72, Cau_50_24, Cau_50_48, Cau_50_72)) %>% 
  unnest(cols = c(data)) 

ggplot(alldata) + 
  geom_point(aes(x = Temperature, y = FvFm), 
             position = "jitter") +
  facet_grid(cols = vars(hour),
             rows = vars(light))


### Creating model as string ###
stan_function =
  "real fvfmmodel (real ps, real ha, real eta, real ktopt, real temperature) {
  real invkelvin = 1.0 / (temperature + 273.15);
  real gas_constant = 8.314/1000.0;
  real tmp0 = (1.0 / ktopt - invkelvin);
  real tmp1 = ha / gas_constant * tmp0 ;
  real tmp2 = (ha * eta) / gas_constant * tmp0;
  tmp1 = (ha * eta) * exp(tmp1);
  tmp2 = 1 - exp(tmp2);
  tmp2 = (ha * eta) - ha * tmp2;
  return ps * tmp1 / tmp2;
  }
"

stanvars = stanvar(scode = stan_function, block = "functions")

### brms model ###
brmsmodel = brmsformula(FvFm ~ fvfmmodel(PS, HA, ET, KT, Temperature),
                        PS ~ 1,
                        HA ~ 1,
                        ET ~ 1,
                        KT ~ 1,
                        zi ~ zTemperature,
                        nl = TRUE,
                        family = brms::zero_inflated_beta(link = "identity",
                                                          link_phi = "log",
                                                          link_zi= "logit"))

alldata %>% filter(near(light, 0)) %>%
  mutate(zTemperature = (Temperature - median(Temperature))^2) %>% 
  get_prior(brmsmodel, data = .)

### build priors ###
# priors = read_csv("../Saccharina_japonica_study/notes/priors_FvFm.csv")
priors = read_csv("priors_FvFm.csv")
priors2 = priors %>% 
  group_by(parameter) %>% 
  summarise(mean = mean(value),
            sd = sd(value))

prs = priors2 %>%
  mutate(mean = ifelse(parameter == "topt", mean + 273.15, mean)) %>%
  mutate(prior = str_glue("normal({mean}, {2*sd})")) %>%
  select(parameter, prior)

PRIORS = set_prior(prs %>% filter(parameter == "eta") %>% pull(prior),
                   class = "b", lb = 0, nlpar = "ET") +
  set_prior(prs %>% filter(parameter == "topt") %>% pull(prior),
            class = "b", lb = 0, nlpar = "KT") +
  set_prior(prs %>% filter(parameter == "ymax") %>% pull(prior),
            class = "b", lb = 0, nlpar = "PS") +
  set_prior(prs %>% filter(parameter == "ha") %>% pull(prior),
            class = "b", lb = 0, nlpar = "HA") 


### Sample from the posterior distribution ###
RESAMPLE = F

ITER = 10000
ADAPT_DELTA = 0.999999
MAX_TREEDEPTH = 14

# Preliminary fit.

stime = lubridate::now()
brmout = brm(brmsmodel,
             data = Cau_0_24 %>% mutate(zTemperature = (Temperature - median(Temperature))^2),
             stanvars = stanvars, prior = PRIORS,
             iter = 100, chains = 1, cores = 1, seed = 2020)


fit_and_save = function(dataset, fname) {
     update(brmout, newdata = dataset,
                  stanvars = stanvars, prior = PRIORS,
                  iter = ITER, chains = 4, cores = 4, seed = 2020,
                  control = list(adapt_delta = ADAPT_DELTA,
                                 max_treedepth = MAX_TREEDEPTH))
}
expose_functions(brmout, vectorize = TRUE)


alldata2 = alldata %>% 
  mutate(zTemperature = (Temperature - median(Temperature))^2) %>% 
  group_nest(light, hour) %>% 
  mutate(fname = paste0("FvFm_Caulerpa_", light, "_", hour, ".RData")) %>% 
  mutate(brmout = map2(data, fname, fit_and_save))

etime = lubridate::now()

save(alldata2,stime, etime, file = "ALLDATA2.RData")



theplots = alldata2 %>% 
  mutate(summary = map(brmout, summary)) %>%
  mutate(ppc_dens_overlay = map(brmout, function(bfit) {
    y = bfit$data %>% pull(FvFm) 
    yrep = posterior_predict(bfit, nsamples = 500)
    ppc_dens_overlay(y, yrep)
  })) %>% 
  mutate(mcmc_rank_overlay = map(brmout, function(bfit){
    posterior_out = as.array(bfit)
    mcmc_rank_overlay(posterior_out)      
  }))

theplots %>% pull(ppc_dens_overlay) %>% 
  ggpubr::ggarrange(plotlist = .)

toplot = alldata2 %>% 
  mutate(predict = map(brmout, function(bfit) {
    bfit$data %>% 
      add_predicted_draws(bfit) %>% 
      group_by(Temperature) %>% 
      mean_hdci(FvFm)
  })) %>% 
  mutate(expect = map(brmout, function(bfit) {
    bfit$data %>% 
      add_linpred_draws(bfit) %>% 
      group_by(Temperature) %>% 
      mean_hdci(FvFm)
  }))





xlabF = expression("Temperature" ~ ({}^degree*C ))
ylabF = "Maximum Quantum Yield (Fv/Fm)" 

ggplot() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature,
                  fill = hour),
              alpha = 0.2,
              data = toplot %>% unnest(predict)) +
  geom_line(aes(x = Temperature, y = FvFm, color = hour),
            data = toplot %>% unnest(expect)) +
  geom_point(aes(x = Temperature, y = FvFm, color = hour),
             data = toplot %>% unnest(data)) +
  facet_grid(rows = vars(hour),
             cols = vars(light))
  


# Original fitting procedure.
# if(!file.exists("FvFm_Caulerpa_0_24.RData")) {
#   brmout1 = brm(brmsmodel, data = Cau_0_24,
#                 stanvars = stanvars, prior = PRIORS,
#                 iter = 5000, chains = 4, cores = 4, seed = 2020,
#                 control = list(adapt_delta = 0.99, max_treedepth = 10))
#   save(brmout1, file="FvFm_Caulerpa_0_24.RData")
# } else {
#   load(file = "FvFm_Caulerpa_0_24.RData", verbose=T)
# }
# 
# if(!file.exists("FvFm_Caulerpa_0_48.RData")) {
#   brmout2 = brm(brmsmodel, data = Cau_0_48,
#                 stanvars = stanvars, prior = PRIORS,
#                 iter = 5000, chains = 4, cores = 4, seed = 2020,
#                 control = list(adapt_delta = 0.9999, max_treedepth = 10))
#   save(brmout2, file="FvFm_Caulerpa_0_48.RData")
# } else {
#   load(file = "FvFm_Caulerpa_0_48.RData", verbose=T)
# }
# 
# if(!file.exists("FvFm_Caulerpa_0_72.RData")) {
#   brmout3 = brm(brmsmodel, data = Cau_0_72,
#                 stanvars = stanvars, prior = PRIORS,
#                 iter = 5000, chains = 4, cores = 4, seed = 2020,
#                 control = list(adapt_delta = 0.99, max_treedepth = 10))
#   save(brmout3, file="FvFm_Caulerpa_0_72.RData")
# } else {
#   load(file = "FvFm_Caulerpa_0_72.RData", verbose=T)
# }
# 
# if(!file.exists("FvFm_Caulerpa_50_24.RData")) {
#   brmout4 = brm(brmsmodel, data = Cau_50_24,
#                 stanvars = stanvars, prior = PRIORS,
#                 iter = 5000, chains = 4, cores = 4, seed = 2020,
#                 control = list(adapt_delta = 0.99, max_treedepth = 10))
#   save(brmout4, file="FvFm_Caulerpa_50_24.RData")
# } else {
#   load(file = "FvFm_Caulerpa_50_24.RData", verbose=T)
# }
# 
# if(!file.exists("FvFm_Caulerpa_50_48.RData")) {
#   brmout5 = brm(brmsmodel, data = Cau_50_48,
#                 stanvars = stanvars, prior = PRIORS,
#                 iter = 5000, chains = 4, cores = 4, seed = 2020,
#                 control = list(adapt_delta = 0.99, max_treedepth = 10))
#   save(brmout5, file="FvFm_Caulerpa_50_48.RData")
# } else {
#   load(file = "FvFm_Caulerpa_50_48.RData", verbose=T)
# }
# 
# if(!file.exists("FvFm_Caulerpa_50_72.RData")) {
#   brmout6 = brm(brmsmodel, data = Cau_50_72,
#                 stanvars = stanvars, prior = PRIORS,
#                 iter = 5000, chains = 4, cores = 4, seed = 2020,
#                 control = list(adapt_delta = 0.99, max_treedepth = 10))
#   save(brmout6, file="FvFm_Caulerpa_50_72.RData")
# } else {
#   load(file = "FvFm_Caulerpa_50_72.RData", verbose=T)
# }

#' Take a look at the samples.
# print(brmout1)
# print(brmout2)
# print(brmout3)
# print(brmout4)
# print(brmout5)
# print(brmout6)

alldata2 %>% pull(brmout)

### Examine HMC ###

## posterior_out1 = as.array(brmout1)
## posterior_out2 = as.array(brmout2)
## posterior_out3 = as.array(brmout3)
## posterior_out4 = as.array(brmout4)
## posterior_out5 = as.array(brmout4)
## posterior_out6 = as.array(brmout6)
## mcmc_rank_overlay(posterior_out1)
## mcmc_rank_overlay(posterior_out2)
## mcmc_rank_overlay(posterior_out3)
## mcmc_rank_overlay(posterior_out4)
## mcmc_rank_overlay(posterior_out5)
## mcmc_rank_overlay(posterior_out6)
expose_functions(brmout, vectorize = TRUE)
# expose_functions(brmout1, vectorize = TRUE)
# expose_functions(brmout2, vectorize = TRUE)
# expose_functions(brmout3, vectorize = TRUE)
# expose_functions(brmout4, vectorize = TRUE)
# expose_functions(brmout5, vectorize = TRUE)
# expose_functions(brmout6, vectorize = TRUE)

## y1 = brmout1$data %>% pull(Yield)
## yrep1 = posterior_predict(brmout1, nsamples = 500)
## ppc_dens_overlay(y1, yrep1)
## ppc_hist(y1, yrep1[1:8,])

## y2 = brmout2$data %>% pull(Yield)
## yrep2 = posterior_predict(brmout2, nsamples = 500)
## ppc_dens_overlay(y2, yrep2)
## ppc_hist(y2, yrep2[1:8,])

## y3 = brmout3$data %>% pull(Yield)
## yrep3 = posterior_predict(brmout3, nsamples = 500)
## ppc_dens_overlay(y3, yrep3)
## ppc_hist(y3, yrep3[1:8,])

## y4 = brmout4$data %>% pull(Yield)
## yrep4 = posterior_predict(brmout4, nsamples = 500)
## ppc_dens_overlay(y4, yrep4)
## ppc_hist(y4, yrep4[1:8,])

## y5 = brmout5$data %>% pull(Yield)
## yrep5 = posterior_predict(brmout5, nsamples = 500)
## ppc_dens_overlay(y5, yrep5)
## ppc_hist(y5, yrep5[1:8,])

# y6 = brmout6$data %>% pull(FvFm)
# yrep6 = posterior_predict(brmout6, nsamples = 500)
# ppc_dens_overlay(y6, yrep6)
# ppc_hist(y6, yrep6[1:8,])
# 
# xval = Cau_50_72$Temperature
# ppc_error_scatter_avg_vs_x(y6, yrep6, xval) +
#   geom_hline(yintercept = 0, linetype = "dashed")

alldata2 %>% 
  mutate(ppcdensoverlay = map(brmout, function(bout) {
    y = bout$data %>% pull(FvFm)
    yrep = posterior_predict(bout, nsamples = 500)
    ppc_dens_overlay(y, yrep)
    })) %>% pull(ppcdensoverlay) %>% 
  ggpubr::ggarrange(plotlist = .)


alldata2 %>% 
  mutate(fitted = map(brmout, function(bout) {
    pdata = bout$data %>% 
      expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
      add_predicted_draws(bout) %>%
      group_by(Temperature) %>%
      mean_hdci()
    
    
    ggplot() + 
      #  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
      #              data = pdata1, alpha = 0.2) +
      geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
                  data = pdata,
                  alpha = 0.2) +
      geom_line(aes(x = Temperature, y = .prediction), data = pdata) +
      geom_point(aes(x = Temperature, y = FvFm), 
                 data = bout$data) +
      labs(x="", y="") +
      scale_x_continuous(breaks = seq(8, 40, by = 4), limits=c(8,40)) +
      scale_y_continuous(breaks = seq(0, 1, length = 5), limits = c(0, 1)) +
      #  annotate("text", x=32, y=1, label="24 h", fontface="bold") +
      annotate("text", x=40, y=1, label="a", fontface="bold")
    
  })) %>% pull(fitted) %>% 
  ggpubr::ggarrange(plotlist = .)



#' # Build the plots
pdata1 = Cau_0_24 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_predicted_draws(brmout1) %>%
  group_by(Temperature) %>%
  mean_hdci()

edata1 = Cau_0_24 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_linpred_draws(brmout1) %>%
  group_by(Temperature) %>%
  mean_hdci()

pdata2 = Cau_0_48 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_predicted_draws(brmout2) %>%
  group_by(Temperature) %>%
  mean_hdci()

edata2 = Cau_0_48 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_linpred_draws(brmout2) %>%
  group_by(Temperature) %>%
  mean_hdci()

pdata3 = Cau_0_72 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_predicted_draws(brmout3) %>%
  group_by(Temperature) %>%
  mean_hdci()

edata3 = Cau_0_72 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_linpred_draws(brmout3) %>%
  group_by(Temperature) %>%
  mean_hdci()

pdata4 = Cau_50_24 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_predicted_draws(brmout4) %>%
  group_by(Temperature) %>%
  mean_hdci()

edata4 = Cau_50_24 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_linpred_draws(brmout4) %>%
  group_by(Temperature) %>%
  mean_hdci()

pdata5 = Cau_50_48 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_predicted_draws(brmout5) %>%
  group_by(Temperature) %>%
  mean_hdci()

edata5 = Cau_50_48 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_linpred_draws(brmout5) %>%
  group_by(Temperature) %>%
  mean_hdci()

pdata6 = Cau_50_72 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_predicted_draws(brmout6) %>%
  group_by(Temperature) %>%
  mean_hdci()

edata6 = Cau_50_72 %>%
  expand(Temperature = seq(min(Temperature), max(Temperature), by = 1)) %>%
  add_linpred_draws(brmout6) %>%
  group_by(Temperature) %>%
  mean_hdci()

### plots ###
h24_Cau_0 =
  ggplot() + 
  #  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
  #              data = pdata1, alpha = 0.2) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
              data = edata1,
              alpha = 0.2) +
  geom_line(aes(x = Temperature, y = .value), data = edata1) +
  geom_point(aes(x = Temperature, y = FvFm), data = Cau_0_24) +
  labs(x="", y="") +
  scale_x_continuous(breaks = seq(8, 40, by = 4), limits=c(8,40)) +
  scale_y_continuous(breaks = seq(0, 1, length = 5), limits = c(0, 1)) +
  #  annotate("text", x=32, y=1, label="24 h", fontface="bold") +
  annotate("text", x=40, y=1, label="a", fontface="bold") +
  theme_iris() 
h24_Cau_0

h48_Cau_0 =
  ggplot() + 
  #  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
  #              data = pdata2, alpha = 0.2) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
              data = edata2,
              alpha = 0.2) +
  geom_line(aes(x = Temperature, y = .value), data = edata2) +
  geom_point(aes(x = Temperature, y = FvFm), data = Cau_0_48) +
  labs(x="", y="") +
  scale_x_continuous(breaks = seq(8, 40, by = 4), limits=c(8,40)) +
  scale_y_continuous(breaks = seq(0, 1, length = 5), limits = c(0, 1)) +
  #  annotate("text", x=32, y=1, label="24 h", fontface="bold") +
  annotate("text", x=40, y=1, label="c", fontface="bold") +
  theme_iris() 
h48_Cau_0

h72_Cau_0 =
  ggplot() + 
  #  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
  #              data = pdata2, alpha = 0.2) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
              data = edata3,
              alpha = 0.2) +
  geom_line(aes(x = Temperature, y = .value), data = edata3) +
  geom_point(aes(x = Temperature, y = FvFm), data = Cau_0_72) +
  labs(x="", y="") +
  scale_x_continuous(breaks = seq(8, 40, by = 4), limits=c(8,40)) +
  scale_y_continuous(breaks = seq(0, 1, length = 5), limits = c(0, 1)) +
  #  annotate("text", x=32, y=1, label="24 h", fontface="bold") +
  annotate("text", x=40, y=1, label="e", fontface="bold") +
  theme_iris() 
h72_Cau_0

h24_Cau_50 =
  ggplot() + 
  #  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
  #              data = pdata2, alpha = 0.2) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
              data = edata4,
              alpha = 0.2) +
  geom_line(aes(x = Temperature, y = .value), data = edata4) +
  geom_point(aes(x = Temperature, y = FvFm), data = Cau_50_24) +
  labs(x="", y="") +
  scale_x_continuous(breaks = seq(8, 40, by = 4), limits=c(8,40)) +
  scale_y_continuous(breaks = seq(0, 1, length = 5), limits = c(0, 1)) +
  #  annotate("text", x=32, y=1, label="24 h", fontface="bold") +
  annotate("text", x=40, y=1, label="b", fontface="bold") +
  theme_iris() 
h24_Cau_50

h48_Cau_50 =
  ggplot() + 
  #  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
  #              data = pdata2, alpha = 0.2) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
              data = edata5,
              alpha = 0.2) +
  geom_line(aes(x = Temperature, y = .value), data = edata5) +
  geom_point(aes(x = Temperature, y = FvFm), data = Cau_50_48) +
  labs(x="", y="") +
  scale_x_continuous(breaks = seq(8, 40, by = 4), limits=c(8,40)) +
  scale_y_continuous(breaks = seq(0, 1, length = 5), limits = c(0, 1)) +
  #  annotate("text", x=32, y=1, label="24 h", fontface="bold") +
  annotate("text", x=40, y=1, label="d", fontface="bold") +
  theme_iris() 
h48_Cau_50

h72_Cau_50 =
  ggplot() + 
  #  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
  #              data = pdata2, alpha = 0.2) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature),
              data = edata6,
              alpha = 0.2) +
  geom_line(aes(x = Temperature, y = .value), data = edata6) +
  geom_point(aes(x = Temperature, y = FvFm), data = Cau_50_72) +
  labs(x="", y="") +
  scale_x_continuous(breaks = seq(8, 40, by = 4), limits=c(8,40)) +
  scale_y_continuous(breaks = seq(0, 1, length = 5), limits = c(0, 1)) +
  #  annotate("text", x=32, y=1, label="24 h", fontface="bold") +
  annotate("text", x=40, y=1, label="f", fontface="bold") +
  theme_iris() 
h72_Cau_50


xlabel1 = expression("Temperature"~(degree*C))
ylabel1 = expression("Maximum quantum yield,"~F[v]/F[m])


globalXlab = textGrob(xlabel1)
globalYlab = textGrob(ylabel1, rot=90)

h24_Cau_0 = ggplotGrob(h24_Cau_0)
h48_Cau_0 = ggplotGrob(h48_Cau_0)
h72_Cau_0 = ggplotGrob(h72_Cau_0)
h24_Cau_50 = ggplotGrob(h24_Cau_50)
h48_Cau_50 = ggplotGrob(h48_Cau_50)
h72_Cau_50 = ggplotGrob(h72_Cau_50)

g1 = rbind(h24_Cau_0, h48_Cau_0, h72_Cau_0,
           size = "first")
g2 = rbind(h24_Cau_50, h48_Cau_50, h72_Cau_50,
           size = "first")
g12 = cbind(g1, g2,
            size = "first")

ncols = length(gtable_width(g12)$arg1)
nrows = length(gtable_height(g12)$arg1)

g12 = gtable_add_rows(g12, heights = unit(1, "lines"), pos=nrows)
g12 = gtable_add_cols(g12, widths = unit(1, "lines"), pos=0)

ncols = length(gtable_width(g12)$arg1)
nrows = length(gtable_height(g12)$arg1)

g12 = gtable_add_grob(g12, globalXlab,
                      l = 4, r = ncols,
                      t = nrows, b = nrows)

g12 = gtable_add_grob(g12, globalYlab,
                      l = 1, r = 1,
                      t = 1, b = nrows-2)

# gtable_show_layout(g1)
grid.newpage()
grid.draw(g12)

#' Save plot
# Phycologia Figure single column and double column width  limits are 86 mm and 179 mm.
width = 179 / 25.4
height = 179 / 25.4

svg("Caulerpa_FvFm.svg",width = width, height = height)
grid.draw(g12)
dev.off()

pdf(file = "Caulerpa_FvFm.pdf",
    width = width,
    height = height)
grid.draw(g12)
dev.off()

##############################
SPO_24_PARS = 
  spread_draws(brmout1, b_zi_Intercept, b_PS_Intercept, b_HA_Intercept, b_ET_Intercept, b_KT_Intercept, b_zi_Temperature) %>%
  mutate(HD = b_HA_Intercept * b_ET_Intercept,
         CT = b_KT_Intercept - 273.15) %>%
  select(b_zi_Intercept:CT) %>%
  gather() %>%
  group_by(key) %>%
  mean_hdci(value) 

##
SPO_48_PARS = 
  spread_draws(brmout2, b_zi_Intercept, b_PS_Intercept, b_HA_Intercept, b_ET_Intercept, b_KT_Intercept, b_zi_Temperature) %>%
  mutate(HD = b_HA_Intercept * b_ET_Intercept,
         CT = b_KT_Intercept - 273.15) %>%
  select(b_zi_Intercept:CT) %>%
  gather() %>%
  group_by(key) %>%
  mean_hdci(value) 

##
SPO_72_PARS = 
  spread_draws(brmout3, b_zi_Intercept, b_PS_Intercept, b_HA_Intercept, b_ET_Intercept, b_KT_Intercept, b_zi_Temperature) %>%
  mutate(HD = b_HA_Intercept * b_ET_Intercept,
         CT = b_KT_Intercept - 273.15) %>%
  select(b_zi_Intercept:CT) %>%
  gather() %>%
  group_by(key) %>%
  mean_hdci(value) 

##
GAM_24_PARS = 
  spread_draws(brmout4, b_zi_Intercept, b_PS_Intercept, b_HA_Intercept, b_ET_Intercept, b_KT_Intercept, b_zi_Temperature) %>%
  mutate(HD = b_HA_Intercept * b_ET_Intercept,
         CT = b_KT_Intercept - 273.15) %>%
  select(b_zi_Intercept:CT) %>%
  gather() %>%
  group_by(key) %>%
  mean_hdci(value) 

##
GAM_48_PARS = 
  spread_draws(brmout5, b_zi_Intercept, b_PS_Intercept, b_HA_Intercept, b_ET_Intercept, b_KT_Intercept, b_zi_Temperature) %>%
  mutate(HD = b_HA_Intercept * b_ET_Intercept,
         CT = b_KT_Intercept - 273.15) %>%
  select(b_zi_Intercept:CT) %>%
  gather() %>%
  group_by(key) %>%
  mean_hdci(value) 

##
GAM_72_PARS = 
  spread_draws(brmout6, b_zi_Intercept, b_PS_Intercept, b_HA_Intercept, b_ET_Intercept, b_KT_Intercept, b_zi_Temperature) %>%
  mutate(HD = b_HA_Intercept * b_ET_Intercept,
         CT = b_KT_Intercept - 273.15) %>%
  select(b_zi_Intercept:CT) %>%
  gather() %>%
  group_by(key) %>%
  mean_hdci(value) 

#' Table of parameter values 
Parameters_h24_SPO = SPO_24_PARS
Parameters_h48_SPO = SPO_48_PARS
Parameters_h72_SPO = SPO_72_PARS
Parameters_h24_GAM = GAM_24_PARS
Parameters_h48_GAM = GAM_48_PARS
Parameters_h72_GAM = GAM_72_PARS

#+ results="asis"
print(xtable(Parameters_h24_SPO, caption="Caulerpa_24"), type="html", include.rownames=TRUE)
print(xtable(Parameters_h48_SPO, caption="Caulerpa_48"), type="html", include.rownames=TRUE)
print(xtable(Parameters_h72_SPO, caption="Caulerpa_72"), type="html", include.rownames=TRUE)

print(xtable(Parameters_h24_GAM, caption="Caulerpa_24"), type="html", include.rownames=TRUE)
print(xtable(Parameters_h48_GAM, caption="Caulerpa_48"), type="html", include.rownames=TRUE)
print(xtable(Parameters_h72_GAM, caption="Caulerpa_72"), type="html", include.rownames=TRUE)

#' ## Session information
sessionInfo()