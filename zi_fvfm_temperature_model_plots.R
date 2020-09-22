library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)

load("Caulerpa_brms_samples_20200917.RData")
brmout = readRDS("zimodel.rds")
expose_functions(brmout, vectorize = T)


# Diagnostic plots -------------------------------------------------------------
diagnostic_plots = outdata %>% 
  mutate(summary = map(brmout, summary)) %>%
  mutate(ppc_dens_overlay = map(brmout, function(bfit) {
    y = bfit$data %>% pull(FvFm) 
    yrep = posterior_predict(bfit, nsamples = 100)
    ppc_dens_overlay(y, yrep)
  })) %>% 
  mutate(mcmc_rank_overlay = map(brmout, function(bfit){
    posterior_out = as.array(bfit)
    mcmc_rank_overlay(posterior_out)      
  }))

diagnostic_plots %>% 
  pull(ppc_dens_overlay) %>% 
  ggpubr::ggarrange(plotlist = .)


diagnostic_plots %>% 
  pull(mcmc_rank_overlay) %>% 
  ggpubr::ggarrange(plotlist = .)

posterior_summary = toplot = outdata %>% 
  mutate(predict = map(brmout, function(bfit) {
    bfit$data %>% 
      expand(Temperature = seq(min(Temperature), 
                               max(Temperature), length = 18)) %>% 
      distinct() %>% 
      mutate(Temperature2 = Temperature^2) %>% 
      add_predicted_draws(bfit) %>% 
      group_by(Temperature) %>% 
      mean_hdci(.prediction)
  })) %>% 
  mutate(expect = map(brmout, function(bfit) {
    bfit$data %>% 
      expand(Temperature = seq(min(Temperature), 
                               max(Temperature), length = 30)) %>% 
      distinct() %>% 
      mutate(Temperature2 = Temperature^2) %>% 
      add_linpred_draws(bfit) %>% 
      group_by(Temperature) %>% 
      mean_hdci(.value)
  })) %>% 
  mutate(zi = map(brmout, function(bfit) {
    bfit$data %>% 
      expand(Temperature = seq(min(Temperature), 
                               max(Temperature), length = 30)) %>% 
      distinct() %>% 
      mutate(Temperature2 = Temperature^2) %>% 
      add_linpred_draws(bfit, dpar = TRUE) %>% 
      group_by(Temperature) %>% 
      mean_hdci(zi)
  }))


xlabF = expression("Temperature" ~ ({}^degree*C ))
ylabF = "Maximum Quantum Yield (Fv/Fm)" 
ltlabF = "Process"
flabF = "Process"
ggplot() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature,
                  fill = "Prediction"),
              alpha = 0.2,
              data = posterior_summary %>% unnest(predict)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature,
                  fill = "Expectation"),
              alpha = 0.4,
              data = posterior_summary %>% unnest(expect)) +
  geom_line(aes(x = Temperature, y = .value,
                linetype = "Expectation"),
            data = posterior_summary %>% unnest(expect)) +
  geom_point(aes(x = Temperature, y = FvFm),
             data = toplot %>% unnest(data),
             position = position_jitter(0.2)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature,
                  fill = "Zero-inflation rate"), 
              alpha = 0.4,
              data = posterior_summary %>% unnest(zi)) +
  geom_line(aes(x = Temperature, y = zi,
                linetype = "Zero-inflation rate"), 
            data = posterior_summary %>% unnest(zi)) +
  scale_x_continuous(xlabF) + 
  scale_y_continuous(ylabF) +
  scale_linetype_discrete(ltlabF) +
  scale_fill_discrete(flabF) +
  facet_grid(rows = vars(Hour),
             cols = vars(Light)) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_blank())

ggsave(filename = "fvfm_plot.png",
       width =  3*80,
       height = 4*80,
       units = "mm")



outdata %>% 
  mutate(parameters = map(brmout, function(df) {
    df %>% 
      spread_draws(phi, `b_*.*`, regex = T) %>% 
      mutate(HD = b_HA_Intercept * b_ET_Intercept,
             CT = b_KT_Intercept - 273.15) %>% 
      dplyr::select(c(starts_with("b"), phi)) %>% 
     gather() %>%
         group_by(key) %>%
         mean_hdci(value) %>% 
      mutate(key = str_remove(key, "_Intercept")) %>% 
      mutate(key = str_remove(key, "b_"))
    })) %>% 
  unnest(parameters) %>% 
  filter(!str_detect(key, "zi")) %>% 
  ggplot() +
  geom_pointrange(aes(x = Hour, 
                      y = value, 
                      ymin = .lower,
                      ymax = .upper,
                      color = Hour),
           position = position_dodge(0.3)) +
  facet_wrap(vars(Light,key),
             nrow = 2,
             scales = "free")

outdata %>% 
  mutate(parameters = map(brmout, function(df) {
    df %>% 
      spread_draws(phi, `b_*.*`, regex = T) %>% 
      mutate(KTOPT = -b_zi_Temperature/b_zi_Temperature2) %>% 
      dplyr::select(c(starts_with("b_zi"), KTOPT)) %>% 
      gather() %>%
      group_by(key) %>%
      mean_hdci(value) 
  })) %>% 
  unnest(parameters) %>% print() %>% 
  ggplot() +
  geom_pointrange(aes(x = Hour, 
                      y = value, 
                      ymin = .lower,
                      ymax = .upper,
                      color = Hour),
                  position = position_dodge(0.3)) +
  facet_wrap(vars(Light,key),
             nrow = 2,
             scales = "free")



x = seq(0, 40, by = 0.2)
y = 31 + -3.1*x + 0.06 * x *x
plot(x,exp(y)/(1+exp(y)))


SPO_24_PARS = 
  spread_draws(brmout1, b_zi_Intercept, b_PS_Intercept, b_HA_Intercept, b_ET_Intercept, b_KT_Intercept, b_zi_Temperature) %>%
  mutate(HD = b_HA_Intercept * b_ET_Intercept,
         CT = b_KT_Intercept - 273.15) %>%
  select(b_zi_Intercept:CT) %>%
  gather() %>%
  group_by(key) %>%
  mean_hdci(value) 

