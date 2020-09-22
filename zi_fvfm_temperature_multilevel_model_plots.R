library(tidyverse)
library(rstan)
library(brms)
library(tidybayes)
library(bayesplot)

load("Caulerpa_brms_samples_20200917_0mumol.RData", verbose = T)
load("Caulerpa_brms_samples_20200917_50mumol.RData", verbose = T)
expose_functions(outdata00, vectorize = T)

# Diagnostic plots -------------------------------------------------------------
y = outdata00$data %>% pull(FvFm) 
yrep = posterior_predict(outdata00, nsamples = 100)
ppc_dens_overlay(y, yrep)

y = outdata50$data %>% pull(FvFm) 
yrep = posterior_predict(outdata50, nsamples = 100)
ppc_dens_overlay(y, yrep)

rstan::check_hmc_diagnostics(outdata00$fit)
rstan::check_hmc_diagnostics(outdata50$fit)

# posterior_out = as.array(outdata00)
# mcmc_rank_overlay(posterior_out)      
# 
# posterior_out = as.array(outdata00)
# mcmc_rank_overlay(posterior_out)      


makeplot = function(brmout, ntemperature = 18*2) {
  
  
  posterior_predict = 
    brmout$data %>% 
    expand(Temperature = seq(min(Temperature), 
                             max(Temperature), length = ntemperature),
           Hour) %>% 
    distinct() %>% 
    mutate(Temperature2 = Temperature^2) %>% 
    add_predicted_draws(brmout) %>% 
    group_by(Hour, Temperature) %>% 
    mean_hdci(.prediction)
  
  posterior_expect = brmout$data %>% 
    expand(Temperature = seq(min(Temperature), 
                             max(Temperature), length = ntemperature),
           Hour) %>% 
    distinct() %>% 
    mutate(Temperature2 = Temperature^2) %>% 
    add_linpred_draws(brmout) %>% 
    group_by(Hour, Temperature) %>% 
    mean_hdci(.value)
  
  zi_expect = brmout$data %>% 
    expand(Temperature = seq(min(Temperature), 
                             max(Temperature), length = ntemperature),
           Hour) %>% 
    distinct() %>% 
    mutate(Temperature2 = Temperature^2) %>% 
    add_linpred_draws(brmout, dpar = TRUE) %>% 
    group_by(Hour,Temperature) %>% 
    mean_hdci(zi)
  
  xlabF = expression("Temperature" ~ ({}^degree*C ))
  ylabF = "Maximum Quantum Yield (Fv/Fm)" 
  ltlabF = "Process"
  flabF = "Process"
  ggplot() + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature,
                    fill = "Prediction"),
                alpha = 0.2,
                data = posterior_predict) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature,
                    fill = "Expectation"),
                alpha = 0.4,
                data = posterior_expect) +
    geom_line(aes(x = Temperature, y = .value,
                  linetype = "Expectation"),
              data = posterior_expect) +
    geom_point(aes(x = Temperature, y = FvFm),
               data = brmout$data,
               position = position_jitter(0.2)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, x = Temperature,
                    fill = "Zero-inflation rate"), 
                alpha = 0.4,
                data = zi_expect) +
    geom_line(aes(x = Temperature, y = zi,
                  linetype = "Zero-inflation rate"), 
              data = zi_expect) +
    scale_x_continuous(xlabF) + 
    scale_y_continuous(ylabF) +
    scale_linetype_discrete(ltlabF) +
    scale_fill_discrete(flabF) +
    facet_grid(rows = vars(Hour)) +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.title = element_blank())
}

p1 = makeplot(brmout = outdata00)
p2 = makeplot(brmout = outdata50)

outdata00
outdata50

p12 = ggpubr::ggarrange(p1,p2, ncol = 2)
p12

ggsave(filename = "fvfm_plot.png",
       plot = p12,
       width =  3*80,
       height = 4*80,
       units = "mm")



# outdata %>% 
#   mutate(parameters = map(brmout, function(df) {
#     df %>% 
#       spread_draws(phi, `b_*.*`, regex = T) %>% 
#       mutate(HD = b_HA_Intercept * b_ET_Intercept,
#              CT = b_KT_Intercept - 273.15) %>% 
#       dplyr::select(c(starts_with("b"), phi)) %>% 
#      gather() %>%
#          group_by(key) %>%
#          mean_hdci(value) %>% 
#       mutate(key = str_remove(key, "_Intercept")) %>% 
#       mutate(key = str_remove(key, "b_"))
#     })) %>% 
#   unnest(parameters) %>% 
#   filter(!str_detect(key, "zi")) %>% 
#   ggplot() +
#   geom_pointrange(aes(x = Hour, 
#                       y = value, 
#                       ymin = .lower,
#                       ymax = .upper,
#                       color = Hour),
#            position = position_dodge(0.3)) +
#   facet_wrap(vars(Light,key),
#              nrow = 2,
#              scales = "free")
# 
# outdata %>% 
#   mutate(parameters = map(brmout, function(df) {
#     df %>% 
#       spread_draws(phi, `b_*.*`, regex = T) %>% 
#       mutate(KTOPT = -b_zi_Temperature/b_zi_Temperature2) %>% 
#       dplyr::select(c(starts_with("b_zi"), KTOPT)) %>% 
#       gather() %>%
#       group_by(key) %>%
#       mean_hdci(value) 
#   })) %>% 
#   unnest(parameters) %>% print() %>% 
#   ggplot() +
#   geom_pointrange(aes(x = Hour, 
#                       y = value, 
#                       ymin = .lower,
#                       ymax = .upper,
#                       color = Hour),
#                   position = position_dodge(0.3)) +
#   facet_wrap(vars(Light,key),
#              nrow = 2,
#              scales = "free")
# 
# 
# 
# x = seq(0, 40, by = 0.2)
# y = 31 + -3.1*x + 0.06 * x *x
# plot(x,exp(y)/(1+exp(y)))
# 
# 
# SPO_24_PARS = 
#   spread_draws(brmout1, b_zi_Intercept, b_PS_Intercept, b_HA_Intercept, b_ET_Intercept, b_KT_Intercept, b_zi_Temperature) %>%
#   mutate(HD = b_HA_Intercept * b_ET_Intercept,
#          CT = b_KT_Intercept - 273.15) %>%
#   select(b_zi_Intercept:CT) %>%
#   gather() %>%
#   group_by(key) %>%
#   mean_hdci(value) 
# 
