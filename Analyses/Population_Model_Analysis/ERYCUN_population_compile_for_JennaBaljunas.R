####################################################################################
# Purpose: Compile vital rate estimates into IPM for investigation of impacts of seed predation on Eryngium cuneifolium
# Author: Joshua Fowler
# Date: May 23, 2025
####################################################################################
##### Set up #####

library(tidyverse)
library(readxl)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

library(popbio)

library(moments) # to calculate skew
library(patchwork) # to add ggplots together

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }


####################################################################################
###### sourcing in the vital rate model objects ####################################
####################################################################################

erycun.survival <- readRDS("erycun.survival.rds")
erycun.growth <- readRDS("erycun.growth.rds")
erycun.flw_status <- readRDS("erycun.flw_status.rds")
erycun.flw_stem <- readRDS("erycun.flw_stem.rds")
# erycun.flw_head <- readRDS("erycun.flw_head.rds")

# pulling out just the posterior means for each parameter
surv_par <- as.data.frame(erycun.survival) %>% summarize(across(everything(), ~mean(.x)))
grow_par <- as.data.frame(erycun.growth) %>% summarize(across(everything(), ~mean(.x)))
flw_par <- as.data.frame(erycun.flw_status) %>% summarize(across(everything(), ~mean(.x)))
flw_stem_par <- as.data.frame(erycun.flw_stem) %>% summarize(across(everything(), ~mean(.x)))

####################################################################################
###### sourcing in the MPM functions script ####################################
####################################################################################
source("Analyses/Population_Model_Analysis/MPM_functions_for_JennaBaljunas.R")

# calculating max size
# source(paste0(getwd(),"/Analyses/data_processing.R"))
ERYCUN_covariates <- read_csv("ERYCUN_covariates.csv")

ERYCUN <- ERYCUN_covariates %>% 
  mutate(log_size.t = log(ros_diameter.t),
         log_size.t1 = log(ros_diameter.t1),
         time_since_fire = time_since_fire_actual) %>% 
  filter(Site_tag != "royce_ranch" & bald != "200")%>% 
  filter(year.t >1989)



ERYCUN_surv.df <- ERYCUN %>% 
  filter(!is.na(surv.t1) & !is.na(log_size.t)) 
ERYCUN_flw_status.df <- ERYCUN %>% 
  filter(!is.na(flw_status.t) & !is.na(log_size.t)) 
ERYCUN_flw_stem.df <- ERYCUN %>% 
  filter(flw_status.t == 1 & !is.na(flw_stem.t) & flw_stem.t !=0 & !is.na(log_size.t))
ERYCUN_flw_head.df <- ERYCUN %>% 
  filter(flw_status.t == 1 & !is.na(flw_head.t) & flw_head.t>0 & !is.na(log_size.t))
ERYCUN_growth.df <- ERYCUN %>% 
  filter(!is.na(ros_diameter.t1) & !is.na(log_size.t)) 

size_bounds_df <- ERYCUN%>% 
  group_by() %>% 
  filter(!is.na(log_size.t1)) %>%
  summarise(max_size = max(log_size.t1),
            max_size_97 = quantile(log_size.t1,probs=0.975),
            max_size_99 = quantile(log_size.t1,probs=0.99),
            min_size = min(log_size.t1),
            min_size_97 = quantile(log_size.t1,probs=0.025),
            min_size_99 = quantile(log_size.t1,probs=0.01)
            )

# setting up models and params to feed in to MPM_functions

####################################################################################
###### setting up data averaging across balds for predictions across fire and elevation gradients ################################
####################################################################################

# setting up the covariates for each bald for prediction

fire_to2022 <- read_xlsx(path = "~/Dropbox/UofMiami/Balds2009_FireIntensityArea_Through2022.xlsx", sheet = "Balds2009_FireIntensityArea", guess_max = 1048576)# guess_max makes the function look deeper in the columns to assign type

fire_summary <- fire_to2022 %>% 
  group_by(Bald_U, Bald_) %>% 
  filter(INTENSITY != 0) %>% 
  summarize(last_fire = max(as.numeric(Year), na.rm = TRUE),
            time_since_fire = 2023-max(as.numeric(Year), na.rm = TRUE),
            fire_list = toString(unique(as.numeric(Year))),
            fire_frequency = length(unique(as.numeric(Year)))) %>% 
  ungroup() %>% 
  add_row(Bald_U = "1S", Bald_ = 1, last_fire = NA, time_since_fire = NA, fire_list = NA, fire_frequency = 0)

# Now getting the relative elevation from the old fire history file

elev.df <- read_xlsx(path = "~/Dropbox/UofMiami/Experiment Set up/firehistory_thru2018.xlsx", sheet = "Rx_Freq", guess_max = 1048576) %>% # guess_max makes the function look deeper in the columns to assign type
  rename(rel_elev = rel.eve) %>% 
  mutate(
    bald = case_when(bald == "01S" ~ "1S",
                     bald == "01N" ~ "1N",
                     bald == "02" ~ "2",
                     bald == "05E" ~ "5E",
                     bald == "07N" ~ "7N",
                     bald == "35N" ~ "35",
                     bald == "65E" ~ "65",
                     TRUE ~ bald)) %>% 
  select(bald, name, rel_elev)




preddata_fire <- as_tibble(expand.grid(bald = NA, 
                             year.t1 = NA,
                             time_since_fire = seq(from = min(fire_summary$time_since_fire, na.rm = T), to = max(fire_summary$time_since_fire, na.rm = T), length = 6),
                             rel_elev = mean(elev.df$rel_elev, na.rm = T),
                             total_seeds=1) )


preddata_elev <- as_tibble(expand.grid(bald = NA, 
                             year.t1 = NA,
                             time_since_fire = mean(fire_summary$time_since_fire, na.rm = T),
                             rel_elev = seq(from = min(elev.df$rel_elev, na.rm = T), to = max(elev.df$rel_elev, na.rm = T), length = 10),
                             total_seeds=1) )



####################################################################################
###### Calculating population summaries ####################################
####################################################################################
# Taking seed dynamics from Hindle et al. "The two fertility scenarios which produced dynamics best fitting to the observed population dynamics were low first year germination (c_f=0), low germination from the seed bank (c_b=0.005), and low seed mortality (d=0.3) and low first year germination (c_f=0), high germination from the seed bank (c_b=0.04), and relatively high seed mortality (d=0.7; Fig. A5). "
models <- make_mods(grow = erycun.growth, surv = erycun.survival, flw = erycun.flw_status, fert = erycun.flw_stem, 
                    seeds_per_stem = 183, seed_mortality = .3, seed_germ1 = 0, seed_germ2 = .005,
                    seedling_surv = erycun.survival, seedling_size = erycun.growth)
# params <- make_params(bald.rfx = F, bald = NA, year.rfx = F, year = NA, 
#                       post_draws = 1,
#                       iter = 1,
#                       preddata = preddata, 
#                       microbe = 0,  # 0 is alive, and 1 is sterile becuase we start with the microbes in the model, but then could choose to turn off the microbes
#                       surv_par = surv_par,
#                       grow_par = grow_par, 
#                       flw_par = flw_par,
#                       fert_par = flw_stem_par,
#                       size_bounds = size_bounds_df)

# gxy(0,0,models, params)
# 
# plot(fx(1:10, models, params))
# 
# 


# we can calculate lambda, but we might also consider later looking at effects of microbes on quantities like the stable stage distribution etc.
nfire <- nrow(preddata_fire)

# Currently, I pulled the parameters as the posterior mean, but we could do this instead by taking say 500 draws of each posterior to incorporate uncertainty.
ndraws <- 1
post_draws <- 1#sample.int(7500,size=ndraws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.


IPM_fire_iter <- list(NA)
IPM_fire <- list(NA)
# lambda_fire <- array(NA, dim = c(ndraws, nfire))
# the model is now set up to use the p

for(f in 1:nfire){
  preddata <- preddata_fire[f,]
  for(i in 1:ndraws){
    IPM_fire_iter[[i]] <- bigmatrix(params = make_params(bald.rfx = F, bald = NA, year.rfx = F, year = NA, 
                                                                     post_draws = post_draws,
                                                                     iter = i,
                                                                     preddata = preddata, 
                                                                     microbe = 0,  # 0 is alive, and 1 is sterile becuase we start with the microbes in the model, but then could choose to turn off the microbes
                                                                     surv_par = surv_par,
                                                                     grow_par = grow_par, 
                                                                     flw_par = flw_par,
                                                                     fert_par = flw_stem_par,
                                                                     size_bounds = size_bounds_df), 
                                                models = models, matdim = 25, extension = 2)$MPMmat
  }
  IPM_fire[[f]] <- IPM_fire_iter
  # print(paste("fire", f))
}



IPM_fire <- unlist(IPM_fire,recursive=F)
saveRDS(IPM_fire, "IPM_fire.Rds")

# example of calculating lambda from the first matrix, which here represents the lowest level of time since fire
popbio::lambda(IPM_fire[[1]])
plot(unlist(lapply(IPM_fire, FUN = popbio::lambda)))

lambda_fire_df <- data.frame(lambda = unlist(lapply(IPM_fire, FUN = popbio::lambda)),
                          time_since_fire = unique(preddata_fire$time_since_fire))


# saveRDS(lambda, "lambda_fire.Rds")
# lambda_fire <- readRDS("lambda_fire.Rds")
# dimnames(lambda_fire) <- list(Iter = paste0("iter",1:ndraws), 
#                          Fire = unique(preddata_fire$time_since_fire))
# lambda


# write_csv(lambda_fire_df, "prelim_fire_lambdas.csv")

# lambda_summary <- lambda_fire_df %>% 
#   filter(!is.na(lambda)) %>% 
#   group_by(Bald, Microbe) %>% 
#   summarize(lambda_mean = mean(lambda),
#             lambda_97.5 = quantile(lambda, 0.975),
#             lambda_02.5 = quantile(lambda, 0.025))

ggplot()+
  geom_hline(aes(yintercept = 1), linetype = "dashed")+
  geom_jitter(data = lambda_fire_df, aes(x = time_since_fire, y = lambda), width = .05, alpha = .4)+
  # geom_point(data = lambda_summary, aes(x = Bald, y = lambda_mean, color = Microbe, group = "microbe"), size = 4)+
  # geom_linerange(data = lambda_summary, aes(x = Bald, ymin = lambda_02.5, ymax = lambda_97.5, color = Microbe), lwd = 1)+
  # ylim(0,10)+
  theme_classic()



# plot(y,sx(y,models,params),xlab="Size",type="l",
#      ylab="Survival Probability",lwd=12)
# points(y,apply(bigmatrix(params = params, models = models, matdim = 100, extension = 2)$Tmat,2,sum),col="red",lwd=3,cex=.1,pch=19)






# now evaluating across elevation
nelev <- nrow(preddata_elev)

IPM_elev_iter <- list()
IPM_elev <- list()

# The model is set up to calculate the posterior means, but we could re-write to incorporate the uncertainty in the vital rate estimates
for(e in 1:nelev){
  preddata <- preddata_elev[e,]
  for(i in 1:ndraws){
    IPM_elev_iter[[i]] <- bigmatrix(params = make_params(bald.rfx = F, bald = NA, year.rfx = F, year = NA, 
                                                    post_draws = post_draws,
                                                    iter = i,
                                                    preddata = preddata, 
                                                    microbe = 0,  # 0 is alive, and 1 is sterile becuase we start with the microbes in the model, but then could choose to turn off the microbes
                                                    surv_par = surv_par,
                                                    grow_par = grow_par, 
                                                    flw_par = flw_par,
                                                    fert_par = flw_stem_par,
                                                    size_bounds = size_bounds_df), 
                               models = models, matdim = 25, extension = 2)$MPMmat
  }
  IPM_elev[[e]] <- IPM_elev_iter
  # print(paste("fire", f))
}

# saveRDS(IPM_elev, "IPM_elev.Rds")

IPM_elev <- unlist(IPM_elev,recursive=F)
saveRDS(IPM_elev, "IPM_elev.Rds")


# example of calculating lambda from the first matrix, which here represents the lowest level of relative elevation. Plotting them across the gradient shows a bit more squiglyness than I expected so I'm not sure what is up with that.
popbio::lambda(IPM_elev[[1]])
plot(unlist(lapply(IPM_elev, FUN = popbio::lambda)), x = preddata_elev$rel_elev)
