#==============================================================================#
#	 █████╗ ███╗   ██╗ █████╗ ██╗  ██╗   ██╗███████╗██╗███████╗	       #
#	██╔══██╗████╗  ██║██╔══██╗██║  ╚██╗ ██╔╝██╔════╝██║██╔════╝	       #
#	███████║██╔██╗ ██║███████║██║   ╚████╔╝ ███████╗██║███████╗	       #
#	██╔══██║██║╚██╗██║██╔══██║██║    ╚██╔╝  ╚════██║██║╚════██║	       #
#	██║  ██║██║ ╚████║██║  ██║███████╗██║   ███████║██║███████║	       #
#	╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝   ╚══════╝╚═╝╚══════╝	       #
#==============================================================================#

#==============================================================================#
#	Loading Packages						       #
#==============================================================================#	

library(purrr)
library(dplyr)
library(corrplot)
library(Hmisc)
library(texreg)
library(margins)
library(ggplot2)
library(RColorBrewer)
library(lemon)
library(ggrepel)
library(ggeffects)
library(splines)
library(conflicted)
library(MASS)

#==============================================================================#
#	Loading Data							       #
#==============================================================================#	

load("../Data/GeoISDControls.Rdata")

#==============================================================================#
#	Creating New Variables						       #
#==============================================================================#

prio_grid_isd <- prio_grid_isd %>% mutate(
				sqrtDeaths = sqrt(deaths),
				sqrtState_based = sqrt(state_based),
				sqrtPpp = sqrt(gcp_ppp),
				sqrtSpAll = sqrt(sp_os_i_sum),
				sqrtSpAny = sqrt(sp_os_i_sum_any),
				sqrtCapdist = sqrt(capdist),
				sqrtBDist = sqrt(bdist3),
				sqrtPopd = sqrt(popd),
				logDeaths = log(deaths + 1),
				logState_based = log(state_based + 1),
				logPpp = log(gcp_ppp + 1),
				logSpAll = log(sp_os_i_sum + 1),
				logSpAny = log(sp_os_i_sum_any + 1),
				logCapdist = log(capdist + 1),
				logBDist = log(bdist3 + 1),
				logPopd = log(popd + 1),
				sqrtorg3 = sqrt(org3),
				sqrtNonstate = sqrt(non_state),
				SpAllXCapDist = logSpAll*capdist,
				SpAnyXCapDist = logSpAny*capdist)

#==============================================================================#
#	Functions				    	             	       #
#==============================================================================#

GeoISDanalysis <- function(dvs, ivs, controls, data, test_label){
  models_out_j <- NULL
  models_out <- NULL
  margins_out <- NULL
  for(i in 1:length(dvs)) {
    dv <- dvs[i] 
    for(j in 1:length(ivs)) {
      m1 <- lm(as.formula(paste(dv, '~', ivs[j], controls)),
               data=data)
      m2 <- lm(as.formula(paste(dv, '~', ivs[j], controls, extended_controls)),
               data=data)
      m3 <- lm(as.formula(paste(dv, '~', ivs[j], controls, extended_controls, 
				all_controls)), data=data)
      iv_margins <- ivs[j] %>% strsplit(., "\\+") %>% unlist()
      m1_m <- summary(margins(m1, variables=iv_margins, change = "sd")) %>% 
        mutate(test='Baseline model', dep_var=dv)
      m2_m <- summary(margins(m2, variables=iv_margins, change = "sd")) %>% 
        mutate(test='Extended controls', dep_var=dv)
      m3_m <- summary(margins(m3, variables=iv_margins, change = "sd")) %>% 
        mutate(test='Full model', dep_var=dv)
      margins_out <- bind_rows(margins_out, m1_m, m2_m, m3_m)
      models_out_j[[j]] <- list(m1, m2, m3)
    } 
    models_out[[i]] <- models_out_j %>% purrr::flatten()
  }
  models_out <- models_out
  final <- list('margins'=margins_out, 'models'=models_out)
  return(final)
}

GeoISD_NBanalysis <- function(dvs, ivs, controls, data, test_label){
  models_out_j <- NULL
  models_out <- NULL
  margins_out <- NULL
  for(i in 1:length(dvs)) {
    dv <- dvs[i] 
    for(j in 1:length(ivs)) {
#      m1 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls)),
#               data=data)
#      m2 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls, extended_controls)),
#               data=data)
      m3 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls, extended_controls, 
				all_controls)), data=data)
      iv_margins <- ivs[j] %>% strsplit(., "\\+") %>% unlist()
#      m1_m <- summary(margins(m1, variables=iv_margins, change = "sd")) %>% 
#        mutate(test='Baseline model', dep_var=dv)
#      m2_m <- summary(margins(m2, variables=iv_margins, change = "sd")) %>% 
#        mutate(test='Extended controls', dep_var=dv)
      m3_m <- summary(margins(m3, variables=iv_margins, change = "sd")) %>% 
        mutate(test='Full model', dep_var=dv)
      margins_out <- bind_rows(margins_out, m3_m)
      models_out_j[[j]] <- list(m3)
    } 
    models_out[[i]] <- models_out_j %>% purrr::flatten()
  }
  models_out <- models_out
  final <- list('margins'=margins_out, 'models'=models_out)
  return(final)
}

#==============================================================================#
#	Descriptive Statistics                                         	       #
#==============================================================================#

# Correlation Matrix Plot

corrdata <- prio_grid_isd %>% dplyr::select(bdist3, capdist, excluded, temp_sd,
					    temp, prec_sd, prec_gpcc, barren_gc,
					    forest_gc, mountains_mean, water_gc,
					    sp_os_sum_any, sp_os_sum,
					    sp_os_i_sum_any, sp_os_i_sum,
					    distcoast, popd, state_based,
					    non_state, org3, deaths, gcp_mer,
					    gcp_ppp)

corrdata_mat <- cor(corrdata, method = c("spearman"), 
                     use="pairwise.complete.obs")

corrdata_mat <- rcorr(as.matrix(corrdata))

# TODO: Print to file
corrplot(corrdata_mat$r, method='color', diag=F, addCoef.col = "black", 
         p.mat = corrdata_mat$P, insig = "blank", tl.col = "black",
         tl.srt = 45)

# Summary Statistics

# TODO: Print to file
summary(corrdata)

#==============================================================================#
#	Variables							       #
#==============================================================================#

controls <- c('+ logCapdist + mountains_mean + water_gc + distcoast')

sqrt_controls <- c('+ sqrtCapdist + mountains_mean + water_gc + distcoast')

extended_controls <- c('+ temp_sd + temp + prec_sd + prec_gpcc + barren_gc +
		       forest_gc')

all_controls <- c('+ logPopd + logBDist')

sqrt_all_controls <- c('+ sqrtPopd + sqrtBDist')

control_names <-c('Baseline', 'Extetended Controls', 'Full Model',
		 'Baseline', 'Extetended Controls', 'Full Model')

dvs <- c('logDeaths', 'logState_based', 'non_state', 'org3', 'logPpp',
	 'nightlights')

nb_dvs <- c('deaths', 'state_based', 'non_state', 'sqrtorg3', 'logPppp', 'nightlights')

captions <- c('Deaths (logged)', 'State based conflict events (logged)',
	      'Non state conflict events', 'Communal violence events', 'PPP
	      (logged)')

ivs <- c('logSpAll', 'logSpAny')

srqt_ivs <- c('sqrtSpAll', 'sqrtSpAny')

# interactions <- c('SpAllXCapDist', 'SpAnyXCapDist')

#==============================================================================#
#	Analysis							       #
#==============================================================================#

linear_models <- GeoISDanalysis(dvs=dvs, ivs=ivs, controls=controls,
				 data=prio_grid_isd, test_label='Linear
				 Models')

#interaction_models <- GeoISDanalysis(dvs=dvs, ivs=interactions, controls=controls,
#				 data=prio_grid_isd, test_label='Interaction
#				 Models')



for (i in 1:length(dvs)) { 
  name <- dvs[i]
  filename <- paste("../Output/",name,".tex",sep="") 
  texreg(linear_models$models[[i]],
         file = filename,
         custom.model.names = control_names,
         #custom.coef.map = coefs,
         stars = c(0.001, 0.01, 0.05, 0.1), 
         sideways = T, use.packages = F, scalebox = 1,
         #custom.note = "Reference region is Eastern Europe and Central Asia.",
         caption = captions[i],
	 label = name,
         table = T)
}

full_margins <- linear_models[[1]] %>% 
  mutate(factor_label = case_when(
    factor == 'logSpAll' ~ 'State Preseance (all)(logged)',
    factor == 'logSpAny' ~ 'State Preseance (any)(logged)'),
	 facet_label = case_when(
	dep_var == 'logDeaths' ~ 'Total fatalities (logged)',
	dep_var == 'logState_based' ~ 'State based conflict events (logged)',
	dep_var == 'non_state' ~ 'Non state conflict events',
	dep_var == 'org3' ~ 'Communal violence events',
	dep_var == 'logPpp' ~ 'PPP (logged)',
	dep_var == 'nightlights' ~ 'Mean nightlights'))  %>% 
	rename(Test=test)


margins_main_plot <- ggplot(data = full_margins,
  aes(x = factor_label, y = AME, ymin = lower, ymax = upper, colour=Test)) +
  facet_wrap(~ facet_label) +
  xlab('') +
  labs(color = 'Models') +
  scale_color_manual(values=(brewer.pal(n = 4, name = "Dark2"))) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.1) +
  coord_flip() +
  theme(panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) +
  geom_hline(yintercept=0, linetype="dotted") +
  guides(alpha="none") +
  labs(caption = "Points are average marginal effects of a 1 SD increase in the 
 independent variable with 95% confidence intervals.")

# TODO: Print to file
margins_main_plot

# Marginal effects (interaction models)
fit <- lm(logDeaths ~ logSpAll * logCapdist + mountains_mean +
	  water_gc + distcoast + logPopd + logBDist + temp_sd + temp + prec_sd +
	  prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

predplot <- ggpredict(fit, terms = c("logSpAll", "logCapdist"))

ggplot(predplot, aes(x, exp(predicted))) +
  geom_line() +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), alpha = .1)

#=============================================================================#
# Negative binomial models

NB_models <- GeoISD_NBanalysis(dvs=nb_dvs, ivs=ivs, controls=controls,
				 data=prio_grid_isd, test_label='Negative
				 binomial Models')

NB_deaths <- glm.nb(deaths ~ logSpAny + logCapdist + mountains_mean +
 			water_gc + distcoast + logPopd + logBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

NB_state_based <- glm.nb(state_based ~ logSpAny + logCapdist + mountains_mean +
 			water_gc + distcoast + logPopd + logBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

NB_non_state <- glm.nb(non_state ~ logSpAny + logCapdist + mountains_mean +
 			water_gc + distcoast + logPopd + logBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

NB_org3 <- glm.nb(sqrtorg3 ~ logSpAny + logCapdist + mountains_mean +
 			water_gc + distcoast + logPopd + logBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

NB_ppp <- glm.nb(logPpp ~ logSpAny + logCapdist + mountains_mean +
 			water_gc + distcoast + logPopd + logBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

NB_nightlights <- glm.nb(nightlights ~ logSpAny + logCapdist + mountains_mean + 
 			water_gc + distcoast + logPopd + logBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

nb_full_margins <- NB_models[[1]] %>% 
  mutate(factor_label = case_when(
    factor == 'logSpAll' ~ 'State Preseance (all)(logged)',
    factor == 'logSpAny' ~ 'State Preseance (any)(logged)'),
	facet_label = case_when(
	dep_var == 'deaths' ~ 'Total fatalities',
	dep_var == 'state_based' ~ 'State based conflict events',
	dep_var == 'non_state' ~ 'Non state conflict events',
	dep_var == 'sqrtorg3' ~ 'Communal violence events (square root
	transformed)')) %>% 
	rename(Test=test)


margins_nb_plot <- ggplot(data = nb_full_margins,
  aes(x = factor_label, y = AME, ymin = lower, ymax = upper, colour=Test)) +
  facet_wrap(~ facet_label) +
  xlab('') +
  labs(color = 'Models') +
  scale_color_manual(values=(brewer.pal(n = 4, name = "Dark2"))) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.1) +
  coord_flip() +
  theme(panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) +
  geom_hline(yintercept=0, linetype="dotted") +
  guides(alpha="none") +
  labs(caption = "Points are average marginal effects of a 1 SD increase in the 
 independent variable with 95% confidence intervals.")

margins_nb_plot

#=============================================================================#
# Square root negative binomial models


NB_sqrt_deaths_min <- glm.nb(sqrtDeaths ~ sqrtSpAny + logCapdist +
			     mountains_mean + distcoast + barren_gc, 
		data = prio_grid_isd)

summary(NB_sqrt_deaths_min )

NB_sqrt_deaths <- glm.nb(sqrtDeaths ~ sqrtSpAny + logCapdist + mountains_mean +
 			water_gc + distcoast + sqrtPopd + sqrtBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

summary(NB_sqrt_deaths )

NB_sqrt_state_based_min <- glm.nb(sqrtState_based ~ sqrtSpAny + sqrtCapdist + mountains_mean +
 			distcoast + barren_gc, 
		data = prio_grid_isd)

summary(NB_sqrt_state_based_min )

NB_sqrt_state_based <- glm.nb(sqrtState_based ~ sqrtSpAny + logCapdist + mountains_mean +
 			water_gc + distcoast + sqrtPopd + sqrtBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

summary(NB_sqrt_state_based )

NB_org3 <- glm.nb(sqrtorg3 ~ sqrtSpAny + logCapdist + 
		  mountains_mean + distcoast + barren_gc,
		     data = prio_grid_isd)

summary(NB_org3)

NB_org3_all <- glm.nb(sqrtorg3 ~ sqrtSpAny + logCapdist + mountains_mean + barren_gc +
 			water_gc + distcoast + logPopd + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

summary(NB_org3_all)

#=============================================================================#
# Square root negative binomial interaction models 

NB_sqrt_deaths_min_inter <- glm.nb(sqrtDeaths ~ sqrtSpAny * logCapdist + mountains_mean +
 			distcoast + barren_gc, 
		data = prio_grid_isd)

summary(NB_sqrt_deaths_min_inter )

NB_sqrt_deaths_inter <- glm.nb(sqrtDeaths ~ sqrtSpAny * logCapdist + mountains_mean +
 			water_gc + distcoast + sqrtPopd + sqrtBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

summary(NB_sqrt_deaths_inter )

interaction_deaths_out <- summary(margins(NB_sqrt_deaths_inter, variables =
				  "sqrtSpAny", at = list( logCapdist =
							  seq(0,8, by=0.5))))

ggplot(interaction_deaths_out,
                    aes(x = exp(logCapdist), y = AME, ymin = lower, ymax = upper)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.1) +
  #goldenScatterCAtheme +
  theme(panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) +
  geom_hline(yintercept=0, linetype="dotted") +
  xlab('Distance from the Modern Capital') + ylab('AME')

NB_sqrt_state_based_min_inter <- glm.nb(sqrtState_based ~ sqrtSpAny * logCapdist + mountains_mean +
 			distcoast + barren_gc, 
		data = prio_grid_isd)

summary(NB_sqrt_state_based_min_inter )

NB_sqrt_state_based_inter <- glm.nb(sqrtState_based ~ sqrtSpAny * logCapdist + mountains_mean +
 			water_gc + distcoast + sqrtPopd + sqrtBDist + temp_sd + temp + prec_sd +
 			prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

summary(NB_sqrt_state_based_inter )

