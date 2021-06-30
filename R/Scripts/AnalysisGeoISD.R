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

#==============================================================================#
#	Loading Data							       #
#==============================================================================#	

load("../Data/GeoISDControls.Rdata")

#==============================================================================#
#	Creating New Variables						       #
#==============================================================================#

prio_grid_isd <- prio_grid_isd %>% mutate(
				logDeaths = log(deaths + 1),
				logState_based = log(state_based + 1),
				logPpp = log(gcp_ppp + 1),
				logSpAll = log(sp_os_i_sum + 1),
				logSpAny = log(sp_os_i_sum_any + 1),
				logCapdist = log(capdist + 1),
				logBDist = log(bdist3 + 1),
				logPopd = log(popd + 1),
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

      margins_out <- bind_rows(margins_out, m1_m, m2_m,m3_m, m3_m)
      models_out_j[[j]] <- list(m1, m2, m3)
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

extended_controls <- c('+ temp_sd + temp + prec_sd + prec_gpcc + barren_gc +
		       forest_gc')

all_controls <- c('+ logPopd + logBDist')

control_names <-c('Baseline', 'Extetended Controls', 'Full Model',
		 'Baseline', 'Extetended Controls', 'Full Model')

dvs <- c('logDeaths', 'logState_based', 'non_state', 'org3', 'logPpp')

captions <- c('Deaths (logged)', 'State based conflict events (logged)',
	      'Non state conflict events', 'Communal violence events', 'PPP
	      (logged)')

ivs <- c('logSpAll', 'logSpAny')

interactions <- c('SpAllXCapDist', 'SpAnyXCapDist')

#==============================================================================#
#	Analysis							       #
#==============================================================================#

linear_models <- GeoISDanalysis(dvs=dvs, ivs=ivs, controls=controls,
				 data=prio_grid_isd, test_label='Linear
				 Models')

interaction_models <- GeoISDanalysis(dvs=dvs, ivs=interactions, controls=controls,
				 data=prio_grid_isd, test_label='Interaction
				 Models')



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
	dep_var == 'logPpp' ~ 'PPP (logged)'))  %>% 
	rename(Test=test)


margins_main_plot <- ggplot(data = full_margins,
                            aes(x = factor_label, y = AME, 
				ymin = lower, ymax = upper, colour=Test)) +
  facet_wrap(~ facet_label) +
  xlab('') +
  labs(color = 'Models') +
  scale_color_manual(values=(brewer.pal(n = 4, name = "Dark2"))) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.1) +
  coord_flip() +
 # goldenScatterCAtheme +
  theme(panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) +
  geom_hline(yintercept=0, linetype="dotted") +
  guides(alpha="none") +
  labs(caption = 
         "Points are average marginal effects of a 1 SD increase in the independent variable with 95% confidence intervals.")

  margins_main_plot

# Marginal effects (interaction models)
fit <- lm(logDeaths ~ logSpAll * capdist + mountains_mean +
	  water_gc + distcoast + popd + bdist3 + temp_sd + temp + prec_sd +
	  prec_gpcc + barren_gc + forest_gc, data = prio_grid_isd)

predplot <- ggpredict(fit, terms = c("logSpAll", "capdist"))

ggplot(predplot, aes(x, exp(predicted))) +
  geom_line() +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.high)), alpha = .1)

