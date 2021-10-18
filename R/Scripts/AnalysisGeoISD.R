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
library(pscl)
library(readr)
library(MASS)
library(spdep)
library(sf)
library(psych)
library(sidedata)
library(raster)

#==============================================================================#
#	Loading Data and functions					       #
#==============================================================================#	

load("../Data/GeoISDControls.Rdata")
source("goldenScatterCAtheme.r")

#==============================================================================#
#	Solving conflicts 						       #
#==============================================================================#		

conflict_prefer("filter", "dplyr")  
conflict_prefer("select", "dplyr")
conflict_prefer("describe", "psych")

#==============================================================================#
#	Creating New Variables						       #
#==============================================================================#

prio_grid_isd <- prio_grid_isd %>% mutate(
				sqrtDeaths = sqrt(deaths),
				sqrtState_based = sqrt(state_based),
				sqrtPpp = sqrt(gcp_ppp),
				sqrtSpAll = sqrt(sp_os_i_sum),
				sqrtSpAny = sqrt(sp_i_sum_any),
				sqrtCapdist = sqrt(capdist),
				sqrtBDist = sqrt(bdist3),
				sqrtPopd = sqrt(popd),
				logDeaths = log(deaths + 1),
				logState_based = log(state_based + 1),
				logPpp = log(gcp_ppp + 1),
				logSpAll = log(sp_os_i_sum + 1),
				logSpAny = log(sp_i_sum_any + 1),
				logCapdist = log(capdist),
				logBDist = log(bdist3 + 1),
				logPopd = log(popd + 1),
				sqrtorg3 = sqrt(org3),
				sqrtNonstate = sqrt(non_state),
				SpAllXCapDist = logSpAll*capdist,
				SpAnyXCapDist = logSpAny*capdist)

#==============================================================================#
#	Functions				    	             	       #
#==============================================================================#

#==============================================================================#
# Negative binomial models function

GeoISDanalysis <- function(dvs, ivs, controls, data, test_label){
  models_out_j <- NULL
  models_out <- NULL
  margins_out <- NULL
  for(i in 1:length(dvs)) {
    dv <- dvs[i] 
    for(j in 1:length(ivs)) {
      m1 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls)),
               data=data)
      m2 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls, extended_controls)),
               data=data)
      #m3 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls,
				    #extended_controls, climate_controls)),
               #data=data)
      iv_margins <- ivs[j] %>% strsplit(., "\\+") %>% unlist()
      m1_m <- summary(margins(m1, variables=iv_margins, change = "sd")) %>% 
        mutate(test='Baseline model', dep_var=dv)
      m2_m <- summary(margins(m2, variables=iv_margins, change = "sd")) %>% 
        mutate(test='Extended controls', dep_var=dv)
      #m3_m <- summary(margins(m3, variables=iv_margins, change = "sd")) %>% 
        #mutate(test='Climate controls', dep_var=dv)
      margins_out <- bind_rows(margins_out, m1_m, m2_m)
      models_out_j[[j]] <- list(m1, m2)
    } 
    models_out[[i]] <- models_out_j %>% purrr::flatten()
  }
  models_out <- models_out
  final <- list('margins'=margins_out, 'models'=models_out)
  return(final)
}

#==============================================================================#
# Same function, without margins to do the interaction models

GeoISD_interactions <- function(dvs, ivs, controls, data, test_label){
  models_out_j <- NULL
  models_out <- NULL
  margins_out <- NULL
  for(i in 1:length(dvs)) {
    dv <- dvs[i] 
    for(j in 1:length(ivs)) {
      m1 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls)),
               data=data)
      m2 <- glm.nb(as.formula(paste(dv, '~', ivs[j], controls, extended_controls)),
               data=data)
      models_out_j[[j]] <- list(m1, m2)
    } 
    models_out[[i]] <- models_out_j %>% purrr::flatten()
  }
  models_out <- models_out
  final <- list('models'=models_out)
  return(final)
}

#==============================================================================#
#	Descriptive Statistics                                         	       #
#==============================================================================#

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

# TODO: Print to file
summary(corrdata)

# Summary Statistics

skimmed <- describe(select(data.frame(prio_grid_isd), deaths, state_based, sp_os_i_sum, capdist))

latex(skimmed, title = "test", file = "../Output/skimmed.tex")

# stargazer and texreg?

d <- ggplot(prio_grid_isd, aes(deaths))

sb <- ggplot(prio_grid_isd, aes(state_based))

sp <- ggplot(prio_grid_isd, aes(sp_os_i_sum))

cd <- ggplot(prio_grid_isd, aes(capdist))

d + geom_histogram()

d + geom_dotplot(
	 dotsize = 0.5) +
	goldenScatterCAtheme
	
d + geom_area(stat = "bin")

d + geom_density(kernel = "gaussian")

#==============================================================================#
#	Variables							       #
#==============================================================================#

controls <- c('+ mountains_mean + water_gc + barren_gc + distcoast')

extended_controls <- c('+ logPopd + bdist3')

climate_controls <- c('+ temp_sd + temp + prec_sd + prec_gpcc')

coefs_cv <- list('logSpAll' = 'Precolonial state presence (log)', 
		'mountains_mean' = 'Mountainous terrain',
	     	'water_gc' = 'Water (%)', 
	 	'barren_gc' = 'Barren (%)', 
	      	'distcoast' = 'Distance to coast',
	      	'logPopd' = 'Population density (log)', 
		'bdist3' = 'Distance to border', 
		'temp_sd' = 'Temperature (SD)',
	      	'temp' = 'Temperature (mean)', 
		'prec_sd' = 'Precipitation (SD)', 
		'prec_gpcc' = 'Precipitation (mean)', 
		'sqrtSpAll' = 'Precolonial state presence (sqrt)')

control_names_int <-c('Baseline', 'Extended Controls')

control_names_cv <-c('Baseline', 'Extended Controls',
		    	'Baseline', 'Extended Controls')

control_names_cv_full <-c('Baseline', 'Extended Controls', 'Climate',
		    'Baseline', 'Extended Controls', 'Climate')

control_names <-c('Baseline', 'Extended Controls', 
		 'Baseline', 'Extended Controls',
		 'Baseline', 'Extended Controls')

dvs <- c('deaths', 'state_based')

cv_dvs <- c('non_state', 'org3')

captions <- c('Fatalities', 'State based conflict events')

cv_captions <- c('Non-state conflict events', 'Communal violence events')

captions_int <- c('Fatalities * Distance to capital', 'State based conflict events *
		  distance to capital')

ivs <- c('logSpAll', 'sqrtSpAll')

interactions <- c('sqrtSpAll * logCapdist')

#==============================================================================#
#	Set up for spatial analysis				       	       #
#==============================================================================#

prio_grid_shp <- st_read('../Data/PRIO-Grid/priogrid_cell.shp')

prio_grid <- read_csv('../Data/PRIO-Grid/priogridyv50-10.csv') %>% 
  as_tibble() %>% filter( (gwno >= 404 & gwno <= 626) | gwno == 651) %>% 
  group_by(gid) %>% 
  summarise(bdist3 = mean(bdist3),
  capdist = mean(capdist), excluded = mean(excluded), 
  	  temp_sd = sd(temp), gwno = last(gwno), temp = mean(temp), 
	  prec_sd = sd(prec_gpcc, na.rm = TRUE), prec_gpcc = mean(prec_gpcc))

prio_grid_static  <- read_csv('../Data/PRIO-Grid/PRIO-GRID Static Variables - 2021-06-04.csv') 

# Merging static and aggregated prio data
prio_grid  <- left_join(prio_grid, prio_grid_static, by = c("gid")) 

# Merging with the grid shape 
prio_grid <- left_join(prio_grid_shp, prio_grid, by = c("gid")) %>% 
  filter( (gwno >= 404 & gwno <= 626) | gwno == 651)

prio_sp <- as(prio_grid, Class = "Spatial")

nbQueen <- poly2nb(prio_sp, queen = T, row.names = prio_sp$gid)

lw = nb2listw(nbQueen, style = "W", zero.policy = TRUE)

moran.test(prio_grid_isd$deaths, listw = lw, zero.policy = TRUE)

prio_grid_isd$deaths_l <- lag.listw(lw, prio_grid_isd$deaths, zero.policy = T)

prio_grid_isd$state_based_l <- lag.listw(lw, prio_grid_isd$state_based, zero.policy = T)

prio_grid_isd$non_state_l <- lag.listw(lw, prio_grid_isd$non_state, zero.policy = T)

prio_grid_isd$org3_l <- lag.listw(lw, prio_grid_isd$org3, zero.policy = T)

#==============================================================================#
#	Analysis							       #
#==============================================================================#

nb_models <- GeoISDanalysis(dvs=dvs, ivs=ivs, controls=controls,
				 data=prio_grid_isd, test_label='Linear
				 Models')

interaction_models <- GeoISD_interactions(dvs=dvs, ivs=interactions, controls=controls,
				 data=prio_grid_isd, test_label='Interaction
				 Models')

#==============================================================================#
# Communal violence analysis & randomness


cv_nb_models <- GeoISDanalysis(dvs=cv_dvs, ivs=ivs, controls=controls,
				 data=prio_grid_isd, test_label='Linear
				 Models')

non_state_NB_mini <- glm.nb(non_state ~ sqrtSpAll , data = prio_grid_isd)
summary(non_state_NB_mini)

non_state_NB <- glm.nb(non_state ~ logSpAll + mountains_mean + water_gc + barren_gc +
		  distcoast + logPopd + bdist3, data = prio_grid_isd)
summary(non_state_NB)

org3_NB_mini_sqrt <- glm.nb(org3 ~ sqrtSpAll, data = prio_grid_isd)
summary(org3_NB_mini_sqrt)

org3_NB_pop_sqrt <- glm.nb(org3 ~ sqrtSpAll + logPopd, data = prio_grid_isd)
summary(org3_NB_pop_sqrt)

org3_NB_mini_log <- glm.nb(org3 ~ logSpAll, data = prio_grid_isd)
summary(org3_NB_mini_log)

# Test for NB or Poisson
org3_NB_mini<- glm.nb(org3 ~ sp_os_i_sum, data = prio_grid_isd)
summary(org3_NB_mini)

org3_P_mini <- glm(org3 ~ logSpAll, data = prio_grid_isd, family = poisson)
summary(org3_P_mini)

pchisq(2 * (logLik(org3_NB_mini) - logLik(org3_P_mini)), df = 1, lower.tail = FALSE)

org3_NB <- glm.nb(org3 ~ sqrtSpAll + mountains_mean + water_gc + barren_gc +
		  distcoast + logPopd + bdist3, data = prio_grid_isd)
summary(org3_NB)

#==============================================================================#
#	Marginal interacation plots 					       #
#==============================================================================#

lighten <- function (col, pct = 0.75, alpha = .8) 
{
    if (abs(pct) > 1) {
        print("Warning:  Error in Lighten; invalid pct")
        pcol <- col2rgb(col)/255
    }
    else {
        col <- col2rgb(col)/255
        if (pct > 0) {
            pcol <- col + pct * (1 - col)
        }
        else {
            pcol <- col * pct
        }
    }
    pcol <- rgb(pcol[1], pcol[2], pcol[3], alpha)
    pcol
}

ggorg3 <- ggeffect(org3_NB, terms = "sqrtSpAll [0:15 by = .5] ")

org3plot <- ggplot(ggorg3, aes(x^2, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
	      fill = lighten("blue")) +
  geom_line(color = "blue") +
  labs(x = 'Precolonial state presence', y = 'Communial violence events') +
  goldenScatterCAtheme

pdf("../Output/CommunalViolenceMargins.pdf",
    width = 10, height = 10/1.68)
org3plot
dev.off()

deathsMainInt <- glm.nb(deaths ~ sqrtSpAll * logCapdist + mountains_mean + water_gc + barren_gc +
		  distcoast + logPopd + bdist3, data = prio_grid_isd)
summary(deathsMainInt)

ggDeathsInt <- ggeffect(deathsMainInt, terms = c("sqrtSpAll [0:15]", "logCapdist
						 [1.309, 6.27, 7.817]"))

cols <- c("red", "green", "blue")

pastels <- NULL
for (i in 1:length(cols)) {
	pastels[i] <- lighten(cols[i])
}

ggplot(ggDeathsInt, aes(x^2, predicted, color = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group, linetype = NA)) +
  scale_fill_manual(values = pastels) +
  geom_line() +
  goldenScatterCAtheme

#==============================================================================#
# Regression tables
# TODO: Add coeficient maps for propper names in regression tables

xor (i in 1:length(dvs)) { 
  name <- dvs[i]
  filename <- paste("../Output/",name,".tex",sep="") 
  texreg(nb_models$models[[i]],
         file = filename,
         custom.model.names = control_names,
         #custom.coef.map = coefs,
         stars = c(0.001, 0.01, 0.05, 0.1), 
         sideways = T, use.packages = F, scalebox = 1,
         #custom.note = "",
         caption = captions[i],
	 label = name,
         table = T)
}

#==============================================================================#
# Interaction tabels

for (i in 1:length(dvs)) { 
  name <- dvs[i]
  filename <- paste("../Output/interaction_",name,".tex",sep="") 
  intlabel <- paste("interaction_",name,sep="")
  texreg(interaction_models$models[[i]],
         file = filename,
         custom.model.names = control_names_int,
         #custom.coef.map = coefs,
         stars = c(0.001, 0.01, 0.05, 0.1), 
         sideways = T, use.packages = F, scalebox = 1,
         #custom.note = "",
         caption = captions_int[i],
	 label = intlabel,
         table = T)
}

#==============================================================================#
# Communal violence regression tables

for (i in 1:length(cv_dvs)) { 
  name <- cv_dvs[i]
  filename <- paste("../Output/",name,".tex",sep="") 
  texreg(cv_nb_models$models[[i]],
         file = filename,
         custom.model.names = control_names_cv,
         custom.coef.map = coefs_cv,
         stars = c(0.001, 0.01, 0.05, 0.1), 
         sideways = T, use.packages = F, scalebox = 1,
         #custom.note = "",
         caption = cv_captions[i],
	 label = name,
         table = T)
}

#==============================================================================#
# Margin plots

full_margins <- nb_models[[1]] %>% 
  mutate(factor_label = case_when(
    factor == 'sqrtSpAll' ~ 'State Presence (square root)',
    factor == 'logSpAll' ~ 'State Presence (logged)',
    factor == 'sp_os_i_sum' ~ 'State Presence'),
	 facet_label = case_when(
	dep_var == 'deaths' ~ 'Total fatalities',
	dep_var == 'state_based' ~ 'State based conflict events'))  %>% 
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
  goldenScatterCAtheme + 
  geom_hline(yintercept=0, linetype="dotted") +
  guides(alpha="none") +
  labs(caption = "Points are average marginal effects of a 1 SD increase in the 
 independent variable with 95% confidence intervals.")

pdf("../Output/conflictMargins.pdf",
    width = 10, height = 10/1.68)
margins_main_plot
dev.off()

#==============================================================================#
# Communal violence margin plots

cv_full_margins <- cv_nb_models[[1]] %>% 
  mutate(factor_label = case_when(
    factor == 'sqrtSpAll' ~ 'State Presence (square root)',
    factor == 'logSpAll' ~ 'State Presence (logged)'),
	 facet_label = case_when(
	dep_var == 'non_state' ~ 'Non-state conflict events',
	dep_var == 'org3' ~ 'Communal violence events'))  %>% 
	rename(Test=test)

cv_margins_main_plot <- ggplot(data = cv_full_margins,
  aes(x = factor_label, y = AME, ymin = lower, ymax = upper, colour=Test)) +
  facet_wrap(~ facet_label) +
  xlab('') +
  labs(color = 'Models') +
  scale_color_manual(values=(brewer.pal(n = 4, name = "Dark2"))) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.1) +
  coord_flip() +
  theme(panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) +
  goldenScatterCAtheme + 
  geom_hline(yintercept=0, linetype="dotted") +
  guides(alpha="none") +
  labs(caption = "Points are average marginal effects of a 1 SD increase in the 
 independent variable with 95% confidence intervals.")

pdf("../Output/CommunalViolenceMargins.pdf",
    width = 10, height = 10/1.68)
cv_margins_main_plot
dev.off()

#=============================================================================#
# Communal violence negative binomial interaction models 

non_state_NB_inter <- glm.nb(non_state ~ sqrtSpAll * logCapdist + mountains_mean +
 			water_gc + distcoast + logPopd + bdist3, data = prio_grid_isd)

summary(non_state_NB_inter)

non_state_int_out <- summary(margins(non_state_NB_inter, variables =
				  "sqrtSpAll", at = list( logCapdist =
							  seq(0,7, by=0.2))))

non_state_int_plot <- ggplot(non_state_int_out,
                    aes(x = exp(logCapdist), y = AME, ymin = lower, ymax = upper)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.1) +
  goldenScatterCAtheme +
 # theme(panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) +
  geom_hline(yintercept=0, linetype="dotted") +
  xlab('Distance from the Modern Capital') + ylab('AME')

pdf("../Output/non_state_int_plot.pdf",
    width = 10, height = 10/1.68)
non_state_int_plot
dev.off()

org3_int <- glm.nb(org3 ~ sqrtSpAll * logCapdist + mountains_mean +
 			water_gc + distcoast + logPopd + bdist3, data = prio_grid_isd)

summary(org3_int)

org3_int_out <- summary(margins(org3_int, variables =
				  "sqrtSpAll", at = list( logCapdist =
							  seq(0,7.7, by=0.5))))

org3_int_plot <- ggplot(org3_int_out,
                    aes(x = exp(logCapdist), y = AME, ymin = lower, ymax = upper)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.1) +
  goldenScatterCAtheme +
 # theme(panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) +
  geom_hline(yintercept=0, linetype="dotted") +
  xlab('Distance from the Modern Capital') + ylab('AME')

pdf("../Output/org3_int_plot.pdf",
    width = 10, height = 10/1.68)
org3_int_plot
dev.off()

#==============================================================================#
#	ZINB  								       #
#==============================================================================#		

zinb_deaths <- zeroinfl(deaths ~ sqrtSpAll * logCapdist + mountains_mean +
 			water_gc + distcoast + logPopd + bdist3, data =
			prio_grid_isd, dist = "negbin")

summary(zinb_deaths)

summary(zinb_deaths$residuals)

# Because the regular model does not work I create standardized independent
# variables, and rerun the model.

prio_grid_isd <- prio_grid_isd %>% mutate(deaths_s = scale(deaths),
					  sqrtSpAll_s = scale(sqrtSpAll),
					  logCapdist_s = scale(logCapdist),
					  mountains_mean_s =
						  scale(mountains_mean),
					  water_gc_s = scale(water_gc),
					  distcoast_s = scale(distcoast),
					  logPopd_s = scale(logPopd),
					  bdist3_s = scale(bdist3))


zinb_deaths_s <- zeroinfl(deaths ~ sqrtSpAll_s * logCapdist_s +
			  mountains_mean_s +
 			water_gc_s + distcoast_s + logPopd_s + bdist3_s, data =
			prio_grid_isd, dist = "negbin")

summary(zinb_deaths_s)

summary(zinb_deaths$residuals)

#==============================================================================#
#	Transforming the main independent variable			       #
#==============================================================================#

# From a theoretical standpoint, what is the most appropriate way to model the
# influence of state presence on conflict?  Should we expect a linear
# relationship, one that increases rapidly and the tapers off (log), one that
# increases less rapidly but tapers off less quickly (sqrt), or one that
# increases slowly, then rapid before tapering off (s-shaped, like logistic)?

hist(prio_grid_isd$sp_os_i_sum)

hist(prio_grid_isd$logSpAll)

hist(prio_grid_isd$sqrtSpAll)

hist(prio_grid_isd$sp_os_i_sum^2)

x <- (0:100)

plot(sqrt(x))

plot(log(x))

plot(x)

#==============================================================================#
#	Excluding unpopulated cells					       #
#==============================================================================#

# TODO: Decide on whether to use 1500 pop.density or more modern. Am I only
# interested in places where there *were* any pople or, places where there
# *are* pople now? Probably current day.

# Draw plot of population denisty (log) anno 1500 where population is > 0

ggplot() +
	geom_sf(data = filter(prio_grid_isd, popd > 0),
            linetype = 0,
            aes_string(fill = "logPopd"),
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()

# Main Regressions without unpopulated cells

sb0 <- glm.nb(state_based ~ sqrtSpAll * capdist + mountains_mean + water_gc + barren_gc +
		    distcoast + logPopd + bdist3, data =
		    filter(prio_grid_isd, popd > 0))
summary(sb0)


ds0 <- glm.nb(deaths ~ sqrtSpAll * capdist + mountains_mean + water_gc + barren_gc +
		    distcoast + logPopd + bdist3, data =
		    filter(prio_grid_isd, popd > 0))
summary(ds0)

ns0 <- glm.nb(non_state ~ sqrtSpAll + mountains_mean + water_gc + barren_gc +
		    distcoast + logPopd + bdist3, data =
		    filter(prio_grid_isd, popd > 0))
summary(ns0)

org0 <- glm.nb(org3 ~ sqrtSpAll + mountains_mean + water_gc + barren_gc +
		    distcoast + logPopd + bdist3, data =
		    filter(prio_grid_isd, popd > 0))
summary(org0)

# Main Regressions without unpopulated cells and regional dummies

sb0 <- glm.nb(state_based ~ sqrtSpAll * capdist + mountains_mean + water_gc + barren_gc +
		    distcoast + logPopd + bdist3 + factor(region), data =
		    filter(prio_grid_isd, popd > 0))
summary(sb0)


# Throws an error
ds0 <- glm.nb(deaths ~ sqrtSpAll * capdist + mountains_mean + water_gc + barren_gc +
		    distcoast + logPopd + bdist3 + factor(region), data =
		    filter(prio_grid_isd, popd > 0))
summary(ds0)

ns0 <- glm.nb(non_state ~ sqrtSpAll + mountains_mean + water_gc + barren_gc +
		    distcoast + logPopd + bdist3 + factor(region), data =
		    filter(prio_grid_isd, popd > 0))
summary(ns0)

org0 <- glm.nb(org3 ~ sqrtSpAll + mountains_mean + water_gc + barren_gc +
		    distcoast + logPopd + bdist3 + temp_sd + temp + prec_sd +
		    prec_gpcc, data = filter(prio_grid_isd, popd > 0))
summary(org0)

#==============================================================================#
#	Danger Zone!						      	       #
#==============================================================================#

testmodel <- glm.nb(deaths ~ sqrtSpAll * capdist + mountains_mean + water_gc + barren_gc +
		    distcoast + logPopd + bdist3 + factor(region), data =
		    prio_grid_isd)
summary(testmodel)

# # The below command only needs to be run once
# side_download(country = "Uganda", year = 2010, marker = "ethnic", dest.dir =
#	      "../Data/", conv.hull = T)

uga.ethnic <- side_load(country = "Uganda", year = 2010, marker = "ethnic",
			source.dir = "../Data")

uga.ethnic.meta.df <- sidemap2data(uga.ethnic)
head(uga.ethnic.meta.df)

names(uga.ethnic) <- uga.ethnic.meta.df$groupname

plot(uga.ethnic)

head(uga.ethnic)

geoisd <- st_read('../../QGIS/Geo-ISD.shp')

geoisd_data <- read_csv('../../QGIS/Geo-ISD.csv',
                        cols(
                          gid = col_double(),
                          xcoord = col_double(),
                          ycoord = col_double(),
                          col = col_double(),
                          row = col_double(),
                          ISDID = col_character(),
                          COWID = col_double(),
                          name = col_character(),
                          year = col_double(),
                          lyear = col_integer(),
                          hyear = col_integer(),
                          source = col_character(),
                          coder = col_character(),
                          note = col_character(),
                          error = col_double(),
                          layer = col_character(),
                          path = col_character()
                        ),
                        col_names = T,
                        )

crs(uga.ethnic) <- crs(geoisd)

ugaSP <- as(uga.ethnic, 'SpatialGridDataFrame')

