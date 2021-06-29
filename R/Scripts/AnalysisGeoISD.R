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

#==============================================================================#
#	Loading Data							       #
#==============================================================================#	

load("../Data/GeoISDControls.Rdata")

#==============================================================================#
#	Functions				    	             	       #
#==============================================================================#

GeoISDanalysis <- function(dvs, ivs, controls, data, test_label){
  models_out_j <- NULL
  models_out <- NULL
  
  for(i in 1:length(dvs)) {
    
    dv <- dvs[i] 
    
    for(j in 1:length(ivs)) {
      m1 <- lm(as.formula(paste(dv, '~', ivs[j], controls)),
               data=data)
      m2 <- lm(as.formula(paste(dv, '~', ivs[j], controls, extended_controls)),
               data=data)
      m3 <- lm(as.formula(paste(dv, '~', ivs[j], controls, extended_controls, all_controls)),
               data=data)

      models_out_j[[j]] <- list(m1, m2, m3)
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

# Correlation Matrix Plot

corrdata <- prio_grid_isd %>% dplyr::select(bdist3, capdist, excluded, temp_sd,
					    temp, prec_sd, prec_gpcc, barren_gc,
					    forest_gc, mountains_mean, water_gc,
					    sp_os_sum_any, sp_os_sum,
					    sp_os_i_sum_any, sp_os_i_sum,
					    distcoast, popd, state_based,
					    non_state, org3)

corrdata_mat <- cor(corrdata, method = c("spearman"), 
                     use="pairwise.complete.obs")

corrdata_mat <- rcorr(as.matrix(corrdata))

corrplot(corrdata_mat$r, method='color', diag=F, addCoef.col = "black", 
         p.mat = corrdata_mat$P, insig = "blank", tl.col = "black",
         tl.srt = 45)

# Summary Statistics

summary(corrdata)

#==============================================================================#
#	Variables							       #
#==============================================================================#

controls <- c('+ capdist + mountains_mean + water_gc + distcoast')

extended_controls <- c('+ temp_sd + temp + prec_sd + prec_gpcc + barren_gc +
		       forest_gc')

all_controls <- c('+ popd + bdist3')

dvs <- c('state_based', 'non_state', 'org3', 'gcp_mer', 'gcp_ppp')

ivs <- c('sp_os_sum_any', 'sp_os_sum', 'sp_os_i_sum_any', 'sp_os_i_sum')

#==============================================================================#
#	Analysis							       #
#==============================================================================#

initial_models <- GeoISDanalysis(dvs=dvs, ivs=ivs, controls=controls,
				 data=prio_grid_isd, test_label='Initial
				 Models')
#model_names <- c('

for (i in 1:length(dvs)) { 
  name <- dvs[i]
  filename <- paste("../Tables and Figures/",name,".tex",sep="") 
  texreg(initial_models$models[[i]],
         file = filename,
         #custom.model.names = control_names,
         #custom.coef.map = coefs,
         stars = c(0.001, 0.01, 0.05, 0.1), 
         sideways=T, use.packages = F, scalebox = 0.8,
         #custom.note = "Reference region is Eastern Europe and Central Asia.",
         caption="", 
         table = T)
}

