
#==============================================================================#
#	Plotting state presence  					       #
#==============================================================================#

# Load data

# Creating variables vector

vars <- c("sp_sum_any", "sp_sum", "sp_i_sum_any", "sp_i_sum", "sp_b_sum_any",
          "sp_b_sum", "sp_b_i_sum_any", "sp_b_i_sum", "sp_o",
          "sp_o_i", "sp_os_sum_any", "sp_os_sum", "sp_os_i_sum_any", 
          "sp_os_i_sum")

# Creating plots

plot_list = list()
for(i in 1:length(vars)) {
  p = ggplot() +
     geom_sf(data = prio_grid_isd,
              linetype = 0,
              aes_string(fill = vars[i]), 
              show.legend = FALSE) + 
      scale_fill_viridis_c() +
      theme_minimal()
  plot_list[[i]] = p 
}

# Printing plots to file

for (i in 1:length(plot_list)) {
  title <- vars[i]
  filename <- paste(title,".pdf",sep="") 
  pdf(filename, width=15, height=15/1.618)
      print(plot_list[[i]])
  dev.off()
}

# Creating logged plots | currently not working

for(i in 1:length(vars)) {
  ivar = as.name(vars[i])
  var = paste0("prio_grid_isd$",ivar)
  p = ggplot() +
    geom_sf(data = prio_grid_isd,
            linetype = 0,
            aes_string(fill = log(var+1)), 
            show.legend = FALSE) + 
    scale_fill_viridis_c() +
    theme_minimal()
  plot_list[[i]] = p 
}

# Printing logged plots to file

for (i in 1:length(plot_list)) {
  title <- vars[i]
  filename <- paste("ln_",title,".pdf",sep="") 
  pdf(filename, width=15, height=15/1.618)
  print(plot_list[[i]])
  dev.off()
}


