#==============================================================================#
#	Loading Packages						       #
#==============================================================================#	

library(VennDiagram)
library(viridis)


#==============================================================================#
#	Loading Data and functions					       #
#==============================================================================#	

load("../Data/GeoISDControls.Rdata")
source("goldenScatterCAtheme.r")

#==============================================================================#
#	Plotting state presence  					       #
#==============================================================================#

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

# Main plot

sp_i_sum_any_plot <- ggplot() +
     geom_sf(data = prio_grid_isd,
              linetype = 0,
              aes_string(fill = sqrt(prio_grid_isd$sp_i_sum_any)),
              show.legend = FALSE) + 
      scale_fill_viridis_c() +
      theme_minimal()

pdf("../Output/sp_i_sum_any_plot.pdf",
    width = 10, height = 10/1.68)
sp_i_sum_any_plot
dev.off()

sp_sum_any_plot <- ggplot() +
     geom_sf(data = prio_grid_isd,
              linetype = 0,
              aes_string(fill = sqrt(prio_grid_isd$sp_sum_any)),
              show.legend = FALSE) + 
      scale_fill_viridis_c() +
      theme_minimal()

pdf("../Output/sp_sum_any_plot.pdf",
    width = 10, height = 10/1.68)
sp_sum_any_plot
dev.off()

# Create the Venn diagram
# some random data
foo <- c('a','b','c','d')
baa <- c('a','e','f','g')

# Generate plot
v <- venn.diagram(list(foo=foo, baa=baa),
                  fill = viridis(2),
                  alpha = c(0.5, 0.5), 
		  cat.cex = 2, 
		  cat.pos = 1,
		  cex=1.2,
                  filename=NULL)

# Over-write labels (5 to 7 chosen by manual check of labels)
# in left only
v[[5]]$label  <- paste("Non-violent dissent")  
# in right only
v[[6]]$label <- paste("Interstate conflict and
all non-state conflict 
including communal violence")  
# intesection
v[[7]]$label <- paste("Violent dissent: 
intrastate- extrastate-
and internationalized 
conflicts")  
# Left blob
v[[8]]$label <- paste("Maximalist dissent")
# Right blob
v[[9]]$label <- paste("Organized violence")

ggsave("../Output/venn.pdf", width = 10, height = 10/1.683, v)

