############################################################################
########## Setting up ggplot theme and colors ##############################
############################################################################

# Credit to Christopher Adolph (with help from Kai Ping Leung) for creating this
# theme (faculty.washington.edu/cadolph).

library(ggplot2)      ## ggplot graphics package
library(RColorBrewer) ## for colors
library(lemon)        ## to "freshen up" ggplot2
library(ggrepel)      ## to avoid point overlap

goldenScatterCAtheme <- theme(
  ## Removes main plot gray background
  panel.background = element_rect(fill = "white"), 
  ## Golden rectangle plotting area (leave out for square)
  aspect.ratio = ((1 + sqrt(5))/2)^(-1), 
  ## All axes changes
  axis.ticks.length = unit(0.5, "char"),  # longer ticks
  ## Horizontal axis changes
  axis.line.x.top = element_line(linewidth = 0.2),    # thinner axis lines
  axis.line.x.bottom = element_line(linewidth = 0.2), # thinner axis lines
  axis.ticks.x = element_line(linewidth = 0.2),       # thinner ticks
  axis.text.x = element_text(color = "black", size = 12),
  ## match type of axis labels and titles
  axis.title.x = element_text(size = 12,
                              margin = ggplot2::margin(t = 7.5, r = 0, b = 0, l = 0)),
  ## match type; pad space between title and labels
  ## Vertical axis changes
  axis.ticks.y = element_blank(), # no y axis ticks (gridlines suffice)
  axis.text.y = element_text(color = "black", size = 12,
                             margin = ggplot2::margin(t = 0, r = -4, b = 0, l = 0)),
  ## match type of axis labels and titles, pad
  axis.title.y = element_text(size = 12,
                              margin = ggplot2::margin(t = 0, r = 7.5, b = 0, l = 0)),
  ## match type of axis labels and titles, pad
  ## Legend
  legend.key = element_rect(fill = NA, color = NA),
  ## Remove unhelpful gray background
  ## Gridlines (in this case, horizontal from left axis only
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(color = "gray45", linewidth = 0.2),
  ## Faceting (small multiples)
  strip.background = element_blank(),
  ## Remove unhelpful trellis-like shading of titles
  strip.text.x = element_text(size=12),  # Larger facet titles
  strip.text.y = element_blank(),        # No titles for rows of facets
  strip.placement = "outside",           # Place titles outside plot
  panel.spacing.x = unit(1.25, "lines"), # Horizontal space b/w plots
  panel.spacing.y = unit(1, "lines")     # Vertical space b/w plots
)
class(goldenScatterCAtheme)  ## Create this as a class
