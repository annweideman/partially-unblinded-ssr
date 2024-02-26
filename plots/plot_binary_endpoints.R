#-------------------------------------------------------------------------------
# Purpose: Plots to compare the three formulas for estimating sample size for a 
# binary endpoint
# Author: Ann Marie Weideman
#-------------------------------------------------------------------------------

# Set working directory to current location
script_path <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Function to calculate sample sizes for the three formulas
calculate_sample_sizes <- function(k, p0C, p0T) {
  
  z_alpha <- qnorm(0.025, lower.tail = FALSE)
  z_beta <- qnorm(0.2, lower.tail = FALSE)
  delta0<-p0C-p0T
  p0<-(p0C+k*p0T)/(1+k)
  var0<-p0*(1-p0)
  var0C<-p0C*(1-p0C)
  var0T<-p0T*(1-p0T)
  
  # Assumes homogeneous variance 
  n0C_f1 <- (z_alpha + z_beta)^2 * (1 + 1/k) * var0 / delta0^2
  n0_f1 <- (1 + k) * n0C_f1
  
  # Assumes heterogeneous variance
  n0C_f2 <- (z_alpha + z_beta)^2 * (var0C + (1/k) * var0T) / delta0^2
  n0_f2 <- (1 + k) * n0C_f2
  
  # Assumes dual variance
  n0C_f3 <- (z_alpha * sqrt((1 + 1/k) * var0) + 
             z_beta * sqrt(var0C + (1/k) * var0T))^2 / delta0^2
  n0_f3 <- (1 + k) * n0C_f3
  
  return(c(n0_f1, n0_f2, n0_f3))
}

# Define the color palette to be colorblind friendly
palette <- brewer.pal(n = 3, name = "Set2")

k_range <- seq(1/3, 3, length.out = 100)
p0C_values <- 0.21
p0T_values <- c(0.05, 0.088, 0.113, 0.139, 0.2, 0.235, 0.276, 0.327, 0.4)

# Create a data frame with all combinations
sample_data <- expand.grid(k = k_range, p0C = p0C_values, p0T = p0T_values) %>%
  rowwise() %>%
  mutate(
    Homogeneous = calculate_sample_sizes(k, p0C, p0T)[1],
    Heterogeneous = calculate_sample_sizes(k, p0C, p0T)[2],
    Dual = calculate_sample_sizes(k, p0C, p0T)[3]
  ) %>%
  ungroup() %>%
  gather(key = "Variance", value = "sample_size", -k, -p0C, -p0T) %>%
  mutate(formula = factor(Variance, levels = c("Homogeneous", "Heterogeneous", "Dual"), 
                          labels = c("Homogeneous", 
                                     "Heterogeneous", 
                                     "Dual")))

# Function to create individual plots within the panel
make_plot <- function(data, p0C_val, p0T_val, show_x_axis_label = FALSE) {
  # Define custom line types
  line_types <- c("solid", "dashed", "dotted")
  
  p <- ggplot(data, aes(x = k, y = sample_size, group = Variance, color = Variance, linetype = Variance)) +
    geom_line(size = 1.2) +  # Increase line size for visibility
    scale_color_manual(values = palette) +
    scale_linetype_manual(values=line_types) +  # Use the custom line types
    scale_x_continuous(breaks = seq(0, 3, by = 0.5)) +
    theme_bw(base_size = 11) +  # Increase the base size for all text
    theme(
      plot.title = element_text(hjust = 0.5, size=9),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10, vjust = 0.5),
      plot.margin = unit(c(0.25, 0.25, 1, 1), "lines")
    ) +
  ggtitle(bquote(italic(p[0*"C"]) == .(p0C_val) * "," ~ 
                   "\n" ~ 
                   italic(p[0*"T"]) == .(p0T_val)))
  if (show_x_axis_label) {
    p <- p + labs(x = "k", y = expression("Total initial planned sample size ( " 
                                          ~ italic(n[0]) *")")) 
  }
  return(p)
}

# Generate a plot for each combination of p0C and p0T
plot_list <- list()
number_of_rows <- length(unique(sample_data$p0T))
counter <- 0

for (p0T_val in p0T_values) {
  for (p0C_val in  p0C_values) {
    counter <- counter + 1
    
    # Determine if this plot is on the bottom row
    show_x_axis_label <- (counter > (number_of_rows - 1) * 3)
    
    # Subset the data for the current combination of p0T and p0C
    data_subset <- subset(sample_data, p0T == p0T_val & p0C == p0C_val)
    
    # Add plot to the growing list of plots
    plot_list[[paste(p0T_val, p0C_val, sep = "_")]] <- 
      make_plot(data_subset, p0C_val=round(p0C_val,2),
                p0T_val=round(p0T_val,2), show_x_axis_label)
  }
}

# Combine the plots into a 3x3 patchwork
combined_plot <- wrap_plots(plot_list, ncol = 3) + 
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom")
print(combined_plot)

# Define the size and resolution of output image
width_in_inches <- 7
height_in_inches <- 6
dpi <- 600

# Open a PNG device
#jpeg(filename = paste0(script_path, "/output/trio_binary.jpg"), 
#     width = width_in_inches * dpi, height = height_in_inches * dpi, res = dpi)

# Print the combined plot
print(combined_plot)

# Add the common x-axis label using grid.text
require(grid)
grid.text(label = expression("Treatment allocation " * (italic(k))), 
          x = unit(0.55, "npc"),  # centered
          y = unit(0.1, "npc"),  # adjust this value as needed
          just = c("center", "bottom"),
          gp = gpar(fontsize = 10))

# Add the common y-axis label using grid.text
grid.text(label = expression("Total initial planned sample size  " * (italic(n[0]))), 
          x = unit(0.01, "npc"), y = unit(0.55, "npc"), rot = 90,
          gp = gpar(fontsize = 10))

# Close the PNG device, which saves the file
#dev.off()