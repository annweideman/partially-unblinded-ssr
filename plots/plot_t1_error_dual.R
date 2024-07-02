#-------------------------------------------------------------------------------
# Topic: Plots of type 1 error for output from script sims_t1_error_dual.R
# Author: Ann Marie Weideman
#-------------------------------------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)

# Modify this URL according to the output you want to plot (dual only for this script)
url<-"https://raw.githubusercontent.com/annweideman/partially-unblinded-ssr/main/output/results_p00.2_k0.5_seed1_dual.R"
load(url(url))

# Pivot the dataframe to a long format for easier plotting with ggplot2
results_long <- results_df %>%
  pivot_longer(cols = starts_with("t1_"),
               names_to = "analysis_type",
               values_to = "type_1_error") %>%
  mutate(analysis_type = factor(analysis_type, levels = c("t1_fixed_error", 
                                                          "t1_interim_error"),
                                labels = c("Fixed", "Interim")))

# Define a color-blind friendly palette
cb_palette <- RColorBrewer::brewer.pal(n = 3, name = "Set2")

# Plot smoothed lines for interim analysis and dashed line for fixed analysis
plot <- ggplot(results_long, aes(x = n0, y = type_1_error, 
                                 group = interaction(delta0, analysis_type))) +
  geom_smooth(data = results_long %>% filter(analysis_type == "Fixed" & delta0 == 0.05), 
              aes(linetype = analysis_type), 
              method = "loess", formula = y ~ x, se = FALSE, color = "black",
              span=0.4, fullrange=F) +
  geom_smooth(data = results_long %>% filter(analysis_type == "Interim"),
              aes(color = as.factor(delta0), linetype = analysis_type),
              method = "loess", formula = y ~ x, se = FALSE,
              span=0.4, fullrange=F) +
  scale_color_manual(values = cb_palette) +
  scale_linetype_manual(values = c("Fixed" = "twodash", "Interim" = "solid")) +
  labs(x = expression("Total initial planned sample size (" * n[0] * ")"),
       y = "Type 1 error rate",
       color = expression(Delta[0]),
       linetype = "Analysis Type") +
  theme_minimal() +
  guides(
    color = guide_legend(override.aes = list(linetype = "solid")),
    linetype = guide_legend(override.aes = list(color = "black"))
  ) +
  scale_x_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 100)) +
  scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, by = 0.01)) + 
  theme(
    panel.grid.major = element_line(color = "#F0F0F0"),
    panel.grid.minor = element_line(color = "#F8F8F8"),
    legend.position = "bottom",
    axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
    axis.title.y = element_text(margin = margin(r = 10, unit = "pt"))
  )

# Print the plot
print(plot)

# Save plot
#ggsave(filename = paste0(script_path, "/output/plot_t1_p00.2_k0.5_hetero.jpg"), 
#       plot, dpi = 600, width = 4, height = 4)
