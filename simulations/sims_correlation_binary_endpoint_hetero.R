#------------------------------------------------------------------------------
# Topic: Simulation studies to investigate correlation between interim effect 
# size and numerator of sample size re-estimation formula for a binary 
# endpoint that assumes heterogeneity of variances.
# Author: Ann Marie Weideman
#------------------------------------------------------------------------------

# Parameter ranges
p1_range <- seq(0.1, 0.9, by = 0.01) 
k_values <- seq(0.5,2,0.1)
n_intC_values <- 10000
iterations <- 100

# Data frame to store results
results <- expand.grid(p1 = p1_range, k = k_values, n_intC = n_intC_values)
results$correlation <- NA # Add a column to store correlations

# Run simulation for each combination of parameters
for (idx in 1:nrow(results)) {
  p_1C <- results$p1[idx]
  p_1T <- results$p1[idx]
  k <- results$k[idx]
  n_intC <- results$n_intC[idx]
  n_intT <- k * n_intC
  n_int <- n_intC + n_intT
  
  delta_1_hats <- numeric(iterations)
  gamma_1_hats <- numeric(iterations)
  
  for (i in 1:iterations) {
    Y_1C <- rbinom(n = n_intC, size = 1, prob = p_1C)
    Y_1T <- rbinom(n = n_intT, size = 1, prob = p_1T)
    
    p_1C_hat <- mean(Y_1C)
    p_1T_hat <- mean(Y_1T)
    
    delta_1_hats[i] <- p_1C_hat - p_1T_hat
    gamma_1_hats[i] <- p_1C_hat * (1 - p_1C_hat) + (1 / k) * p_1T_hat * (1 - p_1T_hat)
  }
  
  results$correlation[idx] <- cor(delta_1_hats, gamma_1_hats)
}

library(plotly)

color_within_range <- 'rgb(255,165,0)'  # Orange for -0.05 < r < 0.05
color_outside_range <- 'rgb(0,0,255)'  # Blue for r ≤ -0.05 and r ≥ 0.05

# 3D Plot
fig <- plot_ly() %>%
  add_trace(
    data = subset(results, correlation >= -0.05 & correlation <= 0.05),
    x = ~p1, y = ~k, z = ~correlation,
    type = 'scatter3d', mode = 'markers',
    marker = list(size = 5, opacity = 0.5, color = color_within_range),
    name = '-0.05 < r < 0.05'  
  ) %>%
  add_trace(
    data = subset(results, correlation < -0.05 | correlation > 0.05),
    x = ~p1, y = ~k, z = ~correlation,
    type = 'scatter3d', mode = 'markers',
    marker = list(size = 5, opacity = 0.5, color = color_outside_range),
    name = 'r ≤ -0.05, r ≥ 0.05'  
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = 'p<sub>1C</sub> = p<sub>1T</sub>',
                   tickvals = seq(0, 1, by = 0.25), range = c(0, 1)), 
      yaxis = list(title = 'k', tickvals = seq(0.5, 2, by = 0.5), range = c(0.5, 2)),
      zaxis = list(title = 'Correlation (r)')
    ),
    legend = list(
      orientation = "h", 
      x = 0.5,  
      xanchor = 'center',  
      y = -0.1,  
      yanchor = 'top'  
    )
  )

fig
