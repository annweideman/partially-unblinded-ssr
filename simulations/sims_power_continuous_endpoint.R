
#------------------------------------------------------------------------------
# Purpose : Simulation studies to investigate sample size and power for a 
#           continuous endpoint under two methods: equal variances for a blinded 
#           interim and dual variances for a partially-unblinded interim.
# Author  : Ann Marie Weideman
# Date    : 27-Oct-2024
#------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(parallel)
library(knitr)

alpha <- 0.05
beta <- 0.8 #nominal power
z_alpha2 <- qnorm(1 - alpha / 2)
z_beta <- qnorm(beta)

# Function to calculate sample size for homogeneity of variances under a blinded interim
sample_size_equal_variances <- function(var1, delta0, k) {
  n1C <- ((1 + 1 / k) * var1 * (z_alpha2 + z_beta)^2) / (delta0^2)
  n1T <- k*n1C
  return(c(ceiling(n1C), ceiling(n1T)))  
}

# Function to calculate sample size for dual variances under a partially-unblinded 
# interim
sample_size_dual_variances <- function(varC, var1T, delta0, k) {
  term_1 <- z_alpha2 * sqrt(var1T + (1 / k) * varC)
  term_2 <- z_beta * sqrt(varC + (1 / k) * var1T)
  n1C <- (term_1 + term_2)^2 / delta0^2
  n1T <- k*n1C
  return(c(ceiling(n1C), ceiling(n1T)))  
}

# Function to simulate a total of num_sims estimates of power
run_power_sims <- function(norm1C, norm1T, muC, muT, 
                           varC, var1T, delta0, n1C, 
                           nint, num_sims, k, se_type) {
  n1T <- k * n1C
  successes <- numeric(num_sims)
  
  # Compute the remaining sample size to enroll post-interim.
  # If re-estimated sample size is over twice original sample size, then allow
  # the remaining sample size to be equal to the interim sample size.
  # Otherwise, compute the difference between re-estimated sample size and interim 
  # sample size to determine remaining sample size.
  n2C_size <- if(n1T+n1C>=4*nint){nint/(k+1)}
              else{max(0, ceiling(n1C - nint / (1 + k)))}
  n2T_size <- if(n1T+n1C>=4*nint){k*nint/(k+1)}
              else{max(0, ceiling(n1T - k * nint / (1 + k)))}

  norm2C_mat <- if (n2C_size > 0) matrix(rnorm(num_sims * n2C_size, mean = muC,
                                               sd = sqrt(varC)), ncol = num_sims) 
  else matrix(0, nrow = 0, ncol = num_sims)
  norm2T_mat <- if (n2T_size > 0) matrix(rnorm(num_sims * n2T_size, mean = muT,
                                               sd = sqrt(var1T)), ncol = num_sims) 
  else matrix(0, nrow = 0, ncol = num_sims)
  
  # Compute estimate of power
  for (i in seq_len(num_sims)) {
    normC <- c(norm1C, norm2C_mat[, i])
    normT <- c(norm1T, norm2T_mat[, i])
    
    muC <- mean(normC)
    muT <- mean(normT)
    
    var_hat <- (var(normC)+k*var(normT))/(1+k)
    se_hat<-sqrt(var_hat * (1 / n1C + 1 / n1T))
    
    df <- n1C + n1T - 2
    
    t_stat <- (muT - muC) / se_hat
    successes[i] <- abs(t_stat) >= qt(1 - alpha / 2, df = df)
  }
  
  return(mean(successes))
}

# Function to simulate a total of num_studies number of studies 
run_simulation_parallel <- function(muC, muT, var0C, var0T, 
                                    varC, var1T, delta0, num_sims, k) {
  
  nint<- sum(sample_size_dual_variances(var0C, var0T, delta0, k))/2
  
  norm1C <- rnorm(nint / (k + 1), muC, sqrt(varC))
  norm1T <- rnorm(k * nint / (k + 1), muT, sqrt(var1T))
  muC_hat<-mean(norm1C)
  muT_hat<-mean(norm1T)
  var1C_hat <- var(norm1C)
  var1T_hat <- var(norm1T)
  var1<-var(c(norm1C,norm1T))
  
  n_method1 <- sample_size_equal_variances(var1, delta0, k)
  n_method2 <- sample_size_dual_variances(var1C_hat, var1T_hat, delta0, k)
  
  power_method1 <- run_power_sims(norm1C, norm1T, muC, muT, varC, 
                                   var1T, delta0, n_method1[1], nint, num_sims, k, 1)
  power_method2<- run_power_sims(norm1C, norm1T, muC, muT, varC, 
                                var1T, delta0, n_method2[1], nint, num_sims, k, 2)
  
  results <- data.table(
    Method = c("Homogeneity of Variances (blinded interim)", 
               "Duality of Variances (partially unblinded interim)"),
    Variance.Ratio = rep(var1T / varC, 2),
    k = rep(k, 2),
    Total.Sample.Size = c(sum(n_method1), sum(n_method2)),
    Power1 = c(power_method1, power_method2)
  )
  
  return(results)
}

# Parallelization to speed up sims
test_variances_and_ks_parallel <- function(variance_ratios, k_values, var0C, varC, 
                                           muC, muT, 
                                           delta0, num_sims, num_studies) {
  cl <- makeCluster(detectCores() - 1) 
  clusterSetRNGStream(cl, 12345)  # Ensures reproducibility across parallel workers
  clusterExport(cl, list("run_simulation_parallel", "sample_size_equal_variances", 
                         "sample_size_dual_variances", "run_power_sims", "rnorm", 
                         "alpha", "qt", "z_alpha2", "z_beta", "num_sims", 
                         "var0C","varC", "delta0", "k_values", "variance_ratios",
                         "muC", "muT"))  
  clusterEvalQ(cl, library(data.table))
  
  all_results <- parLapply(cl, seq_len(num_studies), function(rep) {
    results_list <- list()
    for (k in k_values) {
      for (var_ratio in variance_ratios) {
        var0T <- var0C * var_ratio
        var1T <- varC * var_ratio
        results_list[[length(results_list) + 1]] <- 
          run_simulation_parallel(muC, muT, var0C, var0T,
                                  varC, var1T, delta0, num_sims, k)
      }
    }
    return(rbindlist(results_list))
  })
  
  stopCluster(cl)
  return(rbindlist(all_results))
}

variance_ratios <- c(1, 2, 4, 9)
k_values <- c(1, 2, 3)
var0C <- 1 #anticipated variance in control arm
varC <- 1 #true variance in control arm
muC <- 0.25 #true mean in control arm
muT <- 0.75 #true mean in tx arm
delta0 <- 0.5 #anticipated effect size

num_studies <- 1000
num_sims <- 1000

results <- test_variances_and_ks_parallel(variance_ratios, k_values, var0C, varC, 
                                          muC, muT, delta0, num_sims, 
                                          num_studies)

library(dplyr)
average_results <- results %>%
  group_by(Variance.Ratio, k, Method) %>%
  summarize(
    Total.Sample.Size = round(mean(Total.Sample.Size)),
    Power = round(mean(Power1), 2),
    .groups = 'drop'
  ) %>%
  arrange(Variance.Ratio, k, Method)

# Arrange the data row-wise
rowwise_table <- average_results %>%
  arrange(Variance.Ratio, k, Method)

# Save the row-wise formatted table to a file
library(writexl)
write_xlsx(rowwise_table, 
           "C:\\Users\\anndo\\Dropbox\\Dissertation\\manuscripts\\paper3\\stats in med submission\\reviewer feedback\\output_sims_continuous.xlsx")

# Create a nested table
library(flextable)
nested_table <- flextable(rowwise_table) %>%
  merge_v(j = c("Variance.Ratio", "k")) %>% 
  valign(j = c("Variance.Ratio", "k"), valign = "top") %>% 
  align(j = c("Variance.Ratio", "k"), align = "left", part="all") %>% 
  align(j = c("Total.Sample.Size", "Power"), align = "left", part = "all") %>% 
  align(j = "Method", align = "left", part = "all") %>% 
  autofit() %>% 
  set_header_labels(
    Variance.Ratio = "Variance Ratio",
    Total.Sample.Size = "Total Sample Size"
  )

nested_table
