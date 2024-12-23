#------------------------------------------------------------------------------
# Purpose : Simulation studies to investigate sample size and power for a binary
#           endpoint under two methods: equal variances for a blinded 
#           interim and dual variances for a partially-unblinded interim.
# Author  : Ann Marie Weideman
# Date    : 27-Oct-2024
#------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(knitr)
library(kableExtra)
library(parallel)

# Function to calculate sample size for homogeneity of variances under a blinded interim
sample_size_equal_variances <- function(var1, delta0, k) {
  n1C <- ((1 + 1 / k) * var1 * (z_alpha2 + z_beta)^2) / (delta0^2) + (1+k)/(k*delta0)
  n1T <- k * n1C
  return(c(ceiling(n1C), ceiling(n1T)))
}

# Function to calculate sample size for duality of variances under a partially-unblinded 
# interim
sample_size_dual_variances <- function(var1, var1C, var1T, delta0, k) {
  term_1 <- z_alpha2 * sqrt((1 + 1 / k) * var1)
  term_2 <- z_beta * sqrt(var1C + (1 / k) * var1T)
  n1C <- (term_1 + term_2)^2 / delta0^2 + (1+k)/(k*delta0)
  n1T <- k * n1C
  return(c(ceiling(n1C), ceiling(n1T)))
}

# Function to simulate a total of num_sims estimates of power
run_power_sims <- function(pC, pT, var1C, var1T, delta0, binom1C, binom1T,
                           nint, n1, k, num_sims) {
  n1C <- n1/(1+k)
  n1T <- n1 - n1C
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
  
  binom2C_mat <- if (n2C_size > 0) matrix(rbinom(num_sims * n2C_size, size = 1,
                                                 prob=pC), ncol = num_sims) 
  else matrix(0, nrow = 0, ncol = num_sims)
  binom2T_mat <- if (n2T_size > 0) matrix(rbinom(num_sims * n2T_size, size = 1,
                                                 prob=pT), ncol = num_sims) 
  else matrix(0, nrow = 0, ncol = num_sims)
  
  # Compute estimate of power
  for (i in seq_len(num_sims)) {
  
    binomC <- c(binom1C, binom2C_mat[, i])
    binomT <- c(binom1T, binom2T_mat[, i])
    
    pC <- mean(binomC)
    pT <- mean(binomT)
    
    p_hat <- (pC+k*pT)/(1+k)
    se_hat<-sqrt(p_hat*(1-p_hat) * (1 / n1C + 1 / n1T))
    
    chisqr_stat <- ((pC - pT) / se_hat)^2
    successes[i] <- chisqr_stat >= qchisq(p=1 - alpha, df = 1)
  }
  
  return(mean(successes))
}

# Main simulation to simulate a total of num_studies number of studies 
run_simulation <- function(p0C, p0T, pC, pT, delta0, num_studies, num_sims, k) {
  
  results <- data.table(
    Method = character(),
    Var.Ratio = numeric(),
    k = integer(),
    Total.Sample.Size = integer(),
    Power = numeric()
  )
  for(i in 1:length(p0C)){
    var0C<-p0C[i]*(1-p0C[i])
    var0T<-p0T[i]*(1-p0T[i])
    var0<-(p0C[i]+k*p0T[i])/(1+k)
    nint<- sum(sample_size_dual_variances(var0, var0C, var0T, delta0, k))/2
    
    n_equal<-n_dual<-power_equal<-power_dual<-c()
    
    for(j in 1:num_studies){
      var1C <- pC[i]*(1-pC[i])
      var1T <- pT[i]*(1-pT[i])
      binom1C <- rbinom(nint / (k + 1), 1, pC)
      binom1T <- rbinom(k * nint / (k + 1), 1, pT)
      binom1 <- c(binom1C,binom1T)
      p1C_hat <- mean(binom1C)
      p1T_hat <- mean(binom1T)
      p1_hat <- mean(binom1)
      var1C_hat <- p1C_hat*(1-p1C_hat)
      var1T_hat <- p1T_hat*(1-p1T_hat)
      var1_hat <- p1_hat*(1-p1_hat)
      
      if(sum(binom1C)>2){
      n_equal <- c(n_equal, sum(sample_size_equal_variances(var1_hat, delta0, k)))
      n_dual <- c(n_dual, sum(sample_size_dual_variances(var1_hat, 
                                                     var1C_hat, var1T_hat, delta0, k)))
      }else{n_equal<-c(n_equal, nint*2)
            n_dual<-c(n_dual, nint*2)
      }
      
      power_equal <- c(power_equal, run_power_sims(pC[i], pT[i], var1C, var1T, delta0, 
                                                   binom1C, binom1T, nint,
                                                   n_equal[j], k, num_sims))
      
      power_dual <- c(power_dual, run_power_sims(pC[i], pT[i], var1C, var1T, delta0, 
                                                 binom1C, binom1T, nint,
                                                 n_dual[j], k, num_sims))
    }
    
    temp_results <- data.table(
      Method = c("Homogeneity of variances (blinded interim)", "Duality of variances (partially-unblinded interim)"),
      Var.Ratio = rep(var0T/var0C, 2), 
      k = rep(k, 2), 
      Total.Sample.Size = c(mean(n_equal), mean(n_dual)),
      Power = c(mean(power_equal), mean(power_dual))
    )
    
    results <- rbind(results, temp_results)
  }
  return(results)
}

k_values <- c(1, 2, 3)
alpha <- 0.05
beta <- 0.8 #nominal power
z_alpha2 <- qnorm(1 - alpha / 2)
z_beta <- qnorm(beta)

# Anticipated event rates
p0C <- c(0.05, 0.075, 0.1)
p0T <- c(0.15, 0.175, 0.2)
# True event rates
pC <- c(0.075, 0.1, 0.125)
pT <- c(0.175, 0.2, 0.225)
# Anticipated effect size
delta0 <- 0.1

num_studies <- 1000
num_sims <- 1000

#Parallelization to speed up sims
num_cores <- detectCores() - 1 
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, 12345)  # Ensures reproducibility across parallel workers
clusterExport(cl, c("run_simulation", "sample_size_equal_variances", 
                    "sample_size_dual_variances",
                    "run_power_sims", "p0C","p0T", "pC", "pT", "delta0", "num_studies", 
                    "num_sims", "z_alpha2", "z_beta", "alpha"))
clusterEvalQ(cl, library(data.table))

results <- rbindlist(parLapply(cl, k_values, function(k) {
  run_simulation(p0C, p0T, pC, pT, delta0, num_studies, num_sims, k)
}))

stopCluster(cl)

average_results <- results %>%
  mutate(Var.Ratio = round(Var.Ratio, 2)) %>%  # Adjusting Var.Ratio to 1 decimal place
  group_by(Var.Ratio, k, Method) %>%
  summarize(
    Total.Sample.Size = round(mean(Total.Sample.Size)),
    Power = round(mean(Power), 2),
    .groups = 'drop'
  ) %>%
  arrange(Var.Ratio, k, Method)

# Arrange the data row-wise
rowwise_table <- average_results %>% arrange(Var.Ratio, k, Method)

# Save the row-wise formatted table to a file
write_xlsx(rowwise_table, 
           "C:\\Users\\anndo\\Dropbox\\Dissertation\\manuscripts\\paper3\\stats in med submission\\reviewer feedback\\output_sims_power_binary_endpoint.xlsx")

# Create a nested table
nested_table <- flextable(rowwise_table) %>%
  merge_v(j = c("Var.Ratio", "k")) %>% 
  valign(j = c("Var.Ratio", "k"), valign = "top") %>% 
  align(j = c("Var.Ratio", "k"), align = "left", part="all") %>% 
  align(j = c("Total.Sample.Size", "Power"), align = "left", part = "all") %>% 
  align(j = "Method", align = "left", part = "all") %>% 
  autofit() %>% 
  set_header_labels(
    Var.Ratio = "Variance Ratio",
    Total.Sample.Size = "Total Sample Size"
  )

print(nested_table)
