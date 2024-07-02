#-------------------------------------------------------------------------------
# Topic: Investigation of type 1 error inflation under a partially-unblinded
#        paradigm
# Author: Ann Marie Weideman
#-------------------------------------------------------------------------------

#set working directory to current location
script_path<-setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

k <- 2 #allocation to tx
n0C.vec <- seq(10, 700, 1)  # vector of initial sample sizes in ctrl arm
n0 <- n0C.vec + k * n0C.vec # vector of total initial sample sizes
delta0.vec <- c(0.05, 0.1, 0.2) # vector of effect sizes
p0 <- 0.4 # pooled event rate
p0T<-p0-delta0.vec/(1+k) # event rate in tx arm
p0C<-p0+k*delta0.vec/(1+k) # event rate in ctrl arm

alpha <- 0.05 # nominal Type 1 error rate
z_alpha <- qnorm(alpha/2, lower.tail = FALSE) # z_alpha ~ 1.96
z_beta <- qnorm(0.2, lower.tail = FALSE) # z_beta ~ 0.84

nsims <- 10000
results <- list()
pval.fixed.df <- matrix(NA, nrow = nsims, ncol = length(n0C.vec))

for (delta0 in delta0.vec) {
  set.seed(123456)
  t1.fixed.error <- t1.interim.error <- numeric(length(n0C.vec))
  pval.interim.df <- matrix(NA, nrow = nsims, ncol = length(n0C.vec))
  colnames(pval.fixed.df) <- colnames(pval.interim.df) <- paste0("n0C=", n0C.vec)
  
  for (j in seq_along(n0C.vec)) {
    set.seed(123456)
    n0C <- n0C.vec[j]
    n0T <- k * n0C # original planned sample size in tx arm

    #------------------------------------------------------------------
    # Type 1 error for fixed sample size design (no interim analysis)
    # Will be identical for each delta0, so only run once
    #------------------------------------------------------------------
    
    if(delta0==min(delta0.vec)){
      for(i in 1:nsims){
        y0C<-rbinom(n0C,1,p0) # events at final analysis in ctrl arm under H0
        y0T<-rbinom(n0T,1,p0) # events at final analysis in tx arm under H0
        p0C<-mean(y0C) # event rate in ctrl arm under H0
        p0T<-mean(y0T) # event rate in tx arm under H0
        
        # Continuity-correct test statistic if sample size < 100
        if(n0[j]>=100){yates_chisqr0<-0}else{yates_chisqr0<-1/2*(1/n0T+1/n0C)}
        chisqr=(abs(p0C-p0T)-yates_chisqr0)^2/(p0C*(1-p0C)/n0C+p0T*(1-p0T)/n0T)
        
        p.val<-pchisq(chisqr,1,lower.tail=F)
        pval.fixed.df[i,j] <- p.val
        }
    }
    
    #----------------------------------------------------------
    # Type 1 error for partially-unblinded interim analysis
    #----------------------------------------------------------
  
    for(i in 1:nsims){
      #set.seed(i)
      nint_C<-round(n0C/2) # number at enrolled at interim analysis in ctrl arm
      nint_T<-k*nint_C
      y1C<-rbinom(nint_C,1,p0) # events at interim analysis in ctrl arm under H0
      y1T<-rbinom(nint_T,1,p0) # events at interim analysis in tx arm under H0
      p1C<-sum(y1C)/nint_C # proportion of events in ctrl arm under H0
      p1T<-sum(y1T)/nint_T # proportion of events in tx arm under H0
  
      # Re-estimate sample size in ctrl arm at interim analysis
      # Leave denominator as a function of stage 0 proportions
      # Utilize continuity correction if n0<100
      if(n0[j]>=100){yates_n1<-0}else{yates_n1<-(1+k)/(k*abs(delta0))}
      n1C<-(z_alpha+z_beta)^2*(p1C*(1-p1C)+(1/k)*p1T*(1-p1T))/delta0^2 + yates_n1
  
      # Re-estimate sample size in tx arm at interim analysis
      n1T<-k*n1C
      n1C<-round(n1C)
      n1T<-round(n1T)
  
      # If the re-estimated sample size exceeds the initial estimate
      if(n1C>n0C){
        # Generate remaining (stage 2) binomial events
        y2C<-rbinom((n1C-nint_C),1,p0)
        y2T<-rbinom((n1T-nint_T),1,p0)
      }else{
        # Generate remaining (stage 2) binomial events
        y2C<-rbinom((n0C-nint_C),1,p0)
        y2T<-rbinom((n0T-nint_T),1,p0)
        # Use original estimate since n1C < n0C
        n1C<-n0C
        n1T<-n0T
        }
  
      yT<-c(y1T,y2T) # vector of stage 1 and 2 events in tx arm
      yC<-c(y1C,y2C) # vector of stage 1 and 2 events in ctrl arm
      pT<-mean(yT) # end-of-study event rate in tx arm
      pC<-mean(yC) # end-of-study event rate in ctrl arm
      
      # Continuity-correct test statistic if n0 < 100
      if(n0[j]>=100){yates_chisqr1<-0}else{yates_chisqr1<-1/2*(1/n1T+1/n1C)}
      chisqr=(abs(pC-pT)-yates_chisqr1)^2/(pC*(1-pC)/n1C + pT*(1-pT)/n1T)
      p.val<-pchisq(chisqr,1,lower.tail=F)
      pval.interim.df[i,j] <- p.val
    }
  
    id.fixed.nan<-which(!is.nan(pval.fixed.df[,j]))
    id.interim.nan<-which(!is.nan(pval.interim.df[,j]))
    # Store the results for the fixed and interim T1 errors
    t1.fixed.error[j] <- length(which(pval.fixed.df[id.fixed.nan, j] <= alpha)) / 
                         length(id.fixed.nan)
    t1.interim.error[j] <- length(which(pval.interim.df[id.interim.nan, j] <= alpha)) /
                           length(id.interim.nan)
    print(paste("delta0 =", delta0, "j =", j))
    }
  
  results[[as.character(delta0)]] <- list(
    n0 = n0,
    t1.fixed.error = t1.fixed.error,
    t1.interim.error = t1.interim.error
  )
}

results_df <- data.frame(delta0 = numeric(),
                         n0C = numeric(),
                         t1_fixed_error = numeric(),
                         t1_interim_error = numeric())

# Loop across the effect size to fill the dataframe
for (delta0 in names(results)) {
  temp_df <- data.frame(delta0 = as.numeric(delta0),
                        n0 = results[[delta0]]$n0,
                        t1_fixed_error = results[[delta0]]$t1.fixed.error,
                        t1_interim_error = results[[delta0]]$t1.interim.error)
  results_df <- rbind(results_df, temp_df)
}

save(results_df,
     file = paste0(script_path,"/output/results_p00.4_k2_seed123456.R"))
