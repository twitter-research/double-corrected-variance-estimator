# Copyright 2021 Twitter, Inc.
# SPDX-License-Identifier: Apache-2.0

require(dplyr)
require(reshape2)
require(ggplot2)

mean_pairwise_abs_diff <- function(x){
  out <- mean(as.numeric(dist(x, method = "manhattan")))
  return(out)
}

gini <- function(x){
  out <- mean(as.numeric(dist(x, method = "manhattan")))/2*mean(x)
  return(out)
}

mean_abs_diff <- function(x){
  out <- mean(abs(x - mean(x)))
  return(out)
}

maxmin_diff <- function(x){
  return(max(x) - min(x))
}

maxmin_ratio <- function(x, eps = .05){
  return(max(x + eps)/min(x + eps))
}


max_abs_diff <- function(x){
  return(max(abs(x - mean(x))))
}

generalized_entropy <- function(x, alpha = 1){
  K <- length(x)
  if(alpha != 0 & alpha != 1){
    out <- 1/(K*alpha*(alpha - 1))*sum((x/mean(x))^alpha-1)
  }
  if(alpha == 1){
    out <- 1/K*sum(x/mean(x)*log(x/mean(x)))
  }
  if(alpha == 0){
    out <- -1/K*sum(log(x/mean(x)))
  }
  return(out)
}


Ks <- seq(5, 150, 5)
lowers <- seq(.1, .9, .1)
n <- 5000

results <- list()

for(lower in lowers){
  biases <- NULL
  
  for(K in Ks){
    mu_k <- seq(lower, .9, length.out = K)
    nks <- rep(round(n/K), K)
    mmd <- NULL
    mmr <- NULL
    mads <- NULL
    pw_abs <- NULL
    max_ads <- NULL
    vars <- NULL
    geis <- NULL
    ginis <- NULL
    
    for(i in 1:1000){
      Y_k <- rbinom(length(mu_k), nks,  mu_k)/nks
      
      mmd <- c(mmd, maxmin_diff(Y_k))
      mmr <- c(mmr, maxmin_ratio(Y_k))
      mads <- c(mads, mean_abs_diff(Y_k))
      pw_abs <- c(pw_abs,  mean_pairwise_abs_diff(Y_k))
      vars <- c(vars, var(Y_k))
      max_ads <- c(max_ads, max_abs_diff(Y_k))
      geis <- c(geis, generalized_entropy(Y_k, alpha = .5))
      ginis <- c(ginis, gini(Y_k))
    }
    
    bias <-  c(mean(mmd) - maxmin_diff(mu_k), 
               mean(mmr) - maxmin_ratio(mu_k), 
               mean(mads) - mean_abs_diff(mu_k), 
               mean(pw_abs) - mean_pairwise_abs_diff(mu_k), 
               mean(vars) - var(mu_k), 
               mean(max_ads) - max_abs_diff(mu_k),
               mean(geis - generalized_entropy(mu_k, alpha = .5)), 
               mean(ginis) - gini(mu_k))
    
    biases <- rbind(biases, bias)
  }
  
  results[[as.character(lower)]] <- biases
}



all.biases <- do.call("rbind", results)
allKs <- rep(Ks, length(lowers))
allLowers <- rep(lowers, each = length(Ks))

df <- data.frame(maxmin_diff = all.biases[,1],  max_min_ratio = all.biases[,2], 
                 max_abs_diff = all.biases[,6],  mean_abs_dev = all.biases[,3],  var = all.biases[,5], gei = all.biases[,7], K = allKs, lower = allLowers) #

df <- df %>% melt(id.vars=c("K", "lower"))

ggplot(df, aes(x = K, y = value, color = factor(lower), group = factor(lower))) + 
  geom_line() + 
  theme_bw() + 
  ggtitle("Statistical Bias of Model ``Bias'' Estimators") +
  theme(text = element_text(size = 25)) +
  guides(color=guide_legend(title="Lower Bound")) + 
  facet_wrap(vars(variable), scales = 'free') +
  ylab(expression(paste("E[M(Y)] - E[M(", mu, ")]", sep='')))

ggsave("~/Desktop/max-min.png", width = 12, height = 7, units = 'in')
