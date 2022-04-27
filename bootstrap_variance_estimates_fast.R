# Copyright 2021 Twitter, Inc.
# SPDX-License-Identifier: Apache-2.0

require(matrixStats)
require(ggplot2)
require(reshape2)
require(metafor)
require(dplyr)
require(xtable)

#function to generate data
generate_data <- function(scenario = 1, n = 5000, K = 100, replicates = 1000){
  #Y: a matrix of samples with number of rows = replicates and number of columns = K = number of groups
  #each column has "replicates" samples from group k with success probability mu_k
  #nks: a vector of length K with the kth entry giving the number of individuals in group k
  if(scenario ==1){
    #equal group size, equal accuracies 
    nks <- rep(n/K, K)
    mus <- rep(0.8, K)
    Y <- sapply(seq(K),function(j){rbinom(replicates,nks[j],mus[j])})
  }
  
  if(scenario == 2){
    #unequal group-size, equal accuracies
    nks <-round(seq(10, 90, length.out = K))
    mus <- rep(0.8, K)
    Y <- sapply(seq(K),function(j){rbinom(replicates,nks[j],mus[j])})
  }
  
  if(scenario == 3){
    #equal group size, unequal accuracies
    nks <- rep(n/K, K)
    mus <- seq(0.1, 0.9, length.out = K)
    Y <- sapply(seq(K),function(j){rbinom(replicates,nks[j],mus[j])})
  }
  
  if(scenario == 4){
    #unequal group-size, unequal accuracies
    nks <-round(seq(10, 90, length.out = K))
    mus <- seq(0.1, 0.9, length.out = K)
    Y <- sapply(seq(K),function(j){rbinom(replicates,nks[j],mus[j])})
  }
  
  Y <- t(Y)
  mu_hats <- Y/nks
  
  out <- list(Y = Y, nks= nks, mu_hats = mu_hats)
}


apply_corrections <- function(dat, threshold = FALSE){
  #applies three candidate corrections to the simulated binomial data
  uncorrected_var <- apply(dat$mu_hats, 2, var)
  sigma2_hats <- dat$mu_hats*(1-dat$mu_hats)/dat$nks
  corrected_var <- uncorrected_var - apply(sigma2_hats, 2, mean)
  double_corrected_var <-uncorrected_var - 2*apply(sigma2_hats, 2, mean) + apply(sigma2_hats/dat$nks, 2, mean)
  
  
  
  if(threshold){
    corrected_var <- pmax(0, corrected_var)
    double_corrected_var <- pmax(0, double_corrected_var)
  }
  
  out <- data.frame(uncorrected_var, corrected_var, double_corrected_var)
  
}


generate_bootstrap_samples <- function(dat, B, threshold){
  #generate B bootstrap samples for each replicate
  ## need to make sure this is working right
  bs_estimates <- array(0, dim = c(B, 3, ncol(dat$mu_hats)))
  
  
  #loop over replicates
  for(replicate in 1:ncol(dat$mu_hats)){
    if(replicate%%100 == 0){
      print(replicate)
    }
    #take a sample for the jth replicate, using the mu hats from within each group to re-sample
    
    mus.tmp <- dat$mu_hats[,replicate]
    samps <- t(sapply(seq(K),function(j){rbinom(B, dat$nks[j],mus.tmp[j])}))
    
    tmp <- list(Y = samps, nks = dat$nks, mu_hats = samps/dat$nks)
    estimates <- apply_corrections(tmp, threshold)
    bs_estimates[,,replicate] <- as.matrix(estimates)
    
  }
  return(bs_estimates)
}
  

DL_fun <- function(Y, ni){
  #metafor package
  
  Y[Y==0] <- 0.5
  Y[Y==ni] <- ni - 0.5
  
  mu_hats <- Y/ni
  vi <- mu_hats*(1-mu_hats)/ni

  res <- rma(mu_hats, vi, method="DL", digits=2)
  return(res$tau2)
}

non_parametric_bootstrap <- function(dat, B, threshold = FALSE){
  #non-parametric DL bootstrap
    #run this separately so the other one remains fast
    K <- length(dat$nks)
    non_parametric_matrix <- matrix(0, nrow = B, ncol = ncol(dat$Y)) #bootstrap samples by replicates, each one gets an estimate
    for(replicate in 1:ncol(dat$mu_hats)){
      
      if(replicate%%100 == 0){
        print(replicate)
      }
      
      Y.rep <- dat$Y[,replicate]
      inds <- matrix(sample(1:K, K*B, replace = TRUE), nrow = K, ncol = B)
      Y.samps <- sapply(seq(B), function(j){Y.rep[inds[,j]]})
      nk.samps <- sapply(seq(B), function(j){dat$nks[inds[,j]]})
      muhat.samps <- Y.samps/nk.samps
      tmp <- list(Y = Y.samps, nks = nk.samps, mu_hats = muhat.samps)
      
       #for(b in 1:B){
      #  inds <- sample(1:K, K, replace = TRUE)
       # DL_ests <- DL_fun(Y.rep[inds], dat$nks[inds])
      #  DL_matrix[b,replicate] <- DL_ests
      #}
      
      ests <- apply_corrections(tmp, threshold)
      non_parametric_matrix[,replicate] <- ests$corrected_var
    }
    
    return(non_parametric_matrix)
}



parametric_bootstrap <- function(dat, B, threshold = FALSE){
  #estimate it once for each replicate
  
  #look at thse results
  #Viechtbauer (2007)
  
  K <- length(dat$nks)
  parametric_matrix <- matrix(0, nrow = B, ncol = ncol(dat$Y))
  
  for(replicate in 1:ncol(dat$mu_hats)){
    
    if(replicate%%100 == 0){
      print(replicate)
    }
    
    Y <- dat$Y[,replicate]
    nks <- dat$nks
    
    Y[Y==0] <- 0.5
    Y[Y==nks] <- nks - 0.5
    
    mu_hats <- Y/nks
    vi <- mu_hats*(1-mu_hats)/nks
    
    res <- rma(mu_hats, vi, method="DL", digits=2)
    
   
      
    Y.boot <- sapply(seq(B), function(x){rnorm(K, res$b, sqrt(res$tau2 + vi))}) #this should be a sample of the observed data, i.e. the mu_hats

    tmp <- list(Y = Y.boot, nks = dat$nks, mu_hats = Y.boot)
    ests <- apply_corrections(tmp, threshold)

    #new.est <- rma(Y.boot, vi, method="DL", digits=2)
    parametric_matrix[, replicate] <- ests$corrected_var
    
    }

  return(parametric_matrix)
  
}


##----- Experiment 1: Verify correction works by simulation, no bootstrapping ---- ##
threshold <- TRUE
K <- 100; n <- 5000

scenario_names <- c("Equal Size; Equal Perf",
  "Unequal Size; Equal Perf",
  "Equal Size; Unequal Perf",
  "Unequal Size; Unequal Perf")

truth <- data.frame(scenario = scenario_names, truth = c(0,0, var(seq(0.1, 0.9, length.out = K)), var(seq(0.1, 0.9, length.out = K))))

results <- NULL
for(scenario in 1:length(scenario_names)){
  dat <- generate_data(scenario, n = n, K = K)
  estimates <- apply_corrections(dat, threshold = threshold) 
  estimates <- estimates %>% melt()
  estimates$scenario <- scenario_names[scenario]
  results <- rbind(results, estimates)
  
  dim(results)
  }

ggplot(results %>% filter(variable != "double_corrected_var"), aes(x = value, fill = variable)) + 
  geom_histogram(position = 'identity', alpha = .5) + 
  facet_wrap(vars(scenario), scales = 'free') + 
  geom_vline(data = truth, aes(xintercept = truth)) + 
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  xlab("estimate")

ggsave(paste("~/Desktop/Experiment1-", threshold,"-", K, "-", n, ".png", sep = ''))

##----- Experiment 2: What happens if we bootstrap? ----- ##
B <- 500 #number of bootstrap samples


results <- list()
for(scenario in 1:length(scenario_names)){
  
  dat <- generate_data(scenario, n = n, K = K)
  bs_estimates <- generate_bootstrap_samples(dat, B, threshold = threshold)
  results[[scenario_names[scenario]]] <- bs_estimates

}

pick_replicate <- 20
results_single <- lapply(results, function(x){
  out <- data.frame(x[,,pick_replicate])
  colnames(out) <- c("uncorrected_var", "corrected_var", "double_corrected_var")
  out <- out %>% melt()})

results_single <- do.call("rbind", results_single)
results_single$scenario <- rep(scenario_names, each = nrow(results_single)/length(scenario_names))

interval_single <- results_single %>% group_by(variable, scenario) %>% summarize(lower = quantile(value, prob=c(.025)), upper = quantile(value, prob = .975))
interval_single$y <- 20 + 10*as.numeric(as.factor(interval_single$variable)) 
interval_single <- interval_single %>% filter(variable != "double_corrected_var")

ggplot(results_single %>% filter(variable != "double_corrected_var"), aes(x = value, fill = variable)) + 
  geom_histogram(position = 'identity', alpha = .5)+
  geom_vline(data = truth, aes(xintercept = truth)) +
  geom_errorbarh(data = interval_single, mapping = aes(x = lower, y = y, xmin = lower, xmax = upper, color = variable), size = 1, height = 20) +
  facet_wrap(vars(scenario), scales = 'free_x') + 
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  xlab("estimate")

ggsave(paste("~/Desktop/Experiment2-", threshold,"-", K, "-", n, ".png", sep = ''))

##----- What if we use the double corrected estimator? ----- 
interval_single <- results_single %>% group_by(variable, scenario) %>% summarize(lower = quantile(value, prob=c(.025)), upper = quantile(value, prob = .975))
interval_single$y <- 20 + 20*as.numeric(as.factor(interval_single$variable)) 


ggplot(results_single, aes(x = value, fill = variable)) + 
  geom_histogram(position = 'identity', alpha = .5)+
  geom_vline(data = truth, aes(xintercept = truth)) +
  geom_errorbarh(data = interval_single, mapping = aes(x = lower, y = y, xmin = lower, xmax = upper, color = variable), size = 1, height = 20) +
  facet_wrap(vars(scenario), scales = 'free_x') +
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  xlab("estimate") +
  ylab("count")


ggsave(paste("~/Desktop/Experiment3-", threshold,"-", K, "-", n, ".png", sep = ''))


##----- what is the coverage? -----
intervals <- lapply(results, function(x){apply(x, c(2,3), quantile, probs = c(.025, .975))})
coverage <- lapply(seq(length(scenario_names)), function(j){apply(intervals[[j]][1,,] <= truth$truth[j] & intervals[[j]][2,,] >= truth$truth[j], 1, mean)})
coverage <- do.call("rbind", coverage)
colnames(coverage) <- c("uncorrected_var", "corrected_var", "double_corrected_var")
rownames(coverage) <- scenario_names


print(file = "~/Desktop/coverage_table.tex", xtable(coverage*100, digits = 1))

ci.length <- lapply(seq(length(scenario_names)), function(j){apply(intervals[[j]][2,,] - intervals[[j]][1,,], 1, mean)})
ci.length <- do.call("rbind", ci.length)
colnames(ci.length) <- c("uncorrected_var", "corrected_var", "double_corrected_var")
rownames(ci.length) <- scenario_names



##----- bootstrap medians?
replicates <- 1000
medians <- lapply(results, function(x){z <- data.frame(t(apply(x, c(2,3), median)))
colnames(z) <- c("uncorrected_var", "corrected_var", "double_corrected_var")
z <- z %>% melt()
})
median_dist <- data.frame(do.call("rbind", medians))
median_dist$scenario <- rep(scenario_names, each = replicates*3)

  
ggplot(median_dist, aes(x = value, fill = variable)) + 
  geom_histogram(position = 'identity', alpha = .5)+
  geom_vline(data = truth, aes(xintercept = truth)) +
  facet_wrap(vars(scenario), scales = 'free') 

ggsave(paste("~/Desktop/Experiment4-", threshold,"-", K, "-", n, ".png", sep = ''))


coverage
round(ci.length, 3)

##----- other approaches to bootstrapping -----##
parametric.results <- list()
nonparametric.results <- list()
for(scenario in 1:length(scenario_names)){
  
  dat <- generate_data(scenario)
  #try DL bootstrap
  
 non_parametric_bs_ests <- non_parametric_bootstrap(dat, B) 
 parametric_bs_ests <- parametric_bootstrap(dat, B)  #fix this one up to be HO
 
 parametric.results[[scenario]] <- parametric_bs_ests
 nonparametric.results[[scenario]] <- non_parametric_bs_ests
 
}

parametric.intervals <- lapply(parametric.results, function(x){apply(x, 2, quantile, probs = c(.025, .975))})
parametric.coverage <- lapply(seq(length(scenario_names)), function(j){mean(parametric.intervals[[j]][1,] <= truth$truth[j] & parametric.intervals[[j]][2,] >= truth$truth[j])})

nonparametric.intervals <- lapply(nonparametric.results, function(x){apply(x, 2, quantile, probs = c(.025, .975))})
nonparametric.coverage <- lapply(seq(length(scenario_names)), function(j){mean(nonparametric.intervals[[j]][1,] <= truth$truth[j] & nonparametric.intervals[[j]][2,] >= truth$truth[j])})

nonparametric.coverage
parametric.coverage

coverage <- data.frame(coverage)
coverage$`non-parametric` <- unlist(nonparametric.coverage)
coverage$parametric <- unlist(parametric.coverage)

print(xtable(coverage), file = paste("~/Desktop/Coverage-Table-", threshold,"-", K, "-", n, ".png", sep = ''))



diffs <- NULL
mu_k <- seq(.5, .9, .01)

for(i in 1:100){
Y_k <- rnorm(length(mu_k), mu_k, .05)

diffs <- c(diffs, diff(range(Y_k)))
}

diff(range(mu_k))
mean(diffs)

