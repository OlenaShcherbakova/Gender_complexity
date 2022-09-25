#distribution of r (Bayesian correlation coefficients) values: median; Lower and Upper intervals

source('library.R')

results.dep.world <- 
  list.files(path = "./results_bayestraits_full/", pattern = "\\.dependent.log.Log.txt", 
             full.names = TRUE, recursive = TRUE)

results.indep.world <- 
  list.files(path = "./results_bayestraits_full/", pattern = "\\.independent.log.Log.txt", 
             full.names = TRUE, recursive = TRUE)

log_files_dep <- lapply(results.dep.world, bt_read.log)
log_files_indep <- lapply(results.indep.world, bt_read.log)

desired_length <- length(results.dep.world)
correlation_dep <- vector(mode = "list", length = desired_length)

for (i in 1:length(correlation_dep[])){
  correlation_dep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(correlation_dep[])) {
  correlation_dep[[i]] <- log_files_dep[[i]][["R Trait 1 2"]][-c(1:999)]
}

correlation_indep <- vector(mode = "list", length = desired_length)

for (i in 1:length(correlation_indep[])){
  correlation_indep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(correlation_indep[])) {
  correlation_indep[[i]] <- log_files_indep[[i]][["R Trait 1 2"]][-c(1:999)]
}


correlation_world_dependent <- data.frame(correlation_dep)
correlation_world_dependent <- t(correlation_world_dependent)

median <- data.frame(apply(correlation_world_dependent, 1, median))
names(median)[1] <- "median_correlation_dep"
correlation_world_dependent <- cbind(correlation_world_dependent, median)

intervals <- apply(correlation_world_dependent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_correlation_dep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_correlation_dep"

correlation_world_dependent <- cbind(correlation_world_dependent, lower)
correlation_world_dependent <- cbind(correlation_world_dependent, upper)




correlation_world_independent <- data.frame(correlation_indep)
correlation_world_independent <- t(correlation_world_independent)

intervals <- apply(correlation_world_independent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_correlation_indep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_correlation_indep"

correlation_world_independent <- cbind(correlation_world_independent, lower)
correlation_world_independent <- cbind(correlation_world_independent, upper)

correlations <- cbind(correlation_world_dependent[c("median_correlation_dep", "lower_correlation_dep", "upper_correlation_dep")])