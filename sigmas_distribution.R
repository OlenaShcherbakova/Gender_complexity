#distribution of sigma values: median; Lower and Upper intervals

source('library.R')

results.dep.world <- 
  list.files(path = "./results_bayestraits/", pattern = "\\.dependent.log.Log.txt", 
             full.names = TRUE, recursive = TRUE)

results.indep.world <- 
  list.files(path = "./results_bayestraits/", pattern = "\\.independent.log.Log.txt", 
             full.names = TRUE, recursive = TRUE)

log_files_dep <- lapply(results.dep.world, bt_read.log)
log_files_indep <- lapply(results.indep.world, bt_read.log)

desired_length <- length(results.dep.world)
sigma1_dep <- vector(mode = "list", length = desired_length)

for (i in 1:length(sigma1_dep[])){
  sigma1_dep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(sigma1_dep[])) {
  sigma1_dep[[i]] <- log_files_dep[[i]][["Sigma^2 1"]][-c(1:999)]
}

sigma1_indep <- vector(mode = "list", length = desired_length)

for (i in 1:length(sigma1_indep[])){
  sigma1_indep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(sigma1_indep[])) {
  sigma1_indep[[i]] <- log_files_indep[[i]][["Sigma^2 1"]][-c(1:999)]
}


sigma1_world_dependent <- data.frame(sigma1_dep)
sigma1_world_dependent <- t(sigma1_world_dependent)

median <- data.frame(apply(sigma1_world_dependent, 1, median))
names(median)[1] <- "median_sigma1_dep"
sigma1_world_dependent <- cbind(sigma1_world_dependent, median)

intervals <- apply(sigma1_world_dependent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_sigma1_dep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_sigma1_dep"

sigma1_world_dependent <- cbind(sigma1_world_dependent, lower)
sigma1_world_dependent <- cbind(sigma1_world_dependent, upper)




sigma1_world_independent <- data.frame(sigma1_indep)
sigma1_world_independent <- t(sigma1_world_independent)

median <- data.frame(apply(sigma1_world_independent, 1, median))
names(median)[1] <- "median_sigma1_indep"
sigma1_world_independent <- cbind(sigma1_world_independent, median)

intervals <- apply(sigma1_world_independent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_sigma1_indep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_sigma1_indep"

sigma1_world_independent <- cbind(sigma1_world_independent, lower)
sigma1_world_independent <- cbind(sigma1_world_independent, upper)

sigma1s <- cbind(sigma1_world_dependent[c("median_sigma1_dep", "lower_sigma1_dep", "upper_sigma1_dep")], sigma1_world_independent[c("median_sigma1_indep", "lower_sigma1_indep", "upper_sigma1_indep")])














desired_length <- length(results.dep.world)
sigma2_dep <- vector(mode = "list", length = desired_length)

for (i in 1:length(sigma2_dep[])){
  sigma2_dep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(sigma2_dep[])) {
  sigma2_dep[[i]] <- log_files_dep[[i]][["Sigma^2 2"]][-c(1:999)]
}

sigma2_indep <- vector(mode = "list", length = desired_length)

for (i in 1:length(sigma2_indep[])){
  sigma2_indep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(sigma2_indep[])) {
  sigma2_indep[[i]] <- log_files_indep[[i]][["Sigma^2 2"]][-c(1:999)]
}


sigma2_world_dependent <- data.frame(sigma2_dep)
sigma2_world_dependent <- t(sigma2_world_dependent)

median <- data.frame(apply(sigma2_world_dependent, 1, median))
names(median)[1] <- "median_sigma2_dep"
sigma2_world_dependent <- cbind(sigma2_world_dependent, median)

intervals <- apply(sigma2_world_dependent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_sigma2_dep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_sigma2_dep"

sigma2_world_dependent <- cbind(sigma2_world_dependent, lower)
sigma2_world_dependent <- cbind(sigma2_world_dependent, upper)




sigma2_world_independent <- data.frame(sigma2_indep)
sigma2_world_independent <- t(sigma2_world_independent)

median <- data.frame(apply(sigma2_world_independent, 1, median))
names(median)[1] <- "median_sigma2_indep"
sigma2_world_independent <- cbind(sigma2_world_independent, median)

intervals <- apply(sigma2_world_independent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_sigma2_indep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_sigma2_indep"

sigma2_world_independent <- cbind(sigma2_world_independent, lower)
sigma2_world_independent <- cbind(sigma2_world_independent, upper)

sigma2s <- cbind(sigma2_world_dependent[c("median_sigma2_dep", "lower_sigma2_dep", "upper_sigma2_dep")], sigma2_world_independent [c("median_sigma2_indep", "lower_sigma2_indep", "upper_sigma2_indep")])

sigmas <- cbind(sigma1s, sigma2s)