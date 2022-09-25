#distribution of alpha values: median; Lower and Upper intervals

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
alpha1_dep <- vector(mode = "list", length = desired_length)

for (i in 1:length(alpha1_dep[])){
  alpha1_dep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(alpha1_dep[])) {
  alpha1_dep[[i]] <- log_files_dep[[i]][["Alpha 1"]][-c(1:999)]
}

alpha1_indep <- vector(mode = "list", length = desired_length)

for (i in 1:length(alpha1_indep[])){
  alpha1_indep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(alpha1_indep[])) {
  alpha1_indep[[i]] <- log_files_indep[[i]][["Alpha 1"]]
}

alpha1_world_dependent <- as.data.frame(alpha1_dep)
alpha1_world_dependent <- t(alpha1_world_dependent)

median <- data.frame(apply(alpha1_world_dependent, 1, median))
names(median)[1] <- "median_alpha1_dep"
alpha1_world_dependent <- cbind(alpha1_world_dependent, median)

intervals <- apply(alpha1_world_dependent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_alpha1_dep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_alpha1_dep"

alpha1_world_dependent <- cbind(alpha1_world_dependent, lower)
alpha1_world_dependent <- cbind(alpha1_world_dependent, upper)




alpha1_world_independent <- data.frame(alpha1_indep)
alpha1_world_independent <- t(alpha1_world_independent)

median <- data.frame(apply(alpha1_world_independent, 1, median))
names(median)[1] <- "median_alpha1_indep"
alpha1_world_independent <- cbind(alpha1_world_independent, median)

intervals <- apply(alpha1_world_independent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_alpha1_indep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_alpha1_indep"

alpha1_world_independent <- cbind(alpha1_world_independent, lower)
alpha1_world_independent <- cbind(alpha1_world_independent, upper)

alpha1s <- cbind(alpha1_world_dependent[c("median_alpha1_dep", "lower_alpha1_dep", "upper_alpha1_dep")], alpha1_world_independent[c("median_alpha1_indep", "lower_alpha1_indep", "upper_alpha1_indep")])







desired_length <- length(results.dep.world)
alpha2_dep <- vector(mode = "list", length = desired_length)

for (i in 1:length(alpha2_dep[])){
  alpha2_dep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(alpha2_dep[])) {
  alpha2_dep[[i]] <- log_files_dep[[i]][["Alpha 2"]]
}

alpha2_indep <- vector(mode = "list", length = desired_length)

for (i in 1:length(alpha2_indep[])){
  alpha2_indep[[i]]<- vector(mode='numeric', length = 9000)
}

for (i in 1:length(alpha2_indep[])) {
  alpha2_indep[[i]] <- log_files_indep[[i]][["Alpha 2"]]
}


alpha2_world_dependent <- data.frame(alpha2_dep)
alpha2_world_dependent <- t(alpha2_world_dependent)

median <- data.frame(apply(alpha2_world_dependent, 1, median))
names(median)[1] <- "median_alpha2_dep"
alpha2_world_dependent <- cbind(alpha2_world_dependent, median)

intervals <- apply(alpha2_world_dependent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_alpha2_dep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_alpha2_dep"

alpha2_world_dependent <- cbind(alpha2_world_dependent, lower)
alpha2_world_dependent <- cbind(alpha2_world_dependent, upper)




alpha2_world_independent <- data.frame(alpha2_indep)
alpha2_world_independent <- t(alpha2_world_independent)

median <- data.frame(apply(alpha2_world_independent, 1, median))
names(median)[1] <- "median_alpha2_indep"
alpha2_world_independent <- cbind(alpha2_world_independent, median)

intervals <- apply(alpha2_world_independent, 1, HPDI, prob=0.95)

lower <- data.frame(intervals[1,])
names(lower)[1] <- "lower_alpha2_indep"
upper <- data.frame(intervals[2,])
names(upper)[1] <- "upper_alpha2_indep"

alpha2_world_independent <- cbind(alpha2_world_independent, lower)
alpha2_world_independent <- cbind(alpha2_world_independent, upper)

alpha2s <- cbind(alpha2_world_dependent[c("median_alpha2_dep", "lower_alpha2_dep", "upper_alpha2_dep")], alpha2_world_independent [c("median_alpha2_indep", "lower_alpha2_indep", "upper_alpha2_indep")])

alphas <- cbind(alpha1s, alpha2s)
