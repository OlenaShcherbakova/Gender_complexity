#extracting marginal likelihoods from dependent and independent models to later calculate Bayes Factors

source('library.R')

stones.dep.world <- 
  list.files(path = "./results_bayestraits/", pattern = "\\.dependent.log.Stones.txt", 
             full.names = TRUE, recursive = TRUE)

stones.indep.world <- 
  list.files(path = "./results_bayestraits/", pattern = "\\.independent.log.Stones.txt", 
             full.names = TRUE, recursive = TRUE)

stones_files_dep_world <- lapply(stones.dep.world, bt_read.stones)
stones_files_indep_world <- lapply(stones.indep.world, bt_read.stones)



dependent_Lh <- 0

for (i in 1:length(stones_files_dep_world)) {
  dependent_Lh[i] <- stones_files_dep_world[[i]][["marginal_likelihood"]]
}

independent_Lh <- 0

for (i in 1:length(stones_files_indep_world)) {
  independent_Lh[i] <- stones_files_indep_world[[i]][["marginal_likelihood"]]
}

likelihoods <- data.frame(independent_Lh, dependent_Lh)

colnames(likelihoods) <- c("Lh_indep", "Lh_dep")