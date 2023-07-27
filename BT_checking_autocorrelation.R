#checking for autocorrelation

results.dep.world <- 
  list.files(path = "./results_bayestraits/", pattern = "\\.dependent.log.Log.txt", 
             full.names = TRUE, recursive = TRUE)

results.indep.world <- 
  list.files(path = "./results_bayestraits/", pattern = "\\.independent.log.Log.txt", 
             full.names = TRUE, recursive = TRUE)

log_files_dep <- lapply(results.dep.world, bt_read.log)
log_files_indep <- lapply(results.indep.world, bt_read.log)

desired_length <- length(log_files_dep)
autocorrelation_dep <- vector(mode = "list", length = desired_length)
autocorrelation_indep <- vector(mode = "list", length = desired_length)


for (i in 1:length(autocorrelation_dep[])) {
  autocorrelation_dep[[i]] <- cor(log_files_dep[[i]]$Iteration, log_files_dep[[i]]$Lh)
}

for (i in 1:length(autocorrelation_indep[])) {
  autocorrelation_indep[[i]] <- cor(log_files_indep[[i]]$Iteration, log_files_indep[[i]]$Lh)
}

autocorrelation_dep <- data.frame(autocorrelation_dep)
autocorrelation_dep <- t(autocorrelation_dep)

autocorrelation_indep <- data.frame(autocorrelation_indep)
autocorrelation_indep <- t(autocorrelation_indep)

autocorrelation <- data.frame(cbind(autocorrelation_dep, autocorrelation_indep))
colnames(autocorrelation) <- c("autocorrelation_dependent", "autocorrelation_independent")
rownames(autocorrelation) <- NULL

txt <- results.dep.world

grep(".*Part Description:(.*)Installs.*", "\\1", txt)

desisered_length <- nrow(autocorrelation)
vector <- rep(c(NA), times=desired_length)

for (i in 1:length(txt)) {
  vector[i] <- str_split(txt, "/")[[i]][4]
}

#providing a column corresponding to the file name (First trait used in the analysis, second trait used in the analysis, and Phylogeny)

autocorrelation <- autocorrelation %>%
  mutate(Phylogeny = vector) %>%
  mutate(Phylogeny=recode(Phylogeny,
                          "gray_et_al2009" = "Austronesian",
                          "bouckaert_et_al2012" = "Indo-European",
                          "kolipakam_et_al2018" = "Dravidian",
                          "grollemund_et_al2015" = "Bantu")) %>% View()

