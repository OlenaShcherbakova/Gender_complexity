#BT files generator
#interdependence between semantic classes and agreement patterns
#code modified from the code originally written by Simon J. Greenhill

library(here)
source(here('library.R'))

ANALYSES <- list(
  "Semantic_classes_Agreement_patterns"= c("sem_classes_10", "agr_patterns_10")
)

gb <- load_data_final_short() 

gb$sem_classes <- as.numeric(as.character(gb$sem_classes))
gb$agr_patterns <- as.numeric(as.character(gb$agr_patterns))

gb <- gb %>%
  filter(sem_classes > 0 | agr_patterns > 0) %>%
  mutate(sem_classes_10 = sem_classes * 10,
         agr_patterns_10 = agr_patterns * 10) %>%
  dplyr::select(Glottocode, sem_classes_10, agr_patterns_10)

COMMAND_DEPENDENT <- "4
2
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.dependent.log
Run
"

COMMAND_INDEPENDENT <- "4
2
TestCorrel
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.independent.log
Run
"

# the template for the cluster scheduler / job runner
SLURM_TEMPLATE <- "#!/bin/sh
#SBATCH --job-name bayestraits
#SBATCH --get-user-env
#SBATCH --output run.%%j.log
#SBATCH --mem 4G
#SBATCH --cpus-per-task 1
BayesTraitsV4 %s %s < %s
"

#phylogenies with respective taxa files
PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)



#Indo-European (Bouckaert et al. 2012)
basedir <- basename(PHYLOGENIES[1])
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[1]) #"data/phylogenies/bouckaert_et_al2012/"

taxa <- read.csv(TAXAS[1])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "sem_classes_10", "agr_patterns_10"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))






#Austronesian (Gray et al. 2009)
basedir <- basename(PHYLOGENIES[2])
basedir <- paste0(basedir, "_reduced", paste="")
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[2]) #"data/phylogenies/gray_et_al2009/"

taxa <- read.csv(TAXAS[2])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "sem_classes_10", "agr_patterns_10"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))





#Bantu (Grollemund et al. 2015)
basedir <- basename(PHYLOGENIES[3])
basedir <- paste0(basedir, "_reduced", paste="")
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[3]) #"data/phylogenies/grollemund_et_al2015/"

taxa <- read.csv(TAXAS[3])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "sem_classes_10", "agr_patterns_10"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))




#Dravidian (Kolipakam et al. 2018)
basedir <- basename(PHYLOGENIES[4])
basedir <- paste0(basedir, "_reduced", paste="")
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[4]) #"data/phylogenies/kolipakam_et_al2018/"

taxa <- read.csv(TAXAS[4])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "sem_classes_10", "agr_patterns_10"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))






#World/global tree (Jäger et al. 2018)
basedir <- basename(PHYLOGENIES[5]) 
basedir <- paste0(basedir, "_reduced", paste="")
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[5]) #"data/phylogenies/world/"

taxa <- read.csv(TAXAS[5])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "sem_classes_10", "agr_patterns_10"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))
