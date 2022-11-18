#measuring phylogenetic signal

addTaskCallback(function(...) {set.seed(123);TRUE}) #adding seed for reproducibility

#global tree
source("library.R")

gb <- load_data_final_short()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/world/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb %>%
  dplyr::select(Glottocode, sem_classes, agr_patterns, taxon) %>%
  filter(!is.na(taxon)) -> gb.subset

gb.subset$sem_classes <- as.numeric(as.character(gb.subset$sem_classes))
gb.subset$agr_patterns <- as.numeric(as.character(gb.subset$agr_patterns))

gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/world/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon

#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)


#measuring phylogenetic signal (global):
tree.sem_classes <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$sem_classes), 'Glottocode']))
sem_classes <- get_trait_vector(tree.sem_classes, gb_geo.subset, 'sem_classes')
physig_sem_classes_world_l <- phytools::phylosig(tree.sem_classes, sem_classes, method="lambda", test=TRUE)
lambda_sem_classes_world_l <- physig_sem_classes_world_l[1][["lambda"]]
LR_sem_classes_world_l <- 2*(physig_sem_classes_world_l$logL-physig_sem_classes_world_l$logL0) #performing likelihood ratio test
P_lambda_sem_classes_world_l <- physig_sem_classes_world_l$P


tree.agr_patterns <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$agr_patterns), 'Glottocode']))
agr_patterns <- get_trait_vector(tree.agr_patterns, gb_geo.subset, 'agr_patterns')
physig_agr_patterns_world_l <- phytools::phylosig(tree.agr_patterns, agr_patterns, method="lambda", test=TRUE)
lambda_agr_patterns_world_l <- physig_agr_patterns_world_l[1][["lambda"]]
LR_agr_patterns_world_l <- 2*(physig_agr_patterns_world_l$logL-physig_agr_patterns_world_l$logL0) #performing likelihood-ratio test
P_lambda_agr_patterns_world_l <- physig_agr_patterns_world_l$P

sem_classes_world <- c(lambda_sem_classes_world_l, physig_sem_classes_world_l$logL, physig_sem_classes_world_l$logL0, LR_sem_classes_world_l, P_lambda_sem_classes_world_l)
agr_patterns_world <- c(lambda_agr_patterns_world_l, physig_agr_patterns_world_l$logL, physig_agr_patterns_world_l$logL0, LR_agr_patterns_world_l, P_lambda_agr_patterns_world_l)



#Indo-European tree (1): Bouckaert et al 2012
source("library.R")
gb <- load_data_final_short()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/bouckaert_et_al2012/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb %>%
  #dplyr::select(Glottocode, sem_classes, agr_patterns, taxon) %>%
  filter(!is.na(taxon)) -> gb.subset

gb.subset$sem_classes <- as.numeric(as.character(gb.subset$sem_classes))
gb.subset$agr_patterns <- as.numeric(as.character(gb.subset$agr_patterns))

gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/bouckaert_et_al2012/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon

#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)


#measuring phylogenetic signal (Indo-European):
tree.sem_classes <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$sem_classes), 'Glottocode']))
sem_classes <- get_trait_vector(tree.sem_classes, gb_geo.subset, 'sem_classes')
physig_sem_classes_ie_l <- phytools::phylosig(tree.sem_classes, sem_classes, method="lambda", test=TRUE)
lambda_sem_classes_ie_l <- physig_sem_classes_ie_l[1][["lambda"]]
LR_sem_classes_ie_l <- 2*(physig_sem_classes_ie_l$logL-physig_sem_classes_ie_l$logL0) #performing likelihood ratio test
P_lambda_sem_classes_ie_l <- physig_sem_classes_ie_l$P



tree.agr_patterns <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$agr_patterns), 'Glottocode']))
agr_patterns <- get_trait_vector(tree.agr_patterns, gb_geo.subset, 'agr_patterns')
physig_agr_patterns_ie_l <- phytools::phylosig(tree.agr_patterns, agr_patterns, method="lambda", test=TRUE)
lambda_agr_patterns_ie_l <- physig_agr_patterns_ie_l[1][["lambda"]]
LR_agr_patterns_ie_l <- 2*(physig_agr_patterns_ie_l$logL-physig_agr_patterns_ie_l$logL0) #performing likelihood-ratio test
P_lambda_agr_patterns_ie_l <- physig_agr_patterns_ie_l$P


sem_classes_ie <- c(lambda_sem_classes_ie_l, physig_sem_classes_ie_l$logL, physig_sem_classes_ie_l$logL0, LR_sem_classes_ie_l, P_lambda_sem_classes_ie_l)
agr_patterns_ie <- c(lambda_agr_patterns_ie_l, physig_agr_patterns_ie_l$logL, physig_agr_patterns_ie_l$logL0, LR_agr_patterns_ie_l, P_lambda_agr_patterns_ie_l)


#Bantu tree
source("library.R")
gb <- load_data_final_short()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/grollemund_et_al2015/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb %>%
  dplyr::select(Glottocode, sem_classes, agr_patterns, taxon) %>%
  filter(!is.na(taxon)) -> gb.subset

gb.subset$sem_classes <- as.numeric(as.character(gb.subset$sem_classes))
gb.subset$agr_patterns <- as.numeric(as.character(gb.subset$agr_patterns))

gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/grollemund_et_al2015/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon

#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)


#measuring phylogenetic signal (Bantu):
tree.sem_classes <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$sem_classes), 'Glottocode']))
sem_classes <- get_trait_vector(tree.sem_classes, gb_geo.subset, 'sem_classes')
physig_sem_classes_b_l <- phytools::phylosig(tree.sem_classes, sem_classes, method="lambda", test=TRUE)
lambda_sem_classes_b_l <- physig_sem_classes_b_l[1][["lambda"]]
LR_sem_classes_b_l <- 2*(physig_sem_classes_b_l$logL-physig_sem_classes_b_l$logL0) #performing likelihood ratio test
P_lambda_sem_classes_b_l <- physig_sem_classes_b_l$P



tree.agr_patterns <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$agr_patterns), 'Glottocode']))
agr_patterns <- get_trait_vector(tree.agr_patterns, gb_geo.subset, 'agr_patterns')
physig_agr_patterns_b_l <- phytools::phylosig(tree.agr_patterns, agr_patterns, method="lambda", test=TRUE)
lambda_agr_patterns_b_l <- physig_agr_patterns_b_l[1][["lambda"]]
LR_agr_patterns_b_l <- 2*(physig_agr_patterns_b_l$logL-physig_agr_patterns_b_l$logL0) #performing likelihood-ratio test
P_lambda_agr_patterns_b_l <- physig_agr_patterns_b_l$P


sem_classes_b <- c(lambda_sem_classes_b_l, physig_sem_classes_b_l$logL, physig_sem_classes_b_l$logL0, LR_sem_classes_b_l, P_lambda_sem_classes_b_l)
agr_patterns_b <- c(lambda_agr_patterns_b_l, physig_agr_patterns_b_l$logL, physig_agr_patterns_b_l$logL0, LR_agr_patterns_b_l, P_lambda_agr_patterns_b_l)



#Dravidian tree
source("library.R")

gb <- load_data_final_short()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/kolipakam_et_al2018/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb %>%
  dplyr::select(Glottocode, sem_classes, agr_patterns, taxon) %>%
  filter(!is.na(taxon)) -> gb.subset

gb.subset$sem_classes <- as.numeric(as.character(gb.subset$sem_classes))
gb.subset$agr_patterns <- as.numeric(as.character(gb.subset$agr_patterns))

gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/kolipakam_et_al2018/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon

#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)


#measuring phylogenetic signal (Dravidian):
tree.sem_classes <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$sem_classes), 'Glottocode']))
sem_classes <- get_trait_vector(tree.sem_classes, gb_geo.subset, 'sem_classes')
physig_sem_classes_dr_l <- phytools::phylosig(tree.sem_classes, sem_classes, method="lambda", test=TRUE)
lambda_sem_classes_dr_l <- physig_sem_classes_dr_l[1][["lambda"]]
LR_sem_classes_dr_l <- 2*(physig_sem_classes_dr_l$logL-physig_sem_classes_dr_l$logL0) #performing likelihood ratio test
P_lambda_sem_classes_dr_l <- physig_sem_classes_dr_l$P


tree.agr_patterns <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$agr_patterns), 'Glottocode']))
agr_patterns <- get_trait_vector(tree.agr_patterns, gb_geo.subset, 'agr_patterns')
physig_agr_patterns_dr_l <- phytools::phylosig(tree.agr_patterns, agr_patterns, method="lambda", test=TRUE)
lambda_agr_patterns_dr_l <- physig_agr_patterns_dr_l[1][["lambda"]]
LR_agr_patterns_dr_l <- 2*(physig_agr_patterns_dr_l$logL-physig_agr_patterns_dr_l$logL0) #performing likelihood-ratio test
P_lambda_agr_patterns_dr_l <- physig_agr_patterns_dr_l$P


sem_classes_dr <- c(lambda_sem_classes_dr_l, physig_sem_classes_dr_l$logL, physig_sem_classes_dr_l$logL0, LR_sem_classes_dr_l, P_lambda_sem_classes_dr_l)
agr_patterns_dr <- c(lambda_agr_patterns_dr_l, physig_agr_patterns_dr_l$logL, physig_agr_patterns_dr_l$logL0, LR_agr_patterns_dr_l, P_lambda_agr_patterns_dr_l)



#Austronesian tree
source("library.R")

gb <- load_data_final_short()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/gray_et_al2009/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb %>%
  dplyr::select(Glottocode, sem_classes, agr_patterns, taxon) %>%
  filter(!is.na(taxon)) -> gb.subset

gb.subset$sem_classes <- as.numeric(as.character(gb.subset$sem_classes))
gb.subset$agr_patterns <- as.numeric(as.character(gb.subset$agr_patterns))

gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/gray_et_al2009/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon

#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)


#measuring phylogenetic signal (Austronesian):
tree.sem_classes <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$sem_classes), 'Glottocode']))
sem_classes <- get_trait_vector(tree.sem_classes, gb_geo.subset, 'sem_classes')
physig_sem_classes_a_l <- phytools::phylosig(tree.sem_classes, sem_classes, method="lambda", test=TRUE)
lambda_sem_classes_a_l <- physig_sem_classes_a_l[1][["lambda"]]
LR_sem_classes_a_l <- 2*(physig_sem_classes_a_l$logL-physig_sem_classes_a_l$logL0) #performing likelihood ratio test
P_lambda_sem_classes_a_l <- physig_sem_classes_a_l$P


tree.agr_patterns <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$agr_patterns), 'Glottocode']))
agr_patterns <- get_trait_vector(tree.agr_patterns, gb_geo.subset, 'agr_patterns')
physig_agr_patterns_a_l <- phytools::phylosig(tree.agr_patterns, agr_patterns, method="lambda", test=TRUE)
lambda_agr_patterns_a_l <- physig_agr_patterns_a_l[1][["lambda"]]
LR_agr_patterns_a_l <- 2*(physig_agr_patterns_a_l$logL-physig_agr_patterns_a_l$logL0) #performing likelihood-ratio test
P_lambda_agr_patterns_a_l <- physig_agr_patterns_a_l$P

sem_classes_a <- c(lambda_sem_classes_a_l, physig_sem_classes_a_l$logL, physig_sem_classes_a_l$logL0, LR_sem_classes_a_l, P_lambda_sem_classes_a_l)
agr_patterns_a <- c(lambda_agr_patterns_a_l, physig_agr_patterns_a_l$logL, physig_agr_patterns_a_l$logL0, LR_agr_patterns_a_l, P_lambda_agr_patterns_a_l)



#Making a table out of two measures of phylogenetic signal
physig <- as.data.frame(rbind(sem_classes_a, agr_patterns_a, sem_classes_b, agr_patterns_b, sem_classes_dr, agr_patterns_dr, sem_classes_ie, agr_patterns_ie, sem_classes_world, agr_patterns_world))
colnames(physig) <- c("lambda", "logL", "logL0", "LR (lambda)", "p-value (lambda)")
physig <- round(physig, digits=2)
#rownames(physig) <- c("Semantic rules (Austronesian)", "Agreement patterns (Austronesian)", "Semantic rules (Bantu)", "Agreement patterns (Bantu)", "Semantic rules (Dravidian)", "Agreement patterns (Dravidian)", "Semantic rules (Indo-European)", "Agreement patterns (Indo-European)", "Semantic rules (World)", "Agreement patterns (World)")
phylogeny <- as.data.frame(c(rep(c("Austronesian"), times=2), rep(c("Bantu"), times=2), rep(c("Dravidian"), times=2), rep(c("Indo-European"), times=2), rep(c("World"), times=2)))
colnames(phylogeny) <- "Phylogeny"
features <- as.data.frame(c(rep(c("Semantic rules", "Agreement patterns"), times=5)))
colnames(features) <- "Feature"
physig <- cbind(phylogeny, features, physig)
write.csv(physig, file=here("output_tables", "Table_SI_phylosig_continuous.csv"), row.names = FALSE)

rownames(physig) <- NULL

physig_latex <- physig %>%
  kbl(caption="Phylogenetic signal",
    format="latex") %>% #,
  kable_minimal(full_width = F) %>%
  kable_styling(latex_options = c("scale_down"))  %>%
  column_spec(1, width = "8em")
