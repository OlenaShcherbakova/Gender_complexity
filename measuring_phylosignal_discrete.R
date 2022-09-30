#measuring phylosignal of discrete features

#Indo-European tree (1): Bouckaert et al 2012
source("library.R")
gb <- load_data_final()

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

data <- comparative.data(phy = tree, data = gb_geo.subset, names.col = taxon, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

#from function description: "The value of D depends on phylogeny size - more sister clades yield higher sums - and so the means of the two sets of simulated data are used as calibrations to scale both observed and simulated values of D to set points of 0 (as phylogenetically conserved as expected under a Brownian threshold model) and 1 (random). The value of D can be both smaller than 0 (highly conserved) and greater than 1 (overdispersed) and the distributions of scaled D from the simulations are used to assess the significance of the observed scaled D."
#DEstimate: The estimated D value
#Pval1: A p value, giving the result of testing whether D is significantly different from one
#Pval0: A p value, giving the result of testing whether D is significantly different from zero

physig_phon_prop_D <- caper::phylo.d(data=data, binvar = phon_prop, permut = 1000)
physig_phon_prop_DEstimate <- physig_phon_prop_D$DEstimate[1][["Obs"]]
physig_phon_prop_Pval1 <- physig_phon_prop_D$Pval1 # D compared to 1
physig_phon_prop_Pval0 <- physig_phon_prop_D$Pval0 # D compared to 0
phon_prop_ie <- c(physig_phon_prop_DEstimate, physig_phon_prop_Pval1, physig_phon_prop_Pval0, "Indo-European")

physig_unpredictable_D <- caper::phylo.d(data=data, binvar = unpredictable, permut = 1000)
physig_unpredictable_DEstimate <- physig_unpredictable_D$DEstimate[1][["Obs"]]
physig_unpredictable_Pval1 <- physig_unpredictable_D$Pval1 # D compared to 1
physig_unpredictable_Pval0 <- physig_unpredictable_D$Pval0 # D compared to 0
unpredictable_ie <- c(physig_unpredictable_DEstimate, physig_unpredictable_Pval1, physig_unpredictable_Pval0, "Indo-European")




#Bantu
source("library.R")
gb <- load_data_final()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/grollemund_et_al2015/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb %>%
  #dplyr::select(Glottocode, sem_classes, agr_patterns, taxon) %>%
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

data <- comparative.data(phy = tree, data = gb_geo.subset, names.col = taxon, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

#from function description: "The value of D depends on phylogeny size - more sister clades yield higher sums - and so the means of the two sets of simulated data are used as calibrations to scale both observed and simulated values of D to set points of 0 (as phylogenetically conserved as expected under a Brownian threshold model) and 1 (random). The value of D can be both smaller than 0 (highly conserved) and greater than 1 (overdispersed) and the distributions of scaled D from the simulations are used to assess the significance of the observed scaled D."
#DEstimate: The estimated D value
#Pval1: A p value, giving the result of testing whether D is significantly different from one
#Pval0: A p value, giving the result of testing whether D is significantly different from zero

physig_phon_prop_D <- caper::phylo.d(data=data, binvar = phon_prop, permut = 1000)
physig_phon_prop_DEstimate <- physig_phon_prop_D$DEstimate[1][["Obs"]]
physig_phon_prop_Pval1 <- physig_phon_prop_D$Pval1 # D compared to 1
physig_phon_prop_Pval0 <- physig_phon_prop_D$Pval0 # D compared to 0
phon_prop_b <- c(physig_phon_prop_DEstimate, physig_phon_prop_Pval1, physig_phon_prop_Pval0, "Bantu")

physig_unpredictable_D <- caper::phylo.d(data=data, binvar = unpredictable, permut = 1000)
physig_unpredictable_DEstimate <- physig_unpredictable_D$DEstimate[1][["Obs"]]
physig_unpredictable_Pval1 <- physig_unpredictable_D$Pval1 # D compared to 1
physig_unpredictable_Pval0 <- physig_unpredictable_D$Pval0 # D compared to 0
unpredictable_b <- c(physig_unpredictable_DEstimate, physig_unpredictable_Pval1, physig_unpredictable_Pval0, "Bantu")






#World
source("library.R")
gb <- load_data_final()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/world/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb %>%
  #dplyr::select(Glottocode, sem_classes, agr_patterns, taxon) %>%
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

data <- comparative.data(phy = tree, data = gb_geo.subset, names.col = taxon, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

#from function description: "The value of D depends on phylogeny size - more sister clades yield higher sums - and so the means of the two sets of simulated data are used as calibrations to scale both observed and simulated values of D to set points of 0 (as phylogenetically conserved as expected under a Brownian threshold model) and 1 (random). The value of D can be both smaller than 0 (highly conserved) and greater than 1 (overdispersed) and the distributions of scaled D from the simulations are used to assess the significance of the observed scaled D."
#DEstimate: The estimated D value
#Pval1: A p value, giving the result of testing whether D is significantly different from one
#Pval0: A p value, giving the result of testing whether D is significantly different from zero

physig_phon_prop_D <- caper::phylo.d(data=data, binvar = phon_prop, permut = 1000)
physig_phon_prop_DEstimate <- physig_phon_prop_D$DEstimate[1][["Obs"]]
physig_phon_prop_Pval1 <- physig_phon_prop_D$Pval1 # D compared to 1
physig_phon_prop_Pval0 <- physig_phon_prop_D$Pval0 # D compared to 0
phon_prop_w <- c(physig_phon_prop_DEstimate, physig_phon_prop_Pval1, physig_phon_prop_Pval0, "World")

physig_unpredictable_D <- caper::phylo.d(data=data, binvar = unpredictable, permut = 1000)
physig_unpredictable_DEstimate <- physig_unpredictable_D$DEstimate[1][["Obs"]]
physig_unpredictable_Pval1 <- physig_unpredictable_D$Pval1 # D compared to 1
physig_unpredictable_Pval0 <- physig_unpredictable_D$Pval0 # D compared to 0
unpredictable_w <- c(physig_unpredictable_DEstimate, physig_unpredictable_Pval1, physig_unpredictable_Pval0, "World")



#Making a table
physig <- as.data.frame(rbind(phon_prop_ie, unpredictable_ie, phon_prop_b, unpredictable_b, phon_prop_w, unpredictable_w))

physig <- physig %>%
  as.data.frame() %>%
  mutate(across(c("V1", "V2", "V3"), as.numeric)) %>%
  mutate_if(is.numeric, ~round(., digits=2)) %>%
  rename(`D value` = V1) %>%
  rename(`p-value: departure from phylogenetic randomness` = V2) %>%
  rename(`p-value: departure from Brownian threshold model` = V3) %>%
  rename(`Phylogeny` = V4) %>%
  relocate(Phylogeny, .before = 1) %>%
  add_column(Feature=rep(c("Phonological rules", "Unpredictable"), times=3), .after = "Phylogeny")

rownames(physig) <- NULL

write.csv(physig, "output_tables/Phylogenetic_signal_D_value.csv", row.names = FALSE)

physig_latex <- physig %>%
  kbl(caption="Phylogenetic signal (discrete)",
      format="latex") %>% #,
  kable_minimal(full_width = F) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  column_spec(1, width = "8em")

