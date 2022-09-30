#custom plots from phylogenetic path analysis output


#Plot 1: model sets tested on three phylogenies

#creating a custom plot of the tested causal models in the phylogenetic path analysis on 1705 languages for three following variables: NC, Vf, FWO

library(here)
source(here('library.R'))

TREESET <- "data/phylogenies/world/"

trees <- load_trees(TREESET)
grambank_phylopath_compl <- load_data_final() 

grambank_phylopath_compl <- grambank_phylopath_compl %>%
  filter(sem_classes != "0" | agr_patterns != "0")

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv("data/phylogenies/world/taxa.csv")
grambank_phylopath_compl$taxon <- taxa$taxon[match(grambank_phylopath_compl$Glottocode, taxa$glottocode)]
grambank_phylopath_compl <- grambank_phylopath_compl %>%
  dplyr::filter(!is.na(taxon)) 

grambank_phylopath_compl <- grambank_phylopath_compl[grambank_phylopath_compl$taxon %in% trees[[1]]$tip.label, ]

# prune data and trees to match tips
tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, grambank_phylopath_compl$taxon)
)

# collapse polytomies if needed
if (!is.binary(tree)) {
  tree <- multi2di(tree)
}

# make sure we have a complete match between the tree and the data
stopifnot(all(tree$tip.label %in% grambank_phylopath_compl$taxon))

grambank_phylopath_compl %>%
  dplyr::select(Glottocode, sem_classes, agr_patterns, phon_prop, unpredictable, taxon) %>%
  mutate_at(c("phon_prop", "unpredictable"), as.factor) -> grambank_phylopath_compl

grambank_phylopath_compl$sem_classes <- as.numeric(as.character(grambank_phylopath_compl$sem_classes))
grambank_phylopath_compl$agr_patterns <- as.numeric(as.character(grambank_phylopath_compl$agr_patterns))


rownames(grambank_phylopath_compl) <- grambank_phylopath_compl$taxon

m <- define_model_set(
  null = c(),
  #simple models with one predictor
  a1 = c(agr_patterns~sem_classes),
  a2 = c(agr_patterns~phon_prop),
  a3 = c(agr_patterns~unpredictable),
  
  #models with several distinct predictors: 2 and 3
  b1 = c(agr_patterns~sem_classes + phon_prop),
  b2 = c(agr_patterns~sem_classes + phon_prop + unpredictable),
  
  #models where agreement patterns depend on only one of the predictors, but two predictors are casually connected among themselves
  c1 = c(agr_patterns~sem_classes, sem_classes~phon_prop),
  c2 = c(agr_patterns~phon_prop, phon_prop~sem_classes),
  
  #extensions of 'c' models: 'unpredictable' is casually linked with one of the other precictors 
  d1 = c(agr_patterns~sem_classes, sem_classes~phon_prop + unpredictable),
  d2 = c(agr_patterns~phon_prop, phon_prop~sem_classes + unpredictable)
)

positions <- data.frame(
  name = c('sem_classes', 'agr_patterns', 'phon_prop', 'unpredictable'),
  x = c(0, 0, 1, 1),
  y = c(1, 0, 0, 1)
)

source("customizing_phylopath_plot_model_set_function.R")

model_sets <- plot_model_set(m, manual_layout = positions, edge_width = 1)

names <- c('sem_classes', 'agr_patterns', 'phon_prop', 'unpredictable')
custom_names <- c("semantic\nrules", "agreement\npatterns", "phonological\nrules", "unpredictable")
named_vector <- setNames(custom_names, names)

#plot_model_set(m, manual_layout = positions)
model_sets <- plot_model_set(m, manual_layout = positions, edge_width = 3, curvature=0,
                             text_size = 21, labels=named_vector,
                             box_x=130, box_y=90,
                             arrow = grid::arrow(type = 'closed', 30, grid::unit(12, 'points')))
#model_sets

ggsave(file="output/phylopath_model_sets_custom.svg", plot=model_sets, width=45, height=30)

#Plot 2: visualized best causal models with confidence intervals for each path
source("phylopath_plots_ie.R")
source("phylopath_plots_b.R")
source("phylopath_plots_w.R")

all <- (p_ie | p_b | p_w) / (coef_ie | coef_b | coef_w)
ggsave(file="output/phylopath_all.svg", plot=all, width=17.5, height=10.5)

