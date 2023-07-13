#Bantu tree: interdependence between semantic classes, agreement patterns, phonological properties as a factor in class assignment, and presence of a class that is phonologically and semantically unpredictable

library(here)
source(here('library.R'))

TREESET <- "data/phylogenies/grollemund_et_al2015/"

trees <- load_trees(TREESET)
grambank_phylopath_compl <- load_data_final_short() 

grambank_phylopath_compl <- grambank_phylopath_compl %>%
  filter(sem_classes != "0" | agr_patterns != "0")

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv("data/phylogenies/grollemund_et_al2015/taxa.csv")
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


model_sets <- plot_model_set(m, manual_layout = positions, edge_width = 1)


p <- phylo_path(m, grambank_phylopath_compl, tree, upper.bound=1, lower.bound=0) #setting upper and lower boundaries as required by the warning

#"Specifically, it reports the model name, the number of independence claims made by the model (k), the number of parameters (q), the C statistic and the accompanying p-value. A significant p-value would indicate that the available evidence rejects the model. It also reports model selection information: the C-statistic information criterion corrected for small sample sizes (CICc), the difference in CICc with the top model (delta_CICc) and finally the associated relative likelihoods (l) and CICc weights (w)" (van der Bijl 2018)
s <- summary(p)
s
plot(s)

table <- cbind(s$model, s$CICc, s$delta_CICc, s$l, s$w, s$k, s$q, s$C, s$p) %>%
  as.data.frame() 
colnames(table) <- c("Model", "CICc", "delta_CICc", "Relative likelihood", "CICc weight", "n of independence claims", "n of parameters", "C statistic", "p-value")

Model <- c("null", "a1", "a2", "a3", "b1", "b2", "c1", "c2", "d1", "d2")
verbose_names <- c("null", 
                   "agreement patterns ~ semantic rules", 
                   "agreement patterns ~ phonological rules", 
                   "agreement patterns ~ unpredictable",
                   
                   "agreement patterns ~ semantic rules + phonological rules",
                   "agreement patterns ~ semantic rules + phonological rules + unpredictable",
                   
                   "agreement patterns ~ semantic rules, semantic rules ~ phonological rules",
                   "agreement patterns ~ phonological rules, phonological rules ~ semantic rules",
                   
                   "agreement patterns ~ semantic rules, semantic rules ~ phonological rules + unpredictable",
                   "agreement patterns ~ phonological rules, phonological rules  ~ semantic rules + unpredictable")

model_names_df <- as.data.frame(cbind(Model, verbose_names))

table_final <- table %>% 
  left_join(model_names_df) %>%
  dplyr::select(-Model) %>%
  relocate(verbose_names, everything()) %>%
  rename(Model = verbose_names) %>%
  mutate_at(c(2:9), as.numeric) %>%
  mutate(across(c(2:8), round, 2)) %>%
  mutate(`p-value` = case_when(
    `p-value` < 0.001 ~ "<0.001",
    `p-value` < 0.01 ~ "<0.01",
    TRUE ~ as.character(round(`p-value`,2))
  )) %>% 
  mutate(Phylogeny = "Bantu") %>%
  relocate(Phylogeny, .after = Model)

write.csv(table_final, file="output_tables/Phylopath_output_table_B.csv")

table_final <- table_final %>%
  flextable() %>%
  autofit() %>%
  merge_v(j=c("Phylogeny"))

save_as_docx(
  "Phylogenetic path analysis (Bantu)" = table_final, 
  path = "output_tables/Phylopath_output_table_B.docx")

table_final_latex <- table %>%
  left_join(model_names_df) %>%
  dplyr::select(-Model) %>%
  relocate(verbose_names, everything()) %>%
  rename(Model = verbose_names) %>%
  mutate_at(c(2:9), as.numeric) %>%
  mutate(across(c(2:8), round, 2)) %>%
  mutate(`p-value` = case_when(
    `p-value` < 0.001 ~ "<0.001",
    `p-value` < 0.01 ~ "<0.01",
    TRUE ~ as.character(round(`p-value`,2))
  )) %>% 
  mutate(Phylogeny = "Bantu") %>%
  relocate(Phylogeny, .after = Model) %>%
  kbl(caption="Phylogenetic path analysis (Bantu)",
      format="latex") %>% #,
  kable_minimal(full_width = F) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  column_spec(1, width = "8em")


#png("output/phylopath_ie_default_summary.png", width=1300, height=850, res = 300)
#plot(s)
#dev.off()

b_ci <- best(p, boot = 500)
b_ci

#png("output/phylopath_ie_default_coef_plot.png", width=1300, height=850, res = 300)
#coef_plot(b_ci, error_bar = "se", order_by = "strength") + ggplot2::coord_flip()
#dev.off()

save(b_ci, 
     file = "output/phylopath_b.RData")
load('output/phylopath_b.RData')


b <- best(p)
plot(b, manual_layout = positions, text_size = 4.5)

coef_plot(b_ci, error_bar = "ci", order_by = "strength") + ggplot2::coord_flip()
