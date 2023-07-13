#double checking the results on several posterior trees

library(here)
source(here('library.R'))

TREESET <- "data/phylogenies/bouckaert_et_al2012/"

trees <- load_trees(TREESET)
grambank_phylopath_compl <- load_data_final()


coefficients_matrix <- matrix(NA, 100, 16)

listcombo <- list(c(1, 1), c(1, 2), c(1, 3), c(1, 4), 
                  c(2, 1), c(2, 2), c(2, 3), c(2, 4),
                  c(3, 1), c(3, 2), c(3, 3), c(3, 4),
                  c(4, 1), c(4, 2), c(4, 3), c(4, 4))
predterms <- lapply(listcombo, function(x) paste(c("1","2","3","4")[x],collapse="_on_"))
predterms <- t(as.data.frame(predterms))
predterms <- gsub("1", "unpredictable", predterms, fixed=TRUE)
predterms <- gsub("2", "phon_prop", predterms, fixed=TRUE)
predterms <- gsub("3", "sem_classes", predterms, fixed=TRUE)
predterms <- gsub("4", "agr_patterns", predterms, fixed=TRUE)

colnames(coefficients_matrix) <- predterms

for (i in 1:length(trees[1:100])) {
  tryCatch({
    # sample one tree from full posterior
    tree <- trees[[i]]
    #tree <- trees[[5]]
    
    
    taxa <- read.csv("data/phylogenies/bouckaert_et_al2012/taxa.csv")
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
      filter(!is.na(Glottocode, sem_classes, agr_patterns, phon_prop, unpredictable, taxon)) %>%
      mutate_at(c("sem_classes", "agr_patterns"), as.numeric) %>%
      mutate_at(c("phon_prop", "unpredictable"), as.factor) -> grambank_phylopath_compl
    
    rownames(grambank_phylopath_compl) <- grambank_phylopath_compl$taxon
    
    m <- define_model_set(
      null = c(),
      a = c(agr_patterns~sem_classes),
      b = c(agr_patterns~unpredictable),
      c = c(agr_patterns~phon_prop),
      d = c(agr_patterns~sem_classes + unpredictable),
      e = c(agr_patterns~unpredictable + phon_prop),
      f = c(agr_patterns~sem_classes + phon_prop),
      g = c(agr_patterns~sem_classes + unpredictable + phon_prop),
      h = c(agr_patterns~sem_classes + unpredictable + phon_prop, phon_prop~unpredictable)
    )
    
    
    positions <- data.frame(
      name = c('sem_classes', 'agr_patterns', 'phon_prop', 'unpredictable'),
      x = c(0, 0, 1, 1),
      y = c(1, 0, 0, 1)
    )
    
    p <- phylo_path(m, grambank_phylopath_compl, tree)
    
    #"Specifically, it reports the model name, the number of independence claims made by the model (k), the number of parameters (q), the C statistic and the accompanying p-value. A significant p-value would indicate that the available evidence rejects the model. It also reports model selection information: the C-statistic information criterion corrected for small sample sizes (CICc), the difference in CICc with the top model (delta_CICc) and finally the associated relative likelihoods (l) and CICc weights (w)" (van der Bijl 2018)
    s <- summary(p)
    s
    plot(s)
    
    avg <- average(p)
    
    coefficients_matrix[i, 1] <- avg$coef["unpredictable", "unpredictable"]
    coefficients_matrix[i, 2] <- avg$coef["unpredictable", "phon_prop"]
    coefficients_matrix[i, 3] <- avg$coef["unpredictable", "sem_classes"]
    coefficients_matrix[i, 4] <- avg$coef["unpredictable", "agr_patterns"]
    
    coefficients_matrix[i, 5] <- avg$coef["phon_prop", "unpredictable"]
    coefficients_matrix[i, 6] <- avg$coef["phon_prop", "phon_prop"]
    coefficients_matrix[i, 7] <- avg$coef["phon_prop", "sem_classes"]
    coefficients_matrix[i, 8] <- avg$coef["phon_prop", "agr_patterns"]
    
    coefficients_matrix[i, 9] <- avg$coef["sem_classes", "unpredictable"]
    coefficients_matrix[i, 10] <- avg$coef["sem_classes", "phon_prop"]
    coefficients_matrix[i, 11] <- avg$coef["sem_classes", "sem_classes"]
    coefficients_matrix[i, 12] <- avg$coef["sem_classes", "agr_patterns"]
    
    coefficients_matrix[i, 13] <- avg$coef["agr_patterns", "unpredictable"]
    coefficients_matrix[i, 14] <- avg$coef["agr_patterns", "phon_prop"]
    coefficients_matrix[i, 15] <- avg$coef["agr_patterns", "sem_classes"]
    coefficients_matrix[i, 16] <- avg$coef["agr_patterns", "agr_patterns"]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}




coefficients_matrix_plot <- as.data.frame(coefficients_matrix) %>%
  select_if(colSums(.) != 0) %>%
  gather(key="paths", value="value", everything()) 

coefficients_matrix_plot$paths <- gsub("_on_", "\U2192", coefficients_matrix_plot$paths)

coefficients_matrix_plot_final <- coefficients_matrix_plot %>%
  ggplot(aes(x=paths, y=value)) + 
  geom_boxplot() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 5, vjust =1, hjust=0.6, size = 8))
coefficients_matrix_plot_final

png("output/double_checking_coefficients_ie.png", width=1500, height=1500, res = 300)
coefficients_matrix_plot_final
dev.off()