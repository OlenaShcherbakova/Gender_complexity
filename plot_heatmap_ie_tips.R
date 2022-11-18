#Indo-European

library(here)
source(here('library.R'))

PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)

#Indo-European (Bouckaert et al. 2012): 1
trees <- load_trees(PHYLOGENIES[1])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[1])

gb <- load_data_final_short()

gb <- gb %>%
  filter(sem_classes != "0" | agr_patterns != "0")

gb$sem_classes <- as.numeric(as.character(gb$sem_classes))
gb$agr_patterns <- as.numeric(as.character(gb$agr_patterns))

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]

gb.subset <- gb %>%
  dplyr::filter(!is.na(taxon))

gb.subset <- gb.subset[gb.subset$taxon %in% trees[[1]]$tip.label, ]

# prune data and trees to match tips
tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, gb.subset$taxon)
)

# collapse polytomies if needed
if (!is.binary(tree)) {
  tree <- multi2di(tree)
}

# make sure we have a complete match between the tree and the data
stopifnot(all(tree$tip.label %in% gb.subset$taxon))

rownames(gb.subset) <- gb.subset$taxon

gb.subset$sem_classes <- as.numeric(as.character(gb.subset$sem_classes))
gb.subset$agr_patterns <- as.numeric(as.character(gb.subset$agr_patterns))

gb.subset <- gb.subset %>% 
  rename(`semantic rules` = sem_classes,
         `agreement patterns` = agr_patterns,
         `phonological rules` = phon_prop)

plain_tree <- ggtree(tree, layout = 'rect', branch.length='none') %<+% gb.subset + geom_tiplab(size=5)

p1 <- gheatmap(plain_tree, gb.subset[, 2:3], offset=5, width=1, colnames_angle=0, 
               colnames_offset_y = .5, font.size=5, hjust=0.5, color = NULL, colnames_position = "top") + #ylim(-3, 1500) +
  scale_fill_viridis_c(option="D", name="continuous values:", direction=-1, limits=c(0, 1))

p2 <- p1 + new_scale_fill()

p2 <- gheatmap(p2, gb.subset[, 4:5], offset=18, width=1, colnames_angle=0, 
               colnames_offset_y = .5, font.size=5, hjust=0.5, color = NULL, colnames_position = "top") + #ylim(-3, 1500) +
  scale_fill_viridis_d(option="D", name="discrete values:", direction=-1) + 
  theme(legend.box = "horizontal", legend.position = "bottom",
        text = element_text(size = 25),
        legend.key.size = unit(1.3, 'cm'))

ggsave(file="output/heatmap_ie_rect.svg", plot=p2, width=14, height=14)

