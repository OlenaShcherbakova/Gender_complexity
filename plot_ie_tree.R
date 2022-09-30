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

gb <- load_data_final()

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

dt_Sem_classes <- as.matrix(gb.subset)[,2]
dt_Agr_rules <- as.matrix(gb.subset)[,3]
fit_Sem_classes <- phytools::fastAnc(tree,dt_Sem_classes,vars=TRUE,CI=TRUE)
fit_Agr_rules <- phytools::fastAnc(tree,dt_Agr_rules,vars=TRUE,CI=TRUE)

td_Sem_classes <- data.frame(node = nodeid(tree, names(dt_Sem_classes)),
                 trait = dt_Sem_classes)
td_Agr_rules <- data.frame(node = nodeid(tree, names(dt_Agr_rules)),
                             trait = dt_Agr_rules)

nd_Sem_classes <- data.frame(node = names(fit_Sem_classes$ace), trait = fit_Sem_classes$ace)
nd_Agr_rules <- data.frame(node = names(fit_Agr_rules$ace), trait = fit_Agr_rules$ace)

d_Sem_classes <- rbind(td_Sem_classes, nd_Sem_classes)
d_Agr_rules <- rbind(td_Agr_rules, nd_Agr_rules)

d_Sem_classes$node <- as.numeric(d_Sem_classes$node)
d_Agr_rules$node <- as.numeric(d_Agr_rules$node)

d_Sem_classes$trait <- as.numeric(d_Sem_classes$trait)
d_Agr_rules$trait <- as.numeric(d_Agr_rules$trait)

tree_Sem_classes <- full_join(tree, d_Sem_classes, by = 'node')
tree_Agr_rules <- full_join(tree, d_Agr_rules, by = 'node')

#gb.subset <- gb.subset %>% 
#  rename(`semantic\nrules` = sem_classes,
#         `agreement\npatterns` = agr_patterns,
#         `phonological\nproperties` = phon_prop)

plot <- ggtree(tree_Sem_classes, aes(color=trait)) + geom_tiplab() +
  scale_color_viridis_c(option="D", name="semantic\nrules", direction=-1, breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1)) +
  theme(legend.position="bottom", legend.direction="horizontal") + 
  theme(text = element_text(size=18), legend.key.size = unit(1, 'cm')) #+ xlim(-1, 1)

plot_mirror <- ggtree(tree_Agr_rules, aes(color=trait)) + geom_tiplab(hjust=1) +
  scale_color_viridis_c(option="D", name="agreement\npatterns", direction=-1, breaks=c(0, 0.20, 0.40, 0.60, 0.80, 1), limits=c(0, 1)) + scale_x_reverse() +
  theme(legend.position="bottom", legend.direction="horizontal") + 
  theme(text = element_text(size=18), legend.key.size = unit(1, 'cm'))

joined_trees <- plot | plot_mirror

ggsave(file="output/plot_IE_tree.svg", plot=plot, width=15, height=12, scale=-0.5)

plain_tree <- ggtree(tree, layout = 'rect', branch.length='none') %<+% gb.subset + geom_tiplab() +
  scale_color_continuous(name='semantic rules', value=sem_classes, low='darkgreen', high='red') +
  theme(legend.position="right")



plain_tree <- ggtree(tree, layout = 'rect', branch.length='none') %<+% gb.subset + geom_tiplab()

p1 <- gheatmap(plain_tree, gb.subset[, 2:3], offset=7, width=.9, colnames_angle=0, 
               colnames_offset_y = .5, font.size=3, hjust=0.5, colnames_position = "top") + #ylim(-3, 1500) +
  scale_fill_viridis_c(option="D", name="continuous values:", direction=-1)

p2 <- p1 + new_scale_fill()

p2 <- gheatmap(p2, gb.subset[, 4:5], offset=17.7, width=.9, colnames_angle=0, 
               colnames_offset_y = .5, font.size=3, hjust=0.5, colnames_position = "top") + #ylim(-3, 1500) +
  scale_fill_viridis_d(option="D", name="discrete values:", direction=-1)

png("output/heatmap_ie_rect.png", res=300, width=2600, height=2600)
plot(p2)
dev.off()

