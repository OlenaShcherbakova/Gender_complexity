#world

library(here)
source(here('library.R'))

PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)

#world
trees <- load_trees(PHYLOGENIES[5])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[5])

gb <- load_data_final_short() 

language_glottolog <- read.csv("https://raw.githubusercontent.com/glottolog/glottolog-cldf/v4.4/cldf/languages.csv")

# merge datasets
gb <- gb %>%
  left_join(language_glottolog %>% dplyr::select(Glottocode, Latitude, Longitude, Family_ID))

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


### Adding colored branches of the biggest families in the dataset

gb.subset %>%
  group_by(Family_ID) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(n)) %>%
  filter(!Family_ID == "") %>%
  top_n(9, freq) -> table

biggest_families <- table$Family_ID
gb.subset$family_status <- NA
gb.subset$family_status <- ifelse(gb.subset$Family_ID %in% biggest_families, gb.subset$Family_ID, "other")

#double-checking if all families indeed converted to names and none is left with "NA"
unique(gb.subset$family_status)

gb.subset <- gb.subset %>% 
  mutate(family = 
           recode(family_status,
                  "aust1307"  = "Austronesian",
                  "aust1305"  = "Austroasiatic",
                  "indo1319"  = "Indo-European",
                  "atla1278" = "Atlantic-Congo",
                  "gong1255" = "Ta-Ne-Omotic",
                  "sino1245" = "Sino-Tibetan",
                  "afro1255"  = "Afro-Asiatic",
                  "otom1299" = "Otomanguean",
                  "drav1251" = "Dravidian",
                  "other" = "other")) 

#double-checking if all families indeed converted to names and none is left with "NA"
unique(gb.subset$family)

tips_lists <- vector(mode = "list", length=9)

for (f in 1:length(biggest_families)) {
  tips_lists[[f]] <- gb.subset[gb.subset$Family_ID == biggest_families[f],]$taxon
}

#the correct order within biggest families is preserved and the Glottocodes are replaced with suitable family name labels
biggest_families_verbose <- recode(biggest_families,
                                   "aust1307"  = "Austronesian",
                                   "aust1305"  = "Austroasiatic",
                                   "indo1319"  = "Indo-European",
                                   "atla1278" = "Atlantic-Congo",
                                   "gong1255" = "Ta-Ne-Omotic",
                                   "sino1245" = "Sino-Tibetan",
                                   "afro1255"  = "Afro-Asiatic",
                                   "otom1299" = "Otomanguean",
                                   "drav1251" = "Dravidian",
                                   "other" = "other")

names(tips_lists) <- biggest_families_verbose

nodes <- vector(mode="character", length=length(biggest_families))

for (tips in 1:length(tips_lists)) {
  nodes[tips] <- getMRCA(tree, tips_lists[[tips]])
}

#test
#nodes <- vector(mode="character", length=1)
#nodes[1] <- getMRCA(tree, tips_lists[[2]])

nodes <- as.numeric(nodes)

coloured_branches <- groupClade(tree, nodes)
coloured_branches <- ggtree(coloured_branches, layout = 'rect', branch.length='none', size=0.5)

#gb.subset <- gb.subset %>% 
#  rename(`semantic\nclasses` = sem_classes,
#         `agreement\npatterns` = agr_patterns,
#         `phonological\nproperties` = phon_prop)

gb.subset <- gb.subset %>% 
  rename(`semantic rules` = sem_classes,
         `agreement patterns` = agr_patterns,
         `phonological rules` = phon_prop)

plain_tree <- ggtree(tree, layout = 'rect', branch.length='none') #%<+% gb.subset + geom_tiplab()

p1 <- gheatmap(plain_tree, gb.subset[, 2:3], offset=-5, width=.9, colnames_angle=20, 
               colnames_offset_y = 30, colnames_offset_x = 10, 
               # adding the following to avoid pale output
               color = NULL,
               font.size=8, hjust=0.5, colnames_position = "top") + 
  #ylim(-3, 85) +
  scale_fill_viridis_c(option="D", name="continuous\nvalues", direction=-1) +
  #scale_x_continuous(expand = c(0, 0)) + 
  #scale_y_continuous(expand = c(0.04, 0.04)) +
  geom_cladelabel(node=nodes[1], label=biggest_families_verbose[1], offset=58, align=TRUE, fontsize=9) +
  geom_cladelabel(node=nodes[2], label=biggest_families_verbose[2], offset=58, align=TRUE, fontsize=9) +
  geom_cladelabel(node=nodes[3], label=biggest_families_verbose[3], offset=58, align=TRUE, fontsize=9) +
  geom_cladelabel(node=nodes[4], label=biggest_families_verbose[4], offset=58, align=TRUE, fontsize=9) +
  geom_cladelabel(node=nodes[5], label=biggest_families_verbose[5], offset=58, align=TRUE, fontsize=9) +
  geom_cladelabel(node=nodes[6], label=biggest_families_verbose[6], offset=58, align=TRUE, fontsize=9) +
  geom_cladelabel(node=nodes[7], label=biggest_families_verbose[7], offset=58, align=TRUE, fontsize=9) +
  geom_cladelabel(node=nodes[8], label=biggest_families_verbose[8], offset=58, align=TRUE, fontsize=9) +
  geom_cladelabel(node=nodes[9], label=biggest_families_verbose[9], offset=58, align=TRUE, fontsize=9)
p1
p2 <- p1 + new_scale_fill() 

p2 <- gheatmap(p2, gb.subset[, 4:5], offset=23, width=.9, colnames_angle=20, 
               colnames_offset_y = 30, colnames_offset_x = 10, 
               # adding the following to avoid pale output
               color = NULL,
               font.size=8, hjust=0.5, colnames_position = "top") + 
  #ylim(-3, 85) +
  scale_fill_viridis_d(option="D", name="discrete\nvalues", direction=-1) + 
  theme(legend.box = "horizontal", legend.position = "bottom",
        text = element_text(size = 25),
        legend.key.size = unit(1.3, 'cm'))+ 
  xlim(-1, 115) + ylim(-3, 520)
p2

ggsave(file="output/heatmap_world_rect.svg", plot=p2, width=11, height=10)

heatmap_world <- p2
