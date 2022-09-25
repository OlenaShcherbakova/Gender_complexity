library(here)
source(here('library.R'))

PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)

#Indo-European (Bouckaert et al. 2012)
trees <- load_trees(PHYLOGENIES[1])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[1])

gb <- load_data()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:19]


rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
               colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 30) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","ie_tree.png"), p, dpi=300, dev="png")


#Pama-Nyungan
trees <- load_trees(PHYLOGENIES[2])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[2])

gb <- load_data()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:19]


rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 80) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","pn_tree.png"), p, dpi=300, dev="png")


#Pama-Nyungan: agreement
trees <- load_trees(PHYLOGENIES[2])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[2])

gb <- load_data_NAs_agreement()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:5]


rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 80) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","pn_tree.png"), p, dpi=300, dev="png")



#Austronesian
trees <- load_trees(PHYLOGENIES[6])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[6])

gb <- load_data()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:19]


rect <- ggtree(tree, branch.length='none', layout='rect') #+ geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 250) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","a_tree.png"), p, dpi=300, dev="png")



#Austronesian:agreement
trees <- load_trees(PHYLOGENIES[6])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[6])

gb <- load_data()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:5]


rect <- ggtree(tree, branch.length='none', layout='rect') #+ geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 250) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","a_tree.png"), p, dpi=300, dev="png")





#Sino-Tibetan
trees <- load_trees(PHYLOGENIES[14])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[14])

gb <- load_data()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:19]


rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 38) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","st_tree.png"), p, dpi=300, dev="png")



#Kitchen: agreement

trees <- load_trees(PHYLOGENIES[10])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[10])

gb <- load_data_NAs_agreement()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:5]


rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 38) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","st_tree.png"), p, dpi=300, dev="png")




#Arawak (walker_and_ribeiro2011)
trees <- load_trees(PHYLOGENIES[16])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[16])

gb <- load_data()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:19]


rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 35) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","arawak_tree.png"), p, dpi=300, dev="png")



#Dravidian (Kolipakam et al. 2018)
trees <- load_trees(PHYLOGENIES[11])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[11])

gb <- load_data()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:19]


rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 20) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","dravidian_tree.png"), p, dpi=300, dev="png")




#World (Jager 2018)
trees <- load_trees(PHYLOGENIES[17])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[17])

gb <- load_data_NAs_semantics()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
gb.subset <- gb.subset %>% select(taxon, everything())
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset %>%
  group_by(Family_ID) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) -> table


gb.subset <- gb.subset %>%
  mutate(family = case_when(
    Family_ID == "aust1307"  ~ "Austronesian", #Austronesian
    Family_ID == "aust1305"  ~ "Austroasiatic", #Austroasiatic
    Family_ID == "indo1319"  ~ "Indo-European", #Indo-European
    Family_ID == "ural1272" ~ "Uralic", #Uralic
    Family_ID == "pama1250"  ~ "Pama-Nyungan", #Pama-Nyungan
    Family_ID == "nucl1709" ~ "Nuclear Trans New Guinea", #Nuclear Trans New Guinea
    Family_ID == "atla1278"  ~ "Atlantic-Congo", #Atlantic-Congo
    Family_ID == "afro1255"  ~ "Afro-Asiatic", #Afro-Asiatic
    Family_ID == "drav1251" ~ "Dravidian", #Arawakan
    Family_ID == "utoa1244" ~ "Uto-Aztecan", #Uto-Aztecan
    Family_ID == "sino1245" ~ "Sino-Tibetan", #Sino-Tibetan
    Family_ID == "cent2225" ~ "Central Sudanic",
    Family_ID != "aust1307" & Family_ID != "aust1305" & Family_ID != "indo1319" & Family_ID != "ural1272" &
      Family_ID != "pama1250" & Family_ID != "nucl1709" & Family_ID != "atla1278" & Family_ID != "afro1255" &
      Family_ID != "drav1251" & Family_ID != "utoa1244" & Family_ID != "sino1245" & Family_ID != "cent2225" ~ "other"
  ))

gb.subset$family <- ordered(gb.subset$family, levels = c("Austronesian", "Austroasiatic",  "Indo-European",  "Uralic",  "Pama-Nyungan",  "Nuclear Trans New Guinea",  "Atlantic-Congo",  "Afro-Asiatic", "Dravidian",   "Sino-Tibetan", "Uto-Aztecan",  "Central Sudanic", "other"))
cols=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#000000")


coloured_branches <- ggtree(tree, layout = 'circular', branch.length='none') %<+% gb.subset + geom_tippoint(aes(color=family))
coloured_branches <- coloured_branches + scale_colour_manual(values = cols)


p <- gheatmap(coloured_branches, gb.subset[, 3:8], offset=2, width=.6, colnames_angle=75, 
              colnames_offset_y = .02, font.size=3, hjust=1) + ylim(-3, 1600) +
  scale_fill_viridis_d(option="D", name="values:")
p
ggsave(filename = here("output","plot_global_tree_sem_circ_fams.png"), p, dpi=300, dev="png")


coloured_branches <- ggtree(tree, layout = 'rect', branch.length='none') %<+% gb.subset + geom_tippoint(aes(color=family))
coloured_branches <- coloured_branches + scale_colour_manual(values = cols)

p <- gheatmap(coloured_branches, gb.subset[, 3:8], offset=2, width=.6, colnames_angle=15, 
              colnames_offset_y = .02, font.size=3, hjust=1) + ylim(-3, 1600) +
  scale_fill_viridis_d(option="D", name="values:")
p
ggsave(filename = here("output","plot_global_tree_sem_rect_fams.png"), p, dpi=300, dev="png")






rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 500) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","dravidian_tree.png"), p, dpi=300, dev="png")





#Jager (world tree) agreement
trees <- load_trees(PHYLOGENIES[17])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[17])

gb <- load_data_NAs_agreement()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
gb.subset <- gb.subset %>% select(taxon, everything())

gb.subset <- gb.subset %>%
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

gb.subset %>%
  group_by(Family_ID) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) -> table

gb.subset <- gb.subset %>%
  mutate(family = case_when(
    Family_ID == "aust1307"  ~ "Austronesian", #Austronesian
    Family_ID == "aust1305"  ~ "Austroasiatic", #Austroasiatic
    Family_ID == "indo1319"  ~ "Indo-European", #Indo-European
    Family_ID == "ural1272" ~ "Uralic", #Uralic
    Family_ID == "pama1250"  ~ "Pama-Nyungan", #Pama-Nyungan
    Family_ID == "nucl1709" ~ "Nuclear Trans New Guinea", #Nuclear Trans New Guinea
    Family_ID == "atla1278"  ~ "Atlantic-Congo", #Atlantic-Congo
    Family_ID == "afro1255"  ~ "Afro-Asiatic", #Afro-Asiatic
    Family_ID == "chib1249" ~ "Chibchan", #Arawakan
    Family_ID == "utoa1244" ~ "Uto-Aztecan", #Uto-Aztecan
    Family_ID == "sino1245" ~ "Sino-Tibetan", #Sino-Tibetan
    Family_ID == "cent2225" ~ "Central Sudanic",
    Family_ID != "aust1307" & Family_ID != "aust1305" & Family_ID != "indo1319" & Family_ID != "ural1272" &
      Family_ID != "pama1250" & Family_ID != "nucl1709" & Family_ID != "atla1278" & Family_ID != "afro1255" &
      Family_ID != "chib1249" & Family_ID != "utoa1244" & Family_ID != "sino1245" & Family_ID != "cent2225" ~ "other"
  ))

gb.subset$family <- ordered(gb.subset$family, levels = c("Austronesian", "Austroasiatic",  "Indo-European",  "Uralic",  "Pama-Nyungan",  "Nuclear Trans New Guinea",  "Atlantic-Congo",  "Afro-Asiatic", "Uto-Aztecan", "Chibchan", "Sino-Tibetan",   "Central Sudanic", "other"))
cols=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#000000")


coloured_branches <- ggtree(tree, layout = 'circular', branch.length='none') %<+% gb.subset + geom_tippoint(aes(color=family))
coloured_branches <- coloured_branches + scale_colour_manual(values = cols)


p <- gheatmap(coloured_branches, gb.subset[, 3:7], offset=2, width=.6, colnames_angle=75, 
              colnames_offset_y = .02, font.size=3, hjust=1) + ylim(-3, 1600) +
  scale_fill_viridis_d(option="D", name="values:")
p
ggsave(filename = here("output","plot_global_tree_agr_simple_circ_fams.png"), p, dpi=300, dev="png")


coloured_branches <- ggtree(tree, layout = 'rect', branch.length='none') %<+% gb.subset + geom_tippoint(aes(color=family))
coloured_branches <- coloured_branches + scale_colour_manual(values = cols)

p <- gheatmap(coloured_branches, gb.subset[, 3:7], offset=2, width=.6, colnames_angle=15, 
              colnames_offset_y = .02, font.size=3, hjust=1) + ylim(-3, 1600) +
  scale_fill_viridis_d(option="D", name="values:")
p
ggsave(filename = here("output","plot_global_tree_agr_simple_rect_fams.png"), p, dpi=300, dev="png")










rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)
rect <- ggtree(tree, branch.length='none', layout='circ') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 500) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","world_tree_agreement.png"), p, dpi=300, dev="png")


#Jager (world tree) semantics
trees <- load_trees(PHYLOGENIES[17])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[17])

gb <- load_data_NAs_semantics()

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
#gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "sex", "shape", "animacy", "plant_status", "phon_prop", "unpredictable", "adj", "dem", "art", "numerals", "pron_1_m_f", "pron_2_m_f", "pron_3_gender", "dem_classifiers", "numeral_classifiers", "poss_classifiers", "augm_meaning", "dimin_meaning"))

gb.subset <- gb.subset %>%
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

gb.subset <- gb.subset[, 2:5]


rect <- ggtree(tree, branch.length='none', layout='rect') + geom_tiplab(size=3)
rect <- ggtree(tree, branch.length='none', layout='circ') + geom_tiplab(size=3)


p <- gheatmap(rect, gb.subset, offset=2, width=.6,
              colnames_angle=85, colnames_offset_y = .02, font.size=3) + ylim(-3, 500) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

ggsave(filename = here("output","world_tree_semantics.png"), p, dpi=300, dev="png")





p1 <- gheatmap(circ, gb.subset, offset=.8, width=.2,
               colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue") +
  geom_tippoint(aes(color=Nominal_words_complexity))

p1 <- gheatmap(rect, gb.subset, offset=.15, width=.3,
               colnames_angle=90, colnames_offset_y = .02) + ylim(-3, 30) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

p1 <- gheatmap(circ, gb.subset, width=) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

circ <- ggtree(tree, layout = "circular")
circ <- ggtree(tree, branch.length='none', layout='circular')
circ <- ggtree(tree, branch.length='none')
gheatmap(p, gb.subset$Free_word_order, offset=5, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("-", "0", "1"), 
                    values=c("steelblue", "firebrick", "darkgreen"), name="Free_word_order")

