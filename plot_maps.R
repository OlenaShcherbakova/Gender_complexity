#plotting two maps visualizing semantic classes and agreement patterns scores

source('library.R')

#subsetting the dataset to the languages present on global tree and used in the analyses 

PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)

trees <- load_trees(PHYLOGENIES[5])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[5])
gb <- load_data_final_short()

gb$sem_classes <- as.numeric(as.character(gb$sem_classes))
gb$agr_patterns <- as.numeric(as.character(gb$agr_patterns))

gb <- gb %>%
  filter(sem_classes > 0 | agr_patterns > 0) %>%
  mutate(sem_classes_10 = sem_classes * 10,
         agr_patterns_10 = agr_patterns * 10)


gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("data/languages_glottolog.csv")

# merge datasets
gb.subset <- gb %>%
  left_join(language_glottolog %>% dplyr::select(Glottocode, Latitude, Longitude)) %>%
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

#four conditions
world <- map_data('world', wrap=c(-25,335), ylim=c(-56,80), margin=T)

lakes <- map_data("lakes", wrap=c(-25,335), col="white", border="gray", ylim=c(-55,65), margin=T)

#gb.subset$sem_classes <- as.factor(gb.subset$sem_classes)
#gb.subset$agr_patterns <- as.factor(gb.subset$agr_patterns)

gb.subset$sem_classes <- as.numeric(as.character(gb.subset$sem_classes))
gb.subset$agr_patterns <- as.numeric(as.character(gb.subset$agr_patterns))

#shifting the longlat of the dataframe to match the pacific centered map
combination <- gb.subset %>%
  mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude))


#Basemap
basemap <- ggplot(combination) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               colour="gray92", 
               fill="gray92", size = 0.5) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
               colour="gray92",
               fill="white", size = 0.3)  +
  theme(
    panel.grid.major = element_blank(), #all of these lines are just removing default things like grid lines, axises etc
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())   +
  coord_map(projection = "vandergrinten", ylim=c(-56,67))

#custom_cols = c(brewer.pal(7, "RdYlBu")); custom_cols <- custom_cols[c(7:1)]
#custom_cols = c(brewer.pal(7, "RdYlBu")); custom_cols <- custom_cols[c(7:1)]
custom_cols = viridis::viridis(40); custom_cols <- custom_cols[c(40:1)]
#custom_cols = c(brewer.pal(11, "Spectral")); custom_cols <- custom_cols[c(11:1)]


m1 <- basemap + geom_point(
  aes(x=Longitude, y=Latitude, color=sem_classes), alpha=0.7,
  stat = "identity", size=1.3
) 
m1 <- m1 + scale_colour_gradientn(colours = custom_cols) + labs(color = "Semantic assignment\nrules score") + theme(text = element_text(size=12)) #legend.title = element_text(size=12)

ggsave(file="output/map1_sem_classes.svg", plot=m1, width=12, height=6)


m2 <- basemap + geom_point(
  aes(x=Longitude, y=Latitude, color=agr_patterns), alpha=0.7,
  stat = "identity", size=1.3
) 
m2 <- m2 + scale_colour_gradientn(colours = custom_cols) + labs(color = "Agreement patterns\nscore") + theme(text = element_text(size=12))

ggsave(file="output/map2_agr_patterns.svg", plot=m2, width=12, height=6)

joined <- m1/m2

ggsave(file="output/combined_map1_map2.svg", plot=joined, width=12, height=12)
