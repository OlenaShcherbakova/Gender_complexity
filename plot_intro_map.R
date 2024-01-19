#plotting two maps visualizing semantic classes and agreement patterns scores

source('library.R')

#subsetting the dataset to the languages present on global tree and used in the analyses

PHYLOGENIES <- list.files(here('data/phylogenies'), full.names = TRUE)
TAXAS <-
  list.files(
    here('data/phylogenies'),
    pattern = "taxa.csv",
    full.names = TRUE,
    recursive = TRUE
  )

trees <- load_trees(PHYLOGENIES[5])

# sample one tree from full posterior
tree <- sample(trees, 1)[[1]]

taxa <- read.csv(TAXAS[5])
gb <- load_data_intro_figure()


gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
#language_glottolog <- read.csv("data/languages_glottolog.csv")
language_glottolog <-
  read.csv(
    "https://raw.githubusercontent.com/glottolog/glottolog-cldf/v4.4/cldf/languages.csv"
  )


# merge datasets
gb.subset <- gb %>%
  left_join(language_glottolog %>% dplyr::select(Glottocode, Latitude, Longitude)) %>%
  dplyr::filter(!is.na(taxon))

gb.subset <- gb.subset[gb.subset$taxon %in% trees[[1]]$tip.label,]

# prune data and trees to match tips
tree <- drop.tip(tree,
                 setdiff(tree$tip.label, gb.subset$taxon))

# collapse polytomies if needed
if (!is.binary(tree)) {
  tree <- multi2di(tree)
}

# make sure we have a complete match between the tree and the data
stopifnot(all(tree$tip.label %in% gb.subset$taxon))

rownames(gb.subset) <- gb.subset$taxon

#four conditions
world <-
  map_data(
    'world',
    wrap = c(-25, 335),
    ylim = c(-56, 80),
    margin = T
  )

lakes <-
  map_data(
    "lakes",
    wrap = c(-25, 335),
    col = "white",
    border = "gray",
    ylim = c(-55, 65),
    margin = T
  )

#shifting the longlat of the dataframe to match the pacific centered map
combination <- gb.subset %>%
  mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude)) %>%
  mutate(Gender =
           recode(Gender,
                  "0" = "absent", #first OLD then NEW
                  "1" = "present"))


#Basemap
basemap <- ggplot(combination) +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    colour = "gray92",
    fill = "gray92",
    size = 0.5
  ) +
  geom_polygon(
    data = lakes,
    aes(x = long, y = lat, group = group),
    colour = "gray92",
    fill = "white",
    size = 0.3
  )  +
  theme(
    panel.grid.major = element_blank(),
    #all of these lines are just removing default things like grid lines, axises etc
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )   +
  coord_map(projection = "vandergrinten", ylim = c(-56, 67))

custom_cols = viridis::viridis(40)
custom_cols <- custom_cols[c(40, 1)]

m1 <- basemap + geom_point(
  aes(x = Longitude, y = Latitude, color = Gender),
  alpha = 0.5,
  stat = "identity",
  size = 1.8
)
m1 <- m1 + scale_color_manual(values = custom_cols) +
  labs(color = "Gender/noun class") + theme(text = element_text(size = 18),
                                            legend.key.size = unit(1, 'cm'))
m1

ggsave(
  file = "output/plot_map_intro.svg",
  plot = m1,
  width = 10,
  height = 5
)
ggsave(
  file = "output/plot_map_intro.pdf",
  plot = m1,
  width = 10,
  height = 5
)
ggsave(
  file = "output/plot_map_intro.jpg",
  plot = m1,
  width = 10,
  height = 5,
  dpi = 300
)
