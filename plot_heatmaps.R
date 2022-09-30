#plotting three heatmaps as one plot

source("plot_heatmap_ie.R")
source("plot_heatmap_b.R")
source("plot_heatmap_world_labelled.R")

all_heatmaps <- heatmap_world | (heatmap_b / heatmap_ie)

all_heatmaps <- all_heatmaps + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 30))

ggsave(file="output/heatmaps.svg", plot=all_heatmaps, width=22, height=16)
