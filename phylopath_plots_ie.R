library(here)
source(here('library.R'))

#customizing the plotting function of phylopath to produce plots of phylopath output

box_x = 25; box_y = 19; edge_width = 1.5; text_size=6; curvature = 0 #curvature = 0.08
colors = c('dodgerblue', 'brown3')

custom_layout <- matrix(c(
  1,2,
  1,1,
  2,1,
  2,2),
  ncol=2,byrow=TRUE)

custom_layout <- as.data.frame(custom_layout)
colnames(custom_layout) <- c("x", "y")
custom_layout$name <- c("semantic\nrules", "agreement\npatterns", "phonological\nrules", "unpredictable")

load('output/phylopath_ie.RData')

colnames(b_ci[["coef"]]) <- c("unpredictable", "phonological\nrules", "semantic\nrules", "agreement\npatterns")
rownames(b_ci[["coef"]]) <- c("unpredictable", "phonological\nrules", "semantic\nrules", "agreement\npatterns")

#creating two separate matrices so that important and unimportant paths are plotted differently (different linetypes): 
#1) setting to zero the coefficients that are crossing zero to obtain only meaningful relationships
#2) temporarily setting to zero the coefficients that are not crossing zero to obtain only non-meaningful relationships

coef_important <- as.matrix(b_ci[["coef"]])
coef_unimportant <- as.matrix(b_ci[["coef"]])
lower <- as.matrix(b_ci[["lower"]])
upper <- as.matrix(b_ci[["upper"]])

coef_important[(lower == 0 & upper == 0 ) | (lower < 0 & upper > 0 )] <- 0
coef_unimportant[(lower == 0 & upper == 0 ) | (lower < 0 & upper < 0 ) | (lower > 0 & upper > 0 )] <- 0

min <- -max(abs(b_ci[["coef"]]))
max <- max(abs(b_ci[["coef"]]))

g <- igraph::graph_from_adjacency_matrix(b_ci$coef, weighted = TRUE)
l <- ggraph::create_layout(g, 'igraph', algorithm = "sugiyama") #algorithm does not work
l$x <- custom_layout$x[match(l$name, custom_layout$name)]
l$y <- custom_layout$y[match(l$name, custom_layout$name)]

v <- colnames(b_ci$coef)
df <- as.data.frame(b_ci$coef)
df$from <- rownames(df)
df <- stats::reshape(df, varying = v, 'coef', direction = 'long')
df$to <- v[df$time]
df$lower <- c(b_ci$lower)
df$upper <- c(b_ci$upper)

df <- df[order(df$coef), ]

df <- df[!df$coef == 0,]

df$important <- ifelse((df$upper > 0 & df$lower > 0) | (df$upper < 0 & df$lower < 0), "1", "0")

arrow = grid::arrow(type = 'open', 18, grid::unit(15, 'points'))

p_ie <- ggplot2::ggplot(l) +
  ggraph::geom_edge_arc(
    ggplot2::aes_(colour = ~weight, label = ~round(weight, 2)),
    edge_width = edge_width,
    strength = curvature, arrow = arrow,
    end_cap = ggraph::rectangle(box_x, box_y, 'mm'),
    start_cap = ggraph::rectangle(box_x, box_y, 'mm'),
    show.legend = TRUE,
    linejoin = c('bevel'),
    angle_calc = 'along',
    label_dodge = grid::unit(12, 'points'), label_size=6
  ) +
  ggraph::geom_node_text(ggplot2::aes_(label = ~name), size = 6) +
  ggraph::scale_edge_color_gradient2(
    'standardized\npath coefficient',
    low = colors[1], high = colors[2],
    limits = c(min, max),
    guide = ggraph::guide_edge_colorbar()
  ) +
  ggraph::theme_graph(base_family = 'sans') + ggplot2::theme(text = element_text(size = 18), legend.text = element_text(size = 18), legend.key.size = unit(1.2, units="cm"), legend.position="bottom") +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1)) +
  ggtitle('Indo-European')
p_ie

ggsave(file="output/phylopath_ie.svg", plot=p_ie, width=6.5, height=6)


colnames(b_ci[["coef"]]) <- c("unpredictable", "phonological\nrules", "semantic\nrules", "agreement\npatterns")
rownames(b_ci[["coef"]]) <- c("unpredictable", "phonological\nrules", "semantic\nrules", "agreement\npatterns")

colnames(b_ci[["coef"]]) <- c("unpredictable", "\nphonological rules", "semantic rules", "\nagreement patterns")
rownames(b_ci[["coef"]]) <- c("unpredictable", "\nphonological rules", "semantic rules", "\nagreement patterns")

v <- colnames(b_ci$coef)
df <- as.data.frame(b_ci$coef)
df$from <- rownames(df)
df <- stats::reshape(df, varying = v, 'coef', direction = 'long')
df$to <- v[df$time]
df$lower <- c(b_ci$lower)
df$upper <- c(b_ci$upper)

df$path <- paste(df$from, df$to, sep = ' \U2192 ')

df <- df[order(df$coef), ]

df <- df[!df$coef == 0,]

coef <- ggplot2::ggplot(df,
                        ggplot2::aes_(~path, ~coef, ymin = ~lower, ymax = ~upper)) +
  ggplot2::geom_hline(yintercept = 0, size = 1, lty = 2) +
  ggplot2::geom_pointrange(size = 0.75) +
  ggplot2::xlab('') +
  ggplot2::ylab('standardized regression coefficient \U00B1 CI') + theme_classic() + theme(axis.text.x = element_text(angle = 70, hjust=1, size = 11))
coef

png("output/coefplot_phylopath_ie_custom.png", width=1800, height=1600, res = 300)
coef
dev.off()

png("output/phylopath_ie_combined.png", width=1700, height=3500, res = 300)
p / coef
dev.off()
