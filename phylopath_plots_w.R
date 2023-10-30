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

load('output/phylopath_w.RData')

colnames(a_ci[["coef"]]) <- c("unpredictable", "phonological\nrules", "semantic\nrules", "agreement\npatterns")
rownames(a_ci[["coef"]]) <- c("unpredictable", "phonological\nrules", "semantic\nrules", "agreement\npatterns")

min <- -max(abs(a_ci[["coef"]]))
max <- max(abs(a_ci[["coef"]]))

#plot of the best causal model

g <- igraph::graph_from_adjacency_matrix(a_ci[["coef"]], weighted = TRUE)
l <- ggraph::create_layout(g, 'igraph', algorithm = "sugiyama") #algorithm does not work
l$x <- custom_layout$x[match(l$name, custom_layout$name)]
l$y <- custom_layout$y[match(l$name, custom_layout$name)]

arrow = grid::arrow(type = 'open', 18, grid::unit(15, 'points'))

p_w <- ggplot2::ggplot(l) +
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
    'standardized\ncoefficient',
    low = colors[1], high = colors[2],
    limits = c(min, max),
    guide = ggraph::guide_edge_colorbar()
  ) +
  ggraph::theme_graph(base_family = 'sans') + ggplot2::theme(text = element_text(size = 18), legend.text = element_text(size = 18), legend.key.size = unit(1.2, units="cm"), legend.position="bottom") +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1)) +
  ggtitle('World')
p_w



#plots of coefficients
colnames(a_ci[["coef"]]) <- c("unpredictable", "phonological\nrules", "semantic\nrules", "agreement\npatterns")
rownames(a_ci[["coef"]]) <- c("unpredictable", "phonological\nrules", "semantic\nrules", "agreement\npatterns")

colnames(a_ci[["coef"]]) <- c("unpredictable", "phonological rules", "semantic rules", "\nagreement patterns")
rownames(a_ci[["coef"]]) <- c("unpredictable", "phonological rules", "semantic rules", "\nagreement patterns")

v <- colnames(a_ci$coef)
df <- as.data.frame(a_ci$coef)
df$from <- rownames(df)
df <- stats::reshape(df, varying = v, 'coef', direction = 'long')
df$to <- v[df$time]
df$lower <- c(a_ci$lower)
df$upper <- c(a_ci$upper)

df$path <- paste(df$from, df$to, sep = ' \U2192 ')

df <- df[order(df$coef), ]

df <- df[!df$coef == 0,]

coef_w <- ggplot2::ggplot(df,
                           ggplot2::aes_(~path, ~coef, ymin = ~lower, ymax = ~upper)) +
  ggplot2::geom_hline(yintercept = 0, size = 1, lty = "dashed", color="gray50") +
  ggplot2::geom_pointrange(size = 0.75) +
  ggplot2::xlab('') +
  #scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1)) +
  ggplot2::ylab('standardized coefficient \U00B1 CI') + theme_classic() + theme(axis.text.x = element_text(angle = 35, hjust=1, size = 18), axis.text.y = element_text(size = 18), text = element_text(size = 18), legend.text = element_text(size = 18))
coef_w

#combined plot
w <- p_w / coef_w

ggsave(file="output/phylopath_w.svg", plot=w, width=7, height=10)