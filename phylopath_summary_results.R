#generating the table out of phylopath results across Indo-European, Bantu, and world tree

ie <- read.csv("output_tables/Phylopath_output_table_IE.csv")
b <- read.csv("output_tables/Phylopath_output_table_B.csv")
w <- read.csv("output_tables/Phylopath_output_table_W.csv")

all <- as.data.frame(rbind(ie, b, w)) %>%
  dplyr::select(-1)

colnames(all) <- gsub("\\.", " ", colnames(all))


all_latex <- all %>%
  kbl(caption="Phylogenetic path analysis on three phylogenies",
      format="latex") %>% #,
  kable_minimal(full_width = F) %>%
  kable_styling(latex_options = c("scale_down")) %>% #font_size = 4 for saving space
  column_spec(1, width = "22em")
