#generating an input file used in the repository from the entire Grambank dataset
library(here)

OUTPUTDIR_data<- here("data")		
# create output dir if it does not exist.		
if (!dir.exists(OUTPUTDIR_data)) { dir.create(OUTPUTDIR_data) }

grambank <- read.csv("data/GB_wide_strict.tsv", header = TRUE, sep = '\t', stringsAsFactors=FALSE)
colnames(grambank)[colnames(grambank)=="Language_ID"] <- "Glottocode"
grambank_compl <- subset(x = grambank, select = c("Glottocode", "GB030", "GB051", "GB052", "GB053", "GB054", "GB170", "GB171", "GB172", "GB192", "GB198", "GB321"))
grambank_compl <- na.omit(grambank_compl)
  
grambank_compl %>%
  write_tsv(file="data/GB_input.tsv")