#library.R file

library(ape)
library(bayestraitr)
library(tidyverse)
library(dplyr)
library(phytools)
library(phylopath)
library(here)
library(viridis)
library(geiger)
library(ggraph)
library(ggplot2)
library(RColorBrewer)
library(ggtree)
library(ggnewscale)
library(caper)
library(patchwork)
library(rethinking)
library(stringr)
library(flextable)
library(kableExtra)



# helper function to convert and load trait data into a named vector for
# plotting.
get_trait_vector <- function(tree, data, variable) {
  x <- data[tree$tip.label, variable]
  x[is.na(x)] <- 0  # set NA's to zero to enable plotting.
  names(x) <- tree$tip.label
  x
}

# Function to load trees from DPLACE-data repository
# can be used either with or without renameto = 'glottocode'
load_trees <- function(dirname, type='posterior', mappingfile='taxa.csv', renameto=NA) {
  # check file type
  if (type == 'summary') {
    treefile <- file.path(dirname, 'summary.trees')
  }
  else if (type == 'posterior') {
    treefile <- file.path(dirname, 'posterior.trees')
  } else {
    stop(paste("Unknown Tree Type:", type))
  }
  
  # check file exists
  if (file.exists(treefile) == FALSE) {
    stop(paste("Invalid file:", treefile))
  }
  
  trees <- ape::read.nexus(treefile)
  if (class(trees) == 'phylo') { trees <- c(trees) ; class(trees) <- 'multiPhylo' }
  
  # make full path if just given taxa.csv
  if (mappingfile == 'taxa.csv') { mappingfile <- file.path(dirname, mappingfile) }
  
  if (file.exists(mappingfile) & is.na(renameto) == FALSE) {
    mapping <- read.csv(mappingfile, header = TRUE, stringsAsFactors = FALSE, na.string="")
    
    # check the required columns exist
    if ('taxon' %in% colnames(mapping) == FALSE) stop(paste('column `taxon` not in', mappingfile))
    
    if (renameto %in% colnames(mapping) == FALSE) stop(paste('colname', renameto, 'not in', mappingfile))
    
    trees <- ape::.uncompressTipLabel(trees)
    
    for (i in 1:length(trees)){
      # remove tips not in `taxon` mapping
      missing <- trees[[i]]$tip.label[trees[[i]]$tip.label %in% mapping[['taxon']] == FALSE]
      if (length(missing) > 0) {
        trees[[i]] <- ape::drop.tip(trees[[i]], missing)
      }
      
      # remove tips not in `renameto` mapping
      missing <- mapping[is.na(mapping[[renameto]]), 'taxon']
      if (length(missing) > 0) {
        trees[[i]] <- ape::drop.tip(trees[[i]], missing)
      }
      
      # handle duplicate rename tips
      dupes <- mapping[duplicated(mapping[[renameto]], incomparables=NA), ]
      if (nrow(dupes)) {
        warning(paste("Removing ", nrow(dupes), "tips that will be duplicated after rename:", paste(dupes[['taxon']], collapse=", ")))
        trees[[i]] <- ape::drop.tip(trees[[i]], dupes[['taxon']])
      }
      
      # rename tips
      matches <- match(trees[[i]]$tip.label, mapping[['taxon']])
      trees[[i]]$tip.label <- mapping[matches, renameto]
    }
    trees <- ape::.compressTipLabel(trees, ref=mapping[matches, renameto])
  }
  trees
}



load_data_intro_figure <- function(filename="data/GB_input.tsv") {
  grambank <- read.csv(filename, header = TRUE, sep = '\t', stringsAsFactors=FALSE)
  colnames(grambank)[colnames(grambank)=="Language_ID"] <- "Glottocode"
  #grambank[is.na(grambank)] <- "-"
  grambank_compl <- subset(x = grambank, select = c("Glottocode", "GB030", "GB051", "GB052", "GB053", "GB054", "GB170", "GB171", "GB172", "GB177", "GB192", "GB198", "GB314", "GB315", "GB321"))
  grambank_compl <- na.omit(grambank_compl)
  
  for(i in 1:nrow(grambank_compl)){
    summ <- sum(c(as.numeric(grambank_compl$GB030[i]), as.numeric(grambank_compl$GB051[i]), as.numeric(grambank_compl$GB052[i]), as.numeric(grambank_compl$GB053[i]), as.numeric(grambank_compl$GB054[i]), grambank_compl$GB170[i]), as.numeric(grambank_compl$GB172[i]), as.numeric(grambank_compl$GB177[i]), as.numeric(grambank_compl$GB192[i]), as.numeric(grambank_compl$GB198[i]), as.numeric(grambank_compl$GB314[i]), as.numeric(grambank_compl$GB315[i]), as.numeric(grambank_compl$GB321[i]), na.rm = T)
    if(summ > 0 ){grambank_compl$Gender[i] <- 1}
    else(grambank_compl$Gender[i] <- 0)
  }
  
  grambank_compl <- subset(x = grambank_compl, select = c("Glottocode", "Gender"))
  
  # make sure we have factors here.
  for (col in colnames(grambank_compl)) {
    grambank_compl[[col]] <- as.factor(grambank_compl[[col]])
  }
  rownames(grambank_compl) <- grambank_compl$Glottocode
  grambank_compl
}


load_data_final_short <- function(filename="data/GB_input.tsv") {
  grambank <- read.csv(filename, header = TRUE, sep = '\t', stringsAsFactors=FALSE)
  colnames(grambank)[colnames(grambank)=="Language_ID"] <- "Glottocode"
  #grambank[is.na(grambank)] <- "-"
  grambank_compl <- subset(x = grambank, select = c("Glottocode", "GB030", "GB051", "GB052", "GB053", "GB054", "GB170", "GB171", "GB172", "GB192", "GB198", "GB321"))
  grambank_compl <- na.omit(grambank_compl)
  
  
  for(i in 1:nrow(grambank_compl)){
    summ <- sum(c(as.numeric(grambank_compl$GB051[i]), as.numeric(grambank_compl$GB052[i]), as.numeric(grambank_compl$GB053[i]), as.numeric(grambank_compl$GB054[i])), na.rm = T)
    if(summ > 0 ){grambank_compl$sem_classes[i] <- summ/4}
    else(grambank_compl$sem_classes[i] <- 0)
  }
  
  
  for(i in 1:nrow(grambank_compl)){
    summ <- sum(c(as.numeric(grambank_compl$GB030[i]), as.numeric(grambank_compl$GB170[i]), as.numeric(grambank_compl$GB171[i]), as.numeric(grambank_compl$GB172[i]), as.numeric(grambank_compl$GB198[i])), na.rm = T)
    if(summ > 0 ){grambank_compl$agr_patterns[i] <- summ/5}
    else(grambank_compl$agr_patterns[i] <- 0)
  }
  
  
  grambank_compl <- subset(x = grambank_compl, select = c("Glottocode", "sem_classes", "agr_patterns"))
  
  # make sure we have factors here.
  for (col in colnames(grambank_compl)) {
    grambank_compl[[col]] <- as.factor(grambank_compl[[col]])
  }
  rownames(grambank_compl) <- grambank_compl$Glottocode
  grambank_compl
}

OUTPUTDIR_output<- here("output")		
# create output dir if it does not exist.		
if (!dir.exists(OUTPUTDIR_output)) { dir.create(OUTPUTDIR_output) }
