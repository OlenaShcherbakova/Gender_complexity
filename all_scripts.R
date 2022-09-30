#The order of running scripts

#Continuous method within BayesTraits


#phylogenetic path analysis
source("phylopath_B.R")
source("phylopath_IE.R")
source("phylopath_world.R")

source("phylopath_summary_results.R")

#estimating phylogentic signal
source("measuring_phylosignal.R")
source("measuring_phylosignal_discrete.R")


#plots
source("plots_phylopath.R")
source("plot_maps.R")
source("plot_heatmaps.R")