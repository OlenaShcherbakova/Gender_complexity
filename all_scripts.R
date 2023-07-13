#The order of running scripts
source("library.R")

#All necessary Grambank files have been already generated with the script below. One can run it only with access to the GB_wide_strict.tsv (only when the database becomes publicaly available). 
#source("generating_GB_input_file.R")

#preparing input files for running Continuous method within BayesTraits
source("preparing_files_for_BT.R") #main results
source("preparing_files_for_BT_full.R") #Supplementary Materials results

#summarizing output from Continuous method within BayesTraits

#main results
source("alphas_distribution.R")
source("correlation_distribution.R")
source("Marginal_Likelihoods_BF.R")
source("sigmas_distribution.R")
source("BT_summary_table.R")

#Supplementary Materials results with full samples (including languages that rank 0 on "semantic rules" and "agreement patterns")
source("alphas_distribution_full.R")
source("correlation_distribution_full.R")
source("Marginal_Likelihoods_BF_full.R")
source("sigmas_distribution_full.R")
source("BT_summary_table_full.R")



#phylogenetic path analysis; visualization and extracting output files
source("customizing_phylopath_plot_model_set_function.R")
source("phylopath_B.R")
source("phylopath_IE.R")
source("phylopath_world.R")

source("phylopath_summary_results.R")

#estimating phylogentic signal
source("measuring_phylosignal.R") #continuous features
source("measuring_phylosignal_discrete.R") #discrete/binary features

#plots
source("plots_phylopath.R")
source("plot_maps.R")
#source("plot_intro_map.R") this can only be run on the full GB sample
source("plot_heatmaps.R")

#robustness check: are results robust when different trees from posterior distribution are used?
source("robustness_check_tree_choice.R")