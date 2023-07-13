source('library.R')
source('alphas_distribution_mini.R.R')
source('sigmas_distribution_mini.R.R')
source('correlation_distribution_mini.R')
source('Marginal_Likelihoods_BF_mini.R')

parameters <- cbind(alphas, sigmas, correlations, likelihoods)
distributions <- parameters
distributions$BF <- 2 * (distributions$Lh_dep - distributions$Lh_indep)
distributions <- round(distributions, digits = 2)

row_names <- as.vector(results.indep.world)
row.names(distributions) <- row_names

txt <- results.dep.world

grep(".*Part Description:(.*)Installs.*", "\\1", txt)

desisered_length <- nrow(distributions)
vector <- rep(c(NA), times=desired_length)

for (i in 1:length(txt)) {
  vector[i] <- str_split(txt, "/")[[i]][4]
}

#providing a column corresponding to the file name (First trait used in the analysis, second trait used in the analysis, and Phylogeny)

distributions <- distributions %>%
  mutate(Traits = "X: Semantic rules Y: Agreement patterns",
         Phylogeny = vector)

distributions <- distributions %>% dplyr::relocate(Traits, Phylogeny)

opening_bracket <- "("
closing_bracket <- ")"

#putting lower and upper intervals into one cell
distributions$correlation_intervals <- paste (distributions$lower_correlation_dep, distributions$upper_correlation_dep, sep = " ", collapse = NULL) 
distributions$correlation_int <- paste (opening_bracket, distributions$correlation_interval, closing_bracket, sep = "", collapse = NULL) 
distributions$correlation <- paste (distributions$median_correlation_dep, distributions$correlation_int, sep = " ", collapse = NULL) 

distributions$alpha1_dep_intervals <- paste (distributions$lower_alpha1_dep, distributions$upper_alpha1_dep, sep = " ", collapse = NULL) 
distributions$alpha1_dep_int <- paste (opening_bracket, distributions$alpha1_dep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$alpha1_dep <- paste (distributions$median_alpha1_dep, distributions$alpha1_dep_int, sep = " ", collapse = NULL) 

distributions$alpha1_indep_intervals <- paste (distributions$lower_alpha1_indep, distributions$upper_alpha1_indep, sep = " ", collapse = NULL) 
distributions$alpha1_indep_int <- paste (opening_bracket, distributions$alpha1_indep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$alpha1_indep <- paste (distributions$median_alpha1_indep, distributions$alpha1_indep_int, sep = " ", collapse = NULL) 

distributions$alpha2_dep_intervals <- paste (distributions$lower_alpha2_dep, distributions$upper_alpha2_dep, sep = " ", collapse = NULL) 
distributions$alpha2_dep_int <- paste (opening_bracket, distributions$alpha2_dep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$alpha2_dep <- paste (distributions$median_alpha2_dep, distributions$alpha2_dep_int, sep = " ", collapse = NULL) 

distributions$alpha2_indep_intervals <- paste (distributions$lower_alpha2_indep, distributions$upper_alpha2_indep, sep = " ", collapse = NULL) 
distributions$alpha2_indep_int <- paste (opening_bracket, distributions$alpha2_indep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$alpha2_indep <- paste (distributions$median_alpha2_indep, distributions$alpha2_indep_int, sep = " ", collapse = NULL) 

distributions$sigma1_dep_intervals <- paste (distributions$lower_sigma1_dep, distributions$upper_sigma1_dep, sep = " ", collapse = NULL) 
distributions$sigma1_dep_int <- paste (opening_bracket, distributions$sigma1_dep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$sigma1_dep <- paste (distributions$median_sigma1_dep, distributions$sigma1_dep_int, sep = " ", collapse = NULL) 

distributions$sigma1_indep_intervals <- paste (distributions$lower_sigma1_indep, distributions$upper_sigma1_indep, sep = " ", collapse = NULL) 
distributions$sigma1_indep_int <- paste (opening_bracket, distributions$sigma1_indep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$sigma1_indep <- paste (distributions$median_sigma1_indep, distributions$sigma1_indep_int, sep = " ", collapse = NULL) 

distributions$sigma2_dep_intervals <- paste (distributions$lower_sigma2_dep, distributions$upper_sigma2_dep, sep = " ", collapse = NULL) 
distributions$sigma2_dep_int <- paste (opening_bracket, distributions$sigma2_dep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$sigma2_dep <- paste (distributions$median_sigma2_dep, distributions$sigma2_dep_int, sep = " ", collapse = NULL) 

distributions$sigma2_indep_intervals <- paste (distributions$lower_sigma2_indep, distributions$upper_sigma2_indep, sep = " ", collapse = NULL) 
distributions$sigma2_indep_int <- paste (opening_bracket, distributions$sigma2_indep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$sigma2_indep <- paste (distributions$median_sigma2_indep, distributions$sigma2_indep_int, sep = " ", collapse = NULL)

distributions_subset <- distributions %>%
  dplyr::select(Traits, Phylogeny, BF, Lh_indep, Lh_dep, alpha1_dep, alpha1_indep, alpha2_dep, alpha2_indep, sigma1_dep, sigma1_indep, sigma2_dep, sigma2_indep, correlation) %>%
  mutate(Phylogeny=recode(Phylogeny,
                          "gray_et_al2009" = "Austronesian",
                          "bouckaert_et_al2012" = "Indo-European",
                          "kolipakam_et_al2018" = "Dravidian",
                          "grollemund_et_al2015" = "Bantu"))

distributions_subset_long <- distributions_subset %>% 
  gather(Parameter, Value, c(Lh_indep:sigma2_indep)) %>% 
  separate(col = Parameter, into = c("Par", "Model"), sep = "_") %>%
  spread(key = "Par",
         value = "Value") %>%
  rename(`Marginal likelihood` = Lh) %>%
  dplyr::select(Traits, Phylogeny, Model, `Marginal likelihood`, BF, alpha1, sigma1, alpha2, sigma2, correlation)

BT_tab <- distributions_subset_long %>%
  mutate(`X alpha (95% HPD)` = alpha1,
         `Y alpha (95% HPD)`= alpha2,
         `X sigma (95% HPD)`= sigma1,
         `Y sigma (95% HPD)`= sigma2) %>%
  dplyr::select(Traits, Phylogeny, Model,`Marginal likelihood`, BF, `X alpha (95% HPD)`, `X sigma (95% HPD)`, `Y alpha (95% HPD)`, `Y sigma (95% HPD)`, correlation) %>%
  flextable() %>%
  autofit() %>%
  merge_v(j=c("Traits", "Phylogeny", "BF", "correlation")) %>%
  fix_border_issues()

latex_tab <- distributions_subset_long %>%
  mutate(`X alpha (95% HPD)` = alpha1,
         `Y alpha (95% HPD)`= alpha2,
         `X sigma (95% HPD)`= sigma1,
         `Y sigma (95% HPD)`= sigma2) %>%
  dplyr::select(Phylogeny, Model,`Marginal likelihood`, BF, `X alpha (95% HPD)`, `X sigma (95% HPD)`, `Y alpha (95% HPD)`, `Y sigma (95% HPD)`, correlation) %>%
  #flextable() %>%
  #autofit() %>%
  #merge_v(j=c("Traits", "Phylogeny", "BF", "correlation")) %>%
  kbl(caption="Continuous method results",
      format="latex") %>% #,
  kable_minimal(full_width = F) %>%
  kable_styling(latex_options = c("scale_down"))  %>%
  column_spec(1, width = "8em")

save_as_docx(
  "BayesTraits results" = BT_tab, 
  path = "output_tables/BayesTraits_results.docx")
