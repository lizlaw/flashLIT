# flashLIT runscript

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# see further information in:
# flashLIT_info.R
# flashLIT_database_info.R
# source R code for each section
# function R code for each section

# setup ==========================================================================
# libraries
library(tidyverse)   # programming
library(scimeetr)    # grouping, identification of reading lists
library(tm)          # text mining (term document)
library(tidytext)    # tidy methods for term document matrices
library(revtools)    # grouping (optional)
library(fuzzyjoin)   # fuzzy matching of text strings
library(tidystringdist) # fuzzy matching of text strings
library(Matrix)
library(ecodist)

#library(synthesisr) # duplicates and handling of bibliography
#library(litsearchr) # automated development of new search terms (optional)
#library(bnlearn)    # mapping nodes

# custom functions 
source("./R/flashLIT_functions_import_and_clean.R")  
source("./R/flashLIT_functions_natural_clusters.R")  
source("./R/flashLIT_functions_group_summaries.R")  
source("./R/flashLIT_functions_term_list_revision.R")
source("./R/flashLIT_functions_categorise_papers.R")

# initial term lists
source("./data/flashLIT_term_lists.R")

# ggplot theme
theme_set(theme_bw() + theme(panel.grid = element_blank()))

# 0. Search term results ===========================================================

# quick analysis of the numerical results of the alternative search terms
source("./R/flashLIT_search_term_results.R") 

# 1. import data ====================================================================

# here we import the bibliography data, clean (including stem and complete), and format for the next step
# we also save the stem dictionary which we use again later to clean the term lists

# in this version we import a database containing (and distinguishing) all our sets
source("./R/flashLIT_import_sets.R")    

saveRDS(DL2, "data/flashLIT_DL2.RDS")

# DL2 <- readRDS("data/flashLIT_DL2.RDS")

source("./R/flashLIT_clean.R")  

saveRDS(stemDict, "data/flashLIT_DB1_stemDict.RDS")
saveRDS(DB1, "data/flashLIT_DB1.RDS")

# 2. cluster data =======================================================================

# clustering the data involves 
# a) developing the input data
# b) developing the distances to cluster with
# c) applying the cluster algotithm

# input data options and diatances ----------------
# 1) TAK-raw, 2) TAK-clean 3)TAK-stem
# 4) TAKKP-stem, 5) Bib, 6) TAKKP-stem-Au and 6) TAKKP-stem-Bib.
# comparing across 1:3 and 3:7 separately.
# all versions have eiclidean (_E) and bray-curtis (_BC) versions.

# DB1 <- readRDS("data/flashLIT_DB1.RDS")

source("./R/flashLIT_alternative_cluster_inputs.R")

saveRDS(list(
  TAK_raw_E, TAK_clean_E, TAK_stem_E, TAKKP_stem_E, Bib_E, TAKKP_stem_Bib_E, TAKKP_stem_Au_E,
  TAK_raw_BC, TAK_clean_BC, TAK_stem_BC, TAKKP_stem_BC, Bib_BC, TAKKP_stem_Bib_BC, TAKKP_stem_Au_BC
), "data/flashLIT_Distances.RDS")

# clustering -------------------------------------

# source("./R/flashLIT_alternative_cluster_methods.R")   *** script in progress


# older version using scimeetr like process
# source("./R/flashLIT_natural_clusters.R")  
# saveRDS(graph_DB1, "data/flashLIT_graph_DB1.RDS")
# saveRDS(coms_DB1, "data/flashLIT_coms_DB1.RDS")
# coms_DB1 %>% group_by(cl_4) %>% summarise(n=n()) %>% ggplot() + geom_histogram(aes(x=n), bins = 10)

# 3. develop group summaries ==========================================================

# here we summarise the groups - summarising the top keywords, and top papers, using different criteria.

# DB1 <- readRDS("data/flashLIT_DB1.RDS")
# graph_DB1 <- readRDS("data/flashLIT_graph_DB1.RDS")
# coms_DB1 <- readRDS("data/flashLIT_coms_DB1.RDS")

source("./R/flashLIT_group_summaries.R")  

saveRDS(coms_summaries, "data/flashLIT_coms_DB1_summmaries.RDS") 

# 4. develop and extract a list of key review papers for 'checking' results ======================

# here we search for papers with 'review' 'systematic review', 'structured review', 'meta analysis', or 'metaanalysis' in the abstract
# select papers with most number of cited papers
# and other keywords 'nest competition' 'floral cometition' 'resource competition' and 'pesticide' or 'pathogen'
# these are not meant to be complete, just an indicator of identification. 

source("./R/flashLIT_extract_key_reviews.R") 

saveRDS(KR2_cites, "data/flashLIT_KR1_cites.RDS") 
saveRDS(KR2_refs, "data/flashLIT_KR1_refs.RDS") 

# additionally ---

# tables and citations from the key reviews were extracted using pdf2complextable functions (see the BeeBN folder within)
# C:\Users\Elizabeth.Law\OneDrive - NINA\Other_work\WORKSHOPS_CONFERENCES\Workshop_EvidenceSynthesis\DataExtraction
# I need to combine this back to pdf2complextable as examples

# then process for each of the key reviews (extracting both all references and references in specific topics)

source("./R/flashLIT_extract_Mallinger.R") 
source("./R/flashLIT_extract_Russo.R") 

# in addition, we could develop a stratified sample of the different categories, and classify them manually 

# source("./R/flashLIT_stratified_sample.R") 
# write_csv(stratsamp, "data/DB1_select_for_manual_categorisation.csv")

# 5. match ===============================================================================


# 6. develop term lists further ==========================================================

# here we take the initial term lists, and look for groups that have high occurences of each term stream within them
# we then look at the top keywords from titles, abstracts, and keywords from these groups, in order to MANUALLY select additional terms
# !!! script should be run manually

# DB1 <- readRDS("data/flashLIT_DB1.RDS")
# graph_DB1 <- readRDS("data/flashLIT_graph_DB1.RDS")
# coms_DB1<- readRDS("data/flashLIT_coms_DB1.RDS")
# stemDict <- readRDS("data/flashLIT_DB1_stemDict.RDS")
# coms_summaries <- readRDS("data/flashLIT_coms_DB1_summmaries.RDS") 

source("./R/flashLIT_revise_termlists.R") 

# chosen term lists for inclusion into the next round
T2_floralcompetition <- c('behavior', 'behaviour', 'communities', 'composition', 'diet breadth', 'floral preference', 'floral resources', 'floral resources partitioning', 'floral visitation', 'flowers constancy', 'food', 'foraging', 'foraging period', 'nectar', 'nectar availability','nectar robbing', 'nectar secretion', 'nectar sugar composition', 'nectar volume', 'niche overlap', 'pollen availability', 'pollen collected', 'pollen removal', 'pollination niche', 'pollination preference', 'robbers', 'visitation', 'visitation rates', 'visitors')
T2_nestcompetition <- c('aggregation', 'boxes', 'colonies',  'competition', 'hive', 'nest', 'nest biology', 'nestboxes', 'niche', 'occupied', 'substrates', 'tubes')
T2_ppp <- c('apicystis', 'bacteria', 'bees paralysis virus', 'black queens cells virus', 'capsid', 'capsid protein vp', 'chronic', 'clothianidin', 'colonies collapse disorder', 'concentration', 'deformed wing virus',  'disease', 'ecotoxicology', 'epidemiology', 'exposed', 'exposure',  'fungicides', 'gene flow', 'genetic contamination',  'genetic pollution',  'herbicides', 'honey bees viruses', 'imidacloprid', 'immunosuppressive', 'infection',  'insecticides', 'insecticides toxicity', 'introgression',  'mite',  'neonicotinoid', 'neonicotinoid', 'neonicotinoid insecticides', 'neonicotinoid pesticides', 'nosema', 'paralysis', 'parasites', 'pathogens', 'pathogens loads', 'pathogens spillover', 'pesticides',  'pesticides residues', 'pesticides risk assess',  'predation', 'residues',  'risk assess', 'rna viruses', 'small hives beetles', 'sublethal',  'susceptibility',  'thiamethoxam', 'toxicity', 'transmission',  'trypanosome', 'varroa', 'viral', 'virion', 'virus',  'viruses')

saveRDS(T2_floralcompetition, "data/flashLIT_T2_floralcompetition.RDS")
saveRDS(T2_nestcompetition, "data/flashLIT_T2_nestcompetition.RDS")
saveRDS(T2_ppp, "data/flashLIT_T2_ppp.RDS")

# categorise groups based on prevalence of terms ==========================================================

source("./R/flashLIT_source_categorise_papers.R") 

