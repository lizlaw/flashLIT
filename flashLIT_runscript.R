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

# 1. import data ====================================================================

# here we import the bibliography data, clean (including stem and complete), and format for the next step
# we also save the stem dictionary which we use again later to clean the term lists

source("./R/flashLIT_import_and_clean.R")            

saveRDS(stemDict, "data/flashLIT_DB1_stemDict.RDS")
saveRDS(DB1, "data/flashLIT_DB1.RDS")

# 2. cluster data =======================================================================

# here we cluster the data, by first developing a graph of the documents, with edges as the similarity between documents
# in this version, we weight the abstract similarity more than the title and keywords, and these more than the authors, cited references, and cited journals
# clutering is automated to iterate until all groups are below a threshold size

# DB1 <- readRDS("data/flashLIT_DB1.RDS")

source("./R/flashLIT_natural_clusters.R")  

saveRDS(graph_DB1, "data/flashLIT_graph_DB1.RDS")
saveRDS(coms_DB1, "data/flashLIT_coms_DB1.RDS")

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
# these are not meant to be complete, just an indicator of the group quality

source("./R/flashLIT_extract_key_reviews.R") 

# then process for each of the key reviews (extracting both all references and references in specific topics)

source("./R/flashLIT_extract_Mallinger.R") 
source("./R/flashLIT_extract_Russo.R") 

# in addition, we can develop a stratified sample of the different categories, and classify them manually 

source("./R/flashLIT_stratified_sample.R") 
write_csv(stratsamp, "data/DB1_select_for_manual_categorisation.csv")


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

possible_floralcompetition_words 
# 'adaptation', 'agricultural', 'alien', 'alpine', 'ancient woodland', 'anthidium manicatum', 'ants', 'apidae', 'apis', 'apis mellifera', 'argentine', 'arvensis', 'attracted', 'b', 'bees', 'bees vision', 'behavior', 'behaviour', 'biological', 'biological invasive', 'biology', 'birds pollination', 'bombus', 'bombus impatiens', 'bombus terrestris', 'border', 'british', 'bumblebees', 'butterflies', 'characteristics', 'chemistry', 'choice', 'clarkia', 'coastal', 'colour', 'communities', 'composition', 'consequences', 'conservation', 'crop', 'deceptive', 'diet breadth', 'discrimination', 'distribution', 'double', 'ecology', 'effective', 'environment', 'epipactis', 'evaluated', 'evolution', 'exotic', 'experimental', 'experiments', 'facilitation', 'flies', 'flora', 'floral', 'floral biology', 'floral evolution', 'floral preference', 'floral resources', 'floral resources partitioning', 'floral scent', 'floral traits', 'floral visitation', 'flowers', 'flowers', 'flowers colour', 'flowers constancy', 'food', 'foraging', 'foraging', 'foraging period', 'fruit', 'fruit set', 'fynbos', 'gardens', 'general', 'germination', 'globulus', 'guides', 'habitat', 'honey', 'honeybees', 'hoverflies', 'illumination', 'insect', 'insect pollination', 'introduced', 'invasive', 'invasive plant', 'iridaceae', 'islands', 'isles', 'l', 'landscape', 'lapeirousia', 'leguminosae', 'length', 'lessons', 'linepithema humile', 'linked', 'lysimachia', 'lythrum salicaria', 'mascula', 'mechanisms', 'mimulus', 'model', 'multimodal', 'myrtaceae', 'native', 'natural', 'nectar', 'nectar', 'nectar availability', 'nectar guides', 'nectar robbing', 'nectar secretion', 'nectar sugar composition', 'nectar volume', 'nest', 'niche overlap', 'occurred', 'officinalis', 'ophrys heldreichii', 'orchid', 'orchidaceae', 'orchis', 'outcrossing', 'p', 'palustris', 'patches', 'patterns', 'phenotypic selection', 'plant', 'plantpollinator interactions', 'pollen', 'pollen availability', 'pollen collected', 'pollen removal', 'pollination', 'pollination effective', 'pollination efficiency', 'pollination niche', 'pollination preference', 'pollinatorfriendly', 'polymorphism', 'populations', 'preference', 'presented', 'prey', 'produced', 'production', 'range', 'rates', 'reciprocal transplanted', 'relationship', 'reproductive', 'resources', 'responses', 'rewards', 'robbers', 'role', 'season', 'secondary', 'seed', 'selection', 'set', 'sexual', 'shape', 'signals', 'solitary bees', 'south', 'southern africa', 'species', 'studied', 'subspecies', 'success', 'sugar', 'system', 'terrestris', 'traits', 'treatments', 'trees', 'tropical', 'tuber', 'urban', 'variation', 'vegetation', 'visitation', 'visitation rates', 'visitors', 'wild', 'wild flowers', 'xantiana', 'years'

possible_nestcompetition_words 
# 'activity', 'adult', 'africanized honey bees', 'aggregation', 'agricultural', 'alfalfa', 'alfalfa leafcutting bees', 'almond', 'analis', 'apidae', 'apoidea', 'areas', 'artificial nest', 'bees', 'biology', 'birds', 'blue', 'blue orchards bees', 'boxes', 'breeding', 'caatinga', 'cavities', 'cells', 'centridini', 'centris', 'centris analis', 'cerrado', 'chalkbrood', 'cherry', 'cm', 'colonies', 'commercial', 'competition', 'conservation', 'crop pollination', 'density', 'dispersal', 'distribution', 'female', 'field', 'floral', 'floral host', 'floral resources', 'flowers', 'foraging', 'fruit set', 'habitat', 'halictidae', 'honey bees', 'hymenoptera', 'hymenopterans', 'interactions', 'landscape', 'leafcutting', 'lignaria', 'lipid', 'm', 'macaw', 'management', 'mason', 'medicago sativa', 'megachile', 'megachile rotundata', 'megachilidae', 'melanderi', 'native', 'native bees', 'natural enemies', 'nest', 'nest', 'nest biology', 'nestboxes', 'niche', 'nomia', 'north', 'o', 'observed', 'occupied', 'oilcollecting', 'orchards', 'orchards pollination', 'osmia', 'osmia bicornis', 'osmia lignaria', 'osmia rufa', 'parasitoids', 'plant', 'pollen', 'pollination', 'populations', 'predation', 'preference', 'previously', 'prunus cerasus', 'red mason bees', 'rediviva', 'relationship', 'reproductive', 'resources', 'retention', 'rotundata', 'savanna', 'sites', 'soil', 'solitary', 'solitary bees', 'species', 'states', 'substrates', 'success', 'suggest', 'taxa', 'toads', 'trapnests', 'traps nest', 'tubes', 'vertebrates', 'wasps', 'whether'

possible_ppp_words 
# 'acute', 'adult', 'agricultural', 'apicystis', 'apidae', 'apis', 'apis mellifera', 'apis mellifera mellifera', 'applications', 'ascosphaera', 'ascosphaera apis', 'assess', 'b', 'bacteria', 'bacteria', 'beekeeping', 'bees', 'bees apis mellifera', 'bees conservation', 'bees declines', 'bees health', 'bees paralysis virus', 'behaviour', 'bicornis', 'biosynthesis', 'black', 'black queens cells virus', 'bombi', 'bombus', 'bombus terrestris', 'bqcv', 'bumble', 'bumble bees', 'bumblebees', 'butterflies', 'capsid', 'capsid protein vp', 'cause', 'cells', 'cerana', 'chronic', 'clothianidin', 'colonies', 'colonies collapse disorder', 'combined', 'commercial', 'compared', 'concentration', 'conservation', 'corn', 'crop', 'crop production', 'dark', 'de', 'declines', 'deformed', 'deformed wing virus', 'detected', 'development', 'differential', 'disease', 'dispersal', 'diversity', 'domain', 'dwv', 'e', 'ecosystem services', 'ecotoxicology', 'ecotoxicology', 'emergence', 'environmental', 'eovata', 'epidemiology', 'eucalypt', 'eucalyptus', 'european', 'evidence', 'exotic', 'exotic species', 'exposed', 'exposure', 'f', 'field', 'fieldrealistic', 'fitness', 'flow', 'foraging', 'form', 'fungicides', 'gene', 'gene flow', 'genetic', 'genetic', 'genetic contamination', 'genetic diversity', 'genetic pollution', 'genome', 'globulus', 'habitat', 'herbicides', 'history', 'hokkaido', 'honey', 'honey bees', 'honey bees conservation', 'honey bees viruses', 'honeybees', 'host', 'hybridisation', 'hybridization', 'hybridization', 'hymenoptera', 'iapv', 'imidacloprid', 'immunosuppressive', 'impatiens', 'infection', 'infection', 'insect', 'insect pollination', 'insecticides', 'insecticides', 'insecticides toxicity', 'interactions', 'introduced', 'introgression', 'introgression', 'invasive', 'islands', 'l', 'la', 'landscape', 'landscape genetic', 'ld', 'learning', 'levels', 'lineages', 'm', 'management', 'markers', 'mason bees', 'melipona', 'mellifera', 'microsatellite', 'mite', 'mitochondrial', 'mitochondrial', 'mitochondrial dna', 'mitochondrial genome', 'n', 'native', 'negative', 'neonicotinoid', 'neonicotinoid', 'neonicotinoid insecticides', 'neonicotinoid pesticides', 'nitens', 'north', 'nosema', 'nosema', 'nosema cerana', 'orchid', 'osmia', 'osmia bicornis', 'ovata', 'paradigm', 'paralysis', 'parasites', 'paternity', 'pathogens', 'pathogens loads', 'pathogens spillover', 'pesticides', 'pesticides', 'pesticides residues', 'pesticides risk assess', 'plant', 'plantations', 'pollen', 'pollen dispersal', 'pollenmediated', 'pollination', 'pollination declines', 'populations', 'predation', 'produced', 'production', 'protection', 'protein', 'quadrifasciata', 'queens', 'rates', 'reared', 'recent', 'reduced', 'reducedrisk', 'replicated', 'reproductive', 'residues', 'review', 'risk', 'risk assess', 'risk assess', 'rna', 'rna viruses', 'sbpv', 'sbv', 'seedlings', 'sensitivity', 'sequence', 'services', 'small hives beetles', 'snps', 'social', 'solitary', 'solitary bees', 'spartina', 'species', 'spinosad', 'stingless', 'stressors', 'structure', 'sublethal', 'subspecies', 'sunflower', 'survival', 'susceptibility', 'synergism', 'synergistic', 'terrestris', 'thiamethoxam', 'toxicity', 'traits', 'transmission', 'trees', 'trypanosome', 'variation', 'varroa', 'viral', 'virion', 'virus', 'virus', 'viruses', 'vp', 'western', 'western honey bees', 'wide', 'wild', 'wild bees', 'wild pollination', 'wing', 'workers'

T2_floralcompetition <- c('behavior', 'behaviour', 'communities', 'composition', 'diet breadth', 'floral preference', 'floral resources', 'floral resources partitioning', 'floral visitation', 'flowers constancy', 'food', 'foraging', 'foraging period', 'nectar', 'nectar availability','nectar robbing', 'nectar secretion', 'nectar sugar composition', 'nectar volume', 'niche overlap', 'pollen availability', 'pollen collected', 'pollen removal', 'pollination niche', 'pollination preference', 'robbers', 'visitation', 'visitation rates', 'visitors')
T2_nestcompetition <- c('aggregation', 'boxes', 'colonies',  'competition', 'hive', 'nest', 'nest biology', 'nestboxes', 'niche', 'occupied', 'substrates', 'tubes')
T2_ppp <- c('apicystis', 'bacteria', 'bees paralysis virus', 'black queens cells virus', 'capsid', 'capsid protein vp', 'chronic', 'clothianidin', 'colonies collapse disorder', 'concentration', 'deformed wing virus',  'disease', 'ecotoxicology', 'epidemiology', 'exposed', 'exposure',  'fungicides', 'gene flow', 'genetic contamination',  'genetic pollution',  'herbicides', 'honey bees viruses', 'imidacloprid', 'immunosuppressive', 'infection',  'insecticides', 'insecticides toxicity', 'introgression',  'mite',  'neonicotinoid', 'neonicotinoid', 'neonicotinoid insecticides', 'neonicotinoid pesticides', 'nosema', 'paralysis', 'parasites', 'pathogens', 'pathogens loads', 'pathogens spillover', 'pesticides',  'pesticides residues', 'pesticides risk assess',  'predation', 'residues',  'risk assess', 'rna viruses', 'small hives beetles', 'sublethal',  'susceptibility',  'thiamethoxam', 'toxicity', 'transmission',  'trypanosome', 'varroa', 'viral', 'virion', 'virus',  'viruses')

saveRDS(T2_floralcompetition, "data/flashLIT_T2_floralcompetition.RDS")
saveRDS(T2_nestcompetition, "data/flashLIT_T2_nestcompetition.RDS")
saveRDS(T2_ppp, "data/flashLIT_T2_ppp.RDS")

# categorise groups based on prevalence of terms ==========================================================

source("./R/flashLIT_source_categorise_papers.R") 

