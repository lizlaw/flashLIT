# flashLIT revise term lists

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

# one of the challenges in systematic reviews is term lists - specifically knowing appropriate term lists to start.
# this process can help to revise them

# we begin with basic term lists imported via flashLIT_term_lists
# the main categories we will focus on here are defining:
# floral competition, nest competition, and ppp (pesticides, pathogens, predators and introgression) pathways

# ==================================================================
# first, we need to align these with our cleaned data, and understand how these might be transformed.

# transform terms into clean, stemmed, completed:
# note if not in dictionary (ie not found in original) will return NA for the completed

T_floralcompetition <- term_clean_stem_complete(terms_floralcompetition)
T_nestcompetition <- term_clean_stem_complete(terms_nestcompetition)
T_ppp <- term_clean_stem_complete(c(terms_pesticides, terms_pathogen, terms_predator, terms_introgression))

T_floralcompetition %>% print(n=nrow(.))
T_nestcompetition %>% print(n=nrow(.))
T_ppp %>% print(n=nrow(.))

# from this, we can determine all of our terms are in the dictionary, are reasonably cleaned, and will be found at least once in our data.

# ==================================================================

# this was my attempt at fuzzy lookup terms. not ideal.
# we also want to look up the words that will be selected by the different word choices from the original data
# to do this we first need to create a lookup dictionary from the columns we want (and clean but not stem)
# then we want to 'fuzzy' lookup the words to see which others may be similar and suitable

# note this is not perfect - only takes ngrams from the keywords if they are ngrams 
TAK_dict <- lookupDict(DB1, .cols = 'title, abstract, keywords', ngrams = 1:2)
TAKO_dict <- lookupDict(DB1, .cols = 'title_orig, abstract_orig, keywords_orig', ngrams = 1:2)


# see options for stringdist methods ?stringdist 
# typos = jw, 
# The Optimal String Alignment distance (method='osa') is like the Levenshtein distance (counts the number of deletions, insertions and substitutions necessary to turn b into a) but also allows transposition of adjacent characters. this is the default.
# The full Damerau-Levenshtein distance (method='dl') is like the optimal string alignment distance except that it allows for multiple edits on substrings.
# The Jaro distance (method='jw', p=0), is a number between 0 (exact match) and 1 (completely dissimilar) measuring dissimilarity between strings. 

fuzzyjoin::stringdist_left_join(tibble(word = T_floralcompetition$completed), TAK_dict, 
                by = "word",
                method = "jw", 
                max_dist = 0.2, 
                distance_col = "dist") %>% print(n=nrow(.))

fuzzyjoin::stringdist_left_join(tibble(word = T_floralcompetition$completed), TAK_dict, 
                                by = "word",
                                method = "osa", 
                                max_dist = 3, 
                                distance_col = "dist") %>% print(n=nrow(.))

# not ideal... discard this track for now.

# ==================================================================

# adding to our term lists via examining other keywords used by groups with high prevalence of our keywords

# identify groups with high occurence of specific terms
# i.e. at least one of the group terms in most papers. 

# explore how this works with one example: 
# term_group_proportions(df = DB1, colnames = c(title, abstract, keywords),
#                        groupname = "2_2_3", communities = coms_summaries, 
#                        dict = T_floralcompetition$completed, type = "documents")
#     
# coms_summaries %>% filter(group_id == "2_2_3") %>% pull(group_top_keywords)

# apply to all groups - here we only select title, abstract, and author-supplied keywords.

coms_summaries <- coms_summaries %>% 
  mutate(pDocs_floralcompetition_in_TAK = map_dbl(group_id, ~term_group_proportions(df = DB1, colnames = c(title, abstract, keywords),
                                                                             groupname = .x, communities = coms_summaries, 
                                                                             dict = T_floralcompetition$completed, type = "documents"))) %>% 
  mutate(pDocs_nestcompetition_in_TAK = map_dbl(group_id, ~term_group_proportions(df = DB1, colnames = c(title, abstract, keywords),
                                                                                    groupname = .x, communities = coms_summaries, 
                                                                                    dict = T_nestcompetition$completed, type = "documents"))) %>% 
  mutate(pDocs_ppp_in_TAK = map_dbl(group_id, ~term_group_proportions(df = DB1, colnames = c(title, abstract, keywords),
                                                                                    groupname = .x, communities = coms_summaries, 
                                                                                    dict = T_ppp$completed, type = "documents"))) %>% 
  ungroup()
  

# visualise what proportions there are for the different groups ---

ggplot(coms_summaries) +
  geom_vline(data = coms_summaries[1,], aes(xintercept = pDocs_nestcompetition_in_TAK, color = group_level)) +
  geom_histogram(aes(x = pDocs_floralcompetition_in_TAK, color = group_level, fill = group_level), 
                 alpha = 0.2, position = "identity", binwidth = 0.05)

ggplot(coms_summaries) +
  geom_vline(data = coms_summaries[1,], aes(xintercept = pDocs_nestcompetition_in_TAK, color = group_level)) +
  geom_histogram(aes(x = pDocs_nestcompetition_in_TAK, color = group_level, fill = group_level), 
                 alpha = 0.2, position = "identity", binwidth = 0.05)

ggplot(coms_summaries) +
  geom_vline(data = coms_summaries[1,], aes(xintercept = pDocs_nestcompetition_in_TAK, color = group_level)) +
  geom_histogram(aes(x = pDocs_ppp_in_TAK, color = group_level, fill = group_level), 
                 alpha = 0.2, position = "identity", binwidth = 0.05)

# from this, we can see that we might need to refine the floral competition stream further

# we want to select the communities for which the proportion of documents including the group terms is high (but want to also include some flexibility to identify useful other keywords)

gFC <- coms_summaries %>% filter(pDocs_floralcompetition_in_TAK > 0.95)
gNC <- coms_summaries %>% filter(pDocs_nestcompetition_in_TAK > 0.5)
gPPP <- coms_summaries %>% filter(pDocs_ppp_in_TAK > 0.8)

# then from these we pull out the top keywords and summarise ==============================================
# manually selecting appropriate keywords
# remember, here we are using the cleaned data, so some of the keywords will not necessarily be gramatically correct.
# we are wanting words that will likely be unique to the category. This might need some revision of these terms below.

# create a string of this with our existing list to select from
possible_floralcompetition_words <- suggest_terms(gFC, T_floralcompetition$completed, keywords, title_words, abstract_words)
possible_nestcompetition_words <- suggest_terms(gNC, T_nestcompetition$completed, keywords, title_words, abstract_words)
possible_ppp_words <- suggest_terms(gPPP, T_ppp$completed, keywords, title_words, abstract_words)



possible_floralcompetition_words 
# 'adaptation', 'agricultural', 'alien', 'alpine', 'ancient woodland', 'anthidium manicatum', 'ants', 'apidae', 'apis', 'apis mellifera', 'argentine', 'arvensis', 'attracted', 'b', 'bees', 'bees vision', 'behavior', 'behaviour', 'biological', 'biological invasive', 'biology', 'birds pollination', 'bombus', 'bombus impatiens', 'bombus terrestris', 'border', 'british', 'bumblebees', 'butterflies', 'characteristics', 'chemistry', 'choice', 'clarkia', 'coastal', 'colour', 'communities', 'composition', 'consequences', 'conservation', 'crop', 'deceptive', 'diet breadth', 'discrimination', 'distribution', 'double', 'ecology', 'effective', 'environment', 'epipactis', 'evaluated', 'evolution', 'exotic', 'experimental', 'experiments', 'facilitation', 'flies', 'flora', 'floral', 'floral biology', 'floral evolution', 'floral preference', 'floral resources', 'floral resources partitioning', 'floral scent', 'floral traits', 'floral visitation', 'flowers', 'flowers', 'flowers colour', 'flowers constancy', 'food', 'foraging', 'foraging', 'foraging period', 'fruit', 'fruit set', 'fynbos', 'gardens', 'general', 'germination', 'globulus', 'guides', 'habitat', 'honey', 'honeybees', 'hoverflies', 'illumination', 'insect', 'insect pollination', 'introduced', 'invasive', 'invasive plant', 'iridaceae', 'islands', 'isles', 'l', 'landscape', 'lapeirousia', 'leguminosae', 'length', 'lessons', 'linepithema humile', 'linked', 'lysimachia', 'lythrum salicaria', 'mascula', 'mechanisms', 'mimulus', 'model', 'multimodal', 'myrtaceae', 'native', 'natural', 'nectar', 'nectar', 'nectar availability', 'nectar guides', 'nectar robbing', 'nectar secretion', 'nectar sugar composition', 'nectar volume', 'nest', 'niche overlap', 'occurred', 'officinalis', 'ophrys heldreichii', 'orchid', 'orchidaceae', 'orchis', 'outcrossing', 'p', 'palustris', 'patches', 'patterns', 'phenotypic selection', 'plant', 'plantpollinator interactions', 'pollen', 'pollen availability', 'pollen collected', 'pollen removal', 'pollination', 'pollination effective', 'pollination efficiency', 'pollination niche', 'pollination preference', 'pollinatorfriendly', 'polymorphism', 'populations', 'preference', 'presented', 'prey', 'produced', 'production', 'range', 'rates', 'reciprocal transplanted', 'relationship', 'reproductive', 'resources', 'responses', 'rewards', 'robbers', 'role', 'season', 'secondary', 'seed', 'selection', 'set', 'sexual', 'shape', 'signals', 'solitary bees', 'south', 'southern africa', 'species', 'studied', 'subspecies', 'success', 'sugar', 'system', 'terrestris', 'traits', 'treatments', 'trees', 'tropical', 'tuber', 'urban', 'variation', 'vegetation', 'visitation', 'visitation rates', 'visitors', 'wild', 'wild flowers', 'xantiana', 'years'

possible_nestcompetition_words 
# 'activity', 'adult', 'africanized honey bees', 'aggregation', 'agricultural', 'alfalfa', 'alfalfa leafcutting bees', 'almond', 'analis', 'apidae', 'apoidea', 'areas', 'artificial nest', 'bees', 'biology', 'birds', 'blue', 'blue orchards bees', 'boxes', 'breeding', 'caatinga', 'cavities', 'cells', 'centridini', 'centris', 'centris analis', 'cerrado', 'chalkbrood', 'cherry', 'cm', 'colonies', 'commercial', 'competition', 'conservation', 'crop pollination', 'density', 'dispersal', 'distribution', 'female', 'field', 'floral', 'floral host', 'floral resources', 'flowers', 'foraging', 'fruit set', 'habitat', 'halictidae', 'honey bees', 'hymenoptera', 'hymenopterans', 'interactions', 'landscape', 'leafcutting', 'lignaria', 'lipid', 'm', 'macaw', 'management', 'mason', 'medicago sativa', 'megachile', 'megachile rotundata', 'megachilidae', 'melanderi', 'native', 'native bees', 'natural enemies', 'nest', 'nest', 'nest biology', 'nestboxes', 'niche', 'nomia', 'north', 'o', 'observed', 'occupied', 'oilcollecting', 'orchards', 'orchards pollination', 'osmia', 'osmia bicornis', 'osmia lignaria', 'osmia rufa', 'parasitoids', 'plant', 'pollen', 'pollination', 'populations', 'predation', 'preference', 'previously', 'prunus cerasus', 'red mason bees', 'rediviva', 'relationship', 'reproductive', 'resources', 'retention', 'rotundata', 'savanna', 'sites', 'soil', 'solitary', 'solitary bees', 'species', 'states', 'substrates', 'success', 'suggest', 'taxa', 'toads', 'trapnests', 'traps nest', 'tubes', 'vertebrates', 'wasps', 'whether'

possible_ppp_words 
# 'acute', 'adult', 'agricultural', 'apicystis', 'apidae', 'apis', 'apis mellifera', 'apis mellifera mellifera', 'applications', 'ascosphaera', 'ascosphaera apis', 'assess', 'b', 'bacteria', 'bacteria', 'beekeeping', 'bees', 'bees apis mellifera', 'bees conservation', 'bees declines', 'bees health', 'bees paralysis virus', 'behaviour', 'bicornis', 'biosynthesis', 'black', 'black queens cells virus', 'bombi', 'bombus', 'bombus terrestris', 'bqcv', 'bumble', 'bumble bees', 'bumblebees', 'butterflies', 'capsid', 'capsid protein vp', 'cause', 'cells', 'cerana', 'chronic', 'clothianidin', 'colonies', 'colonies collapse disorder', 'combined', 'commercial', 'compared', 'concentration', 'conservation', 'corn', 'crop', 'crop production', 'dark', 'de', 'declines', 'deformed', 'deformed wing virus', 'detected', 'development', 'differential', 'disease', 'dispersal', 'diversity', 'domain', 'dwv', 'e', 'ecosystem services', 'ecotoxicology', 'ecotoxicology', 'emergence', 'environmental', 'eovata', 'epidemiology', 'eucalypt', 'eucalyptus', 'european', 'evidence', 'exotic', 'exotic species', 'exposed', 'exposure', 'f', 'field', 'fieldrealistic', 'fitness', 'flow', 'foraging', 'form', 'fungicides', 'gene', 'gene flow', 'genetic', 'genetic', 'genetic contamination', 'genetic diversity', 'genetic pollution', 'genome', 'globulus', 'habitat', 'herbicides', 'history', 'hokkaido', 'honey', 'honey bees', 'honey bees conservation', 'honey bees viruses', 'honeybees', 'host', 'hybridisation', 'hybridization', 'hybridization', 'hymenoptera', 'iapv', 'imidacloprid', 'immunosuppressive', 'impatiens', 'infection', 'infection', 'insect', 'insect pollination', 'insecticides', 'insecticides', 'insecticides toxicity', 'interactions', 'introduced', 'introgression', 'introgression', 'invasive', 'islands', 'l', 'la', 'landscape', 'landscape genetic', 'ld', 'learning', 'levels', 'lineages', 'm', 'management', 'markers', 'mason bees', 'melipona', 'mellifera', 'microsatellite', 'mite', 'mitochondrial', 'mitochondrial', 'mitochondrial dna', 'mitochondrial genome', 'n', 'native', 'negative', 'neonicotinoid', 'neonicotinoid', 'neonicotinoid insecticides', 'neonicotinoid pesticides', 'nitens', 'north', 'nosema', 'nosema', 'nosema cerana', 'orchid', 'osmia', 'osmia bicornis', 'ovata', 'paradigm', 'paralysis', 'parasites', 'paternity', 'pathogens', 'pathogens loads', 'pathogens spillover', 'pesticides', 'pesticides', 'pesticides residues', 'pesticides risk assess', 'plant', 'plantations', 'pollen', 'pollen dispersal', 'pollenmediated', 'pollination', 'pollination declines', 'populations', 'predation', 'produced', 'production', 'protection', 'protein', 'quadrifasciata', 'queens', 'rates', 'reared', 'recent', 'reduced', 'reducedrisk', 'replicated', 'reproductive', 'residues', 'review', 'risk', 'risk assess', 'risk assess', 'rna', 'rna viruses', 'sbpv', 'sbv', 'seedlings', 'sensitivity', 'sequence', 'services', 'small hives beetles', 'snps', 'social', 'solitary', 'solitary bees', 'spartina', 'species', 'spinosad', 'stingless', 'stressors', 'structure', 'sublethal', 'subspecies', 'sunflower', 'survival', 'susceptibility', 'synergism', 'synergistic', 'terrestris', 'thiamethoxam', 'toxicity', 'traits', 'transmission', 'trees', 'trypanosome', 'variation', 'varroa', 'viral', 'virion', 'virus', 'virus', 'viruses', 'vp', 'western', 'western honey bees', 'wide', 'wild', 'wild bees', 'wild pollination', 'wing', 'workers'

# endscript