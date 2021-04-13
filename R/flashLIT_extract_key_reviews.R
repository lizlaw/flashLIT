# flashLIT identify review papers

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

# we want these papers as a list of papers for checking the results
# here we search for papers with 'review' 'systematic review', 'structured review', 'meta analysis', or 'metaanalysis' in the abstract
# select papers with most number of cited papers
# and other keywords 'nest competition' 'floral cometition' 'resource competition' and 'pesticide' or 'pathogen'
# these are not meant to be complete, just an indicator of the group quality

allreviews <- DB1 %>% 
  filter(grepl(pattern ='review|metaanalysis|meta analysis', x = abstract))
# 143 papers from 2400

compreviews <- allreviews %>% 
  filter(grepl(pattern ='competition', x = abstract))
# 22 papers

floralcompreviews <- compreviews %>% 
  filter(grepl(pattern ='floral|flower|food', x = abstract)) %>% 
  pull(ref)
  #pull(abstract_orig)
# 16 papers

nestcompreviews <- compreviews %>% 
  filter(grepl(pattern ='nest', x = abstract)) %>% 
  pull(ref)
 #pull(abstract_orig)
# 5 papers

pppreviews <- allreviews %>% 
  filter(grepl(pattern ='pesticide|insecticide|pathogen|predator|bacteria|virus|disease', x = abstract)) %>%  # 50 papers
  filter(grepl(pattern ='systematic', x = abstract)) %>%   # 7 papers
  pull(ref)
#pull(abstract_orig)

opppreviews <- allreviews %>% 
  filter(grepl(pattern ='pathogen|predator|bacteria|virus', x = abstract)) %>%  # 22 papers
  pull(ref)
#pull(abstract_orig)


# general, multiple categories
"MALLINGER R, et al., 2017, Do managed bees have negative effects on wild bees?: A systematic review of the literature, PLOS ONE"      
"RUSSO L, 2016, Positive and Negative Impacts of Non-Native Bee Species around the World, INSECTS"
"PAINI D, 2004, Impact of the introduced honey bee (Apis mellifera) (Hymenoptera : Apidae) on native bees: A review, AUSTRAL ECOL"
"INOUE M, et al., 2010, Competition for flower resources and nest sites between Bombus terrestris (L.) and Japanese native bumblebees, APPL ENTOMOL ZOOL"
# foral competition specific
"WOJCIK V, et al., 2018, Floral Resource Competition Between Honey Bees and Wild Bees: Is There Clear Evidence and Can We Guide Management and Conservation?, ENVIRON ENTOMOL"
# ppp specific
"LUNDIN O, et al., 2015, Neonicotinoid Insecticides and Their Impacts on Bees: A Systematic Review of Research Approaches and Identification of Knowledge Gaps, PLOS ONE"
"ARENA M, et al., 2014, A meta-analysis comparing the sensitivity of bees to pesticides, ECOTOXICOLOGY"
"CULLEN M, et al., 2019, Fungicides, herbicides and bees: A systematic review of existing research and methods, PLOS ONE" 
"FUNFHAUS A, et al., 2018, Bacterial pathogens of bees, CURR OPIN INSECT SCI" 
"EVISON S, et al., 2018, The biology and prevalence of fungal diseases in managed and wild bees, CURR OPIN INSECT SCI" 
"MANLEY R, et al., 2015, Emerging viral disease risk to pollinating insects: ecological, evolutionary and anthropogenic factors, J APPL ECOL"

# other papers that could be useful
"RODRIGUEZGIRONES M, et al., 2012, Effects of body size and sociality on the anti-predator behaviour of foraging bees, OIKOS" 
"CANE J, et al., 2011, Predicted fates of ground-nesting bees in soil heated by wildfire: Thermal tolerances of life stages and a survey of nesting depths, BIOL CONSERV" 
"WOOD T, et al., NA, Managed honey bees as a radar for wild bee decline?, APIDOLOGIE"  
"OLGUN T, et al., 2020, Comparative analysis of viruses in four bee species collected from agricultural, urban, and natural landscapes, PLOS ONE" 


# for each of these papers, we can download their bibliographies (or the lists within a paper), and identify the papers allocated to each category
# these can be used to test the categorisation.

list.files("data/KeyReferences")

# What we would like to split on -----------------------------

# Nest competition (i.e. specifically citing nest site competition)
# Floral competition (inc visitation rates, foraging behaviour, floral resources, and managed bee driven changes in pollinator networks)
# Pathogens, pesticides, etc (relative sensitivity)
# Plant mediated (i.e. plant driven changes in pollinator networks)
# No mechanism suggested (change in presence, abundance, diversity/composition, fitness (including size, survival and nest-based fitness)

# displacement and niche overlap is likely to indicate floral competition


