# flashLIT term lists

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# =============================

# meaningless words are a combination of scimeetr and revtools
meaningless_words <- c(tm::stopwords("english"), 
                       # scimeetr
                       'use', 'used', 'using', 'uses',
                       'new', 'effect', 'effects', 'affect', 'affects', 'impact',
                       'impacts', 'implication', 'implications', 'potential',
                       'influence', 'influences', 'influenced', 'study', '-',
                       'data', 'can', 'results', 'different', 'similar', 'also',
                       'c', 'may', 'based', 'important', 'within','however',
                       'found', 'analysis', 'changes', 'among', 'large',
                       'number', 'higher', 'well', 'studies', 'total',
                       'increased', 'increases', 'elsevier', 'level', 'many',
                       'rights', 'present', 'will', 'low', 'across', 'showed',
                       'associated', 'approach', 'related', 'provide', 'including',
                       'increase',
                       # revtools
                       "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
                       "eleven", "twelve", "thirteen", "fourteen", "fifteen",
                       "sixteen", "seventeen", "eighteen", "nineteen",
                       "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety",
                       "hundred", "thousand", "million", "billion", "trillion",
                       "first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth",
                       "eleventh", "twelfth", "thirteenth", "fourteenth", "fifteenth",
                       "sixteenth", "seventeenth", "eighteenth", "nineteenth",
                       "twentieth", "thirtieth", "fortieth", "fiftieth", "sixtieth",
                       "seventieth", "eightieth", "ninetieth",
                       "hundredth", "thousandth", "millionth", "billionth")

# =============================

# themes, nodes, and associated terms
# here we develop this really simply, first, for terms that should split the 4 pathways, i.e. terms UNIQUE to these pathways
# we also identify several covariates - many as 'flags' that this information occurs.
# these terms will be refined, so this is just the start.

  terms_floralcompetition <- c('floral resources', 'flower', 'nectar',
                              'floral visits', 'foraging', 'visitation rate', 'foraging period',
                              'nectar availability', 'nectar volume', 'nectar robbing',
                              'pollen availability', 'pollen collection', 
                              'niche overlap', 'floral preference', 'diet breadth',
                              'behaviour')
  terms_nestcompetition <- c('nest')
  terms_pesticides <- c('pesticide', 'LD50', 'insecticide', 'neonicotinoids', 'ecotoxicology', 'risk assessment')
  terms_pathogen <- c('mite', 'varroa', 'small hive beetle', 'bacteria', 'virus', 'transmission', 'pathogen load', 'nosema', 'apicystis', 'infection')
  terms_predator <- c('predation')
  terms_introgression <- c('introgression', 'genetic', 'mitochondrial', 'hybridization')


# covariates - let's ignore these for now ----------------------------------
# terms_plants <- c('pollination', 'seed', 'fruit', 'yield', 'self-incompatability',  'self-compatability','reproductive assurance')
# terms_es <- c('ecosystem services')
# terms_species_native <- c('solitary', 'semi-social', 'ground nest', 'stem nest')
# terms_species_exotic <- c('Apis mellifera', 'Bombus terrestris', 'Bombus ruderatus') 
# terms_habitat <- c('forest', 'savanna', 'savannah', 'shrubland', 'tundra', 'grassland', 'wetland', 'desert', 'arable', 'pasture', 'plantation', 'crop', 'rural', 'urban', 'residential', 'degraded forest', 'grazing', 'karst', 'rocky', 'coastal')
# terms_climate <- c('boreal', 'subarctic', 'subantarctic', 'temperate', 'subtropical', 'tropical', 'dry', 'wet', 'lowland', 'montane', 'mediterranean', 'continental', 'precipitation', 'temperature')
# terms_outcomes <- c('abundance', 'weight', 'width', 'mass', 'size', 'fitness', 'reproduc', 'diversity', 'density')
# terms_researchcontext <- c('experiment', 'lab', 'field', 'observation', 'transect', 'pan trap', 'net', )
# terms_countries <- raster::ccodes()[,"NAME"]
# terms_regions <- raster::ccodes()[,"UNREGION1"]
# terms_continents <- terms_regions <- raster::ccodes()[,"continent"]
