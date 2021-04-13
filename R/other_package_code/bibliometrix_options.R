# bibliometrix version

# install.packages("bibliometrix")
library(tidyverse)
library(bibliometrix)

# upload bib files. Here the txt version doesnt work, and so we use the bib version.
## during the process several additional columns are coded, we can add more
dd <- list.files("./data/raw_WoS_20201105/as_bib/", full.names = TRUE)
M <- convert2df(dd, dbsource = "wos", format = "bibtex")
M <- metaTagExtraction(M, Field = "AU_CO", sep = ";")

# 'correct' PY 
## as the data for early access articles do not have a year, 
## we sub in the early access date, as we don't want to loose these from the analysis
## probably also want to convert all this way.
# M %>% select(TI, PY, early.access.date) %>% filter(is.na(PY))
M <- M %>% 
  mutate_at("PY", ~ pmap_dbl(list(M$PY, M$early.access.date), 
                          function(py,ea) ifelse(!is.na(py), py, ea %>% str_extract('[0-9]{4}') %>% 
                                                   as.numeric()))) 

# Descriptive analysis - simple tallies and plots can be viewed:
results <- biblioAnalysis(M, sep = ";")
S <- summary(object = results, k = 10, pause = FALSE)
plot(x = results, k = 10, pause = FALSE)

# networks
## these create a network by adding connections that count number of co-occurences of (authors, citations, keywords)
## scaled relatively - so as not to benefit those with many and penalise those with few.

# country collaboration network
NetMatrix <- biblioNetwork(M, analysis = "collaboration", network = "countries", sep = ";")
net <- networkPlot(NetMatrix, n = dim(NetMatrix)[1], Title = "Country Collaboration", 
                   type = "circle", size=TRUE, remove.multiple=FALSE,labelsize=0.8)

# co-citation network
## NetMatrix <- biblioNetwork(M, analysis = "co-citation", network = "references", sep = ";")
## net=networkPlot(NetMatrix, n = 30, Title = "Co-Citation Network", type = "fruchterman", 
##                size=T, remove.multiple=FALSE, labelsize=0.7,edgesize = 5)

# Keyword co-occurrences
NetMatrix <- biblioNetwork(M, analysis = "co-occurrences", network = "keywords", sep = ";")
net=networkPlot(NetMatrix, normalize="association", weighted=T, n = 30, Title = "Keyword Co-occurrences", 
                type = "fruchterman", size=T,edgesize = 5,labelsize=0.7)

# co word analysis - conceptual structure of a field
##  *conceptualStructure* that performs a CA or MCA to draw a conceptual structure of the field and
## K-means clustering to identify clusters of documents which express common concepts
## includes natural language processing (NLP) routines (see the function *termExtraction*) to extract terms from titles and abstracts.
## usefully, i think there is scope for adding specific remove terms and also keeping terms (compoind words) that can be put in 
## (the latter should not matter for clustering however). Also synonyms It can implement the Porter's stemming algorithm to reduce inflected 
## (or sometimes derived) words to their word stem, base or root form. k can be specified, or "auto"

CS <- conceptualStructure(M, field="ID", method="MCA", minDegree=4, clust=4 ,k.max=8, stemming=FALSE, labelsize=10, documents=10)

# HISTORICAL DIRECT CITATION NETWORK
## The historiographic map ...a chronological network map of most relevant direct citations resulting from a bibliographic collection. 
histResults <- histNetwork(M, sep = ";")
net <- histPlot(histResults, size = 10)