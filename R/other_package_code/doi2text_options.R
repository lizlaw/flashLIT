# doi2text
library(tidyverse)
library(scimeetr)
library(doi2txt)

# load the scimeetr data, and get the urbanisation subgroup
scisub_bktj30 <- readRDS("data/raw_WoS_20201105/scisub_bktj30.RDS")
sci111 <- dive_to(scisub_bktj30, aim_at = "com1_1_1")

# exploring which types of paper have no doi
sci111[[1]]$dfsci %>% 
  select(DI, D2, PT, DT, PY, TI) %>% 
  filter(DI == "")
## these seem still relevant - would ideally find the DOI post hoc

## screeen for language type and document type, then extract the doi
doidf <- sci111[[1]]$dfsci %>% 
  select(DI) %>% 
  filter(DI != "") %>% 
  as_tibble()

doidf <- doidf %>% 
  mutate(article = map(DI, function(x) try(doi2html(x)))) %>% 
  mutate(abstract = map(article, function(x) try(extract_section(x, section = "abstract", 
                                                                 max_lines = 1, clean=TRUE, min_words = 50, forcestart = TRUE)))) %>% 
  mutate(methods = map(article, function(x) try(extract_section(x, section = "methods", 
                                                                 max_lines = 50, clean=TRUE, min_words = 50, forcestart = FALSE))))
