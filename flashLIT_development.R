# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# Causal inference informed, semi-automated systematic mapping for compex topics									

# Background ---									
# Systematic review challenged for complex topics (large numbers of interacting components, many papers, inconsistent terminology)									
# Automated grouping repeatable, but often not replicable (especially for adding new papers), and still requires interpretation of group themes and interactions.									
# Classic systematic review weighs sensitivity and specificity, narrowness of theme with volume of work. 									
# Term-list allocation problematic when incomplete list of terms known, lack of term sensitivity ot specificity for themes, lack of resources to manually resolve these.									
# More manual work = more chance for human error and bias.									
# Largely automated methods offer the benefits of truly systematic handling of literature for complex topics.									
# We propose a method that leverages the replicability of term-list allocation with automated grouping methods 									
# This can provide a living database that can be further queried and refined.									

# Aim	---								
# Develop a (semi) automated process for systematic mapping of complex topics using a modern systematic review toolbox.									

# Process	---								
# 0a	Define the 'draft' causal inference map for the topic, at least defining main interacting themes as nodes.								
# 1a	Define 'library' using broad, catch-all search terms, compile using revtools								

# 2a	Automated grouping using revtools/scimeetr (inc references, authors). First scimeetr (auto group k) then use group k within revtools. Iteratively reduce to deliver groups with <10-20 members								
# 2b	If using more than one group method, compare group membership at each level between algorithms, identify 'robust' and 'uncertain' groups	
# 2c	(optional: Extract group-level keywords, characteristics, reading lists)								

# 3a	Create a 'node-term library' from nodes & known synonyms, expert knowledge, etc. Needs to be based on the causal paths/themes.								

# 4a	Create a term-document matrix for each of the title, keyword, abstract. These should contain terms! citations could be used in early iterations, but should not be needed as existing terms should lead to groups where other potential terms can be found. 								

# 5a	Map papers with enough support onto causal map nodes based on the node-term list (should start simple, but can develop into a specificity * occurence metric, for example).								
# 5b	For each node with new papers, investigate group-level keywords (2c) to supplement (3a), then repeat process. In more advanced repetitions, use automated methods to extract keyword-combinations								

# 6a	If the group has most members in the same node/theme, consider to classify othersin the group into that node/theme								
# 6b	At last resort, identify likely themes for each group manually via keywords (revtools, scimeetr), reading lists (scimeetr)								

# 7a	Classify likely node specificity and sensitiviy (likelihood of giving good information for that node) based on classification metadata
# 7b	Allocate likelihood for unallocated groups or papers.								

# 8a	Different versions can loop through parts of the process and be documented in metadata.								

# Database structure---			

# Keys:									
# PID 	Paper ID								
# GID	Group ID								
# TID	Term ID								
# NID	Node ID		Nodes must be nested within themes, but nodes can == themes. in simplest form, 1 node is the whole library.						
# Theme	From theme list		Themes should = paths through the causal graph. 						
# RID	Researcher ID								
# AID 	Additional term ID. Additional terms are e.g. review, covariates, etc.							

# Paper details---									
# PID	Authors	Title	KW	Abstract	References	Year	Journal	Citation	DOI    (from revtools/import, need to make tools use same input)
# PID	doi exists,	date,	doi2text,	n sections,	methods extractable			(DOI check)
# PID	method	GID, level		(grouping)					
# PID	node	node certainty	node method	nodecheck method	RID	source	date		
# PID	theme	theme certainty	theme method	themecheck method	RID	source	date		

# Term document matrices---									
# PID - TID	n	(term document matrix)							
# PID - AID	n	(additional term document matrix to potentially identify e.g. lab, field, experiment, review, metaanalysis, etc.)							

# Group characteristics ---									
# GID	n	keywords	some sort of metric indicating overlap between grouping methods						
# GID	theme	RID	method	(groups identified to theme - if done)					
# GID	node	RID	method	(groups identified to node - if done)		
# GID readinglist type, PIDs allocated as this (readinglists e.g. of review, highly cited)

# Term-node identifiers	---
# theme node
# category list (fro Aterms)
# TID	term	node	node specificity	RID	method	source	date	
# AID	Aterm	category	these can be identifying e.g. review papers, or species, locations/habitats etc. Can initiate with e.g. country lists, common habitat types, common bee groups and species. 						
# 
# Node and theme characteristics	--- (created from the above, used for queries)								
# NID	number of papers, reading lists (scimeetr), number of likely review papers.								
# Theme	number of papers, reading lists (scimeetr), number of likely review papers.								

# Version metadata ---									
# VID	Description	included PID, GID, TID, AID, Theme, NID, and resulting PID:NID lists. 							

# ! There are likely to be several versions of many of these tables. Therefore dates and ways to sort and split them are required, including a 'discarded' version of them.									
# I suggest there be a default (pre-compiled) version									

# Classification model for missing PID: nodes can be based on the term-document matrix and/or group membership.									

# Ways to utilize the database ---								
# Search for papers on a specific topic, using nodes and theme information (and combinations of nodes)									
# From a known paper, identify others that are likely similar (group, keywords, terms, nodes, themes, cite, cited by)									
# Identify the number of papers attributed to each node and various node-combinations (arcs)									
# Identify nodes/arcs(i.e. node pairs)/theme knowledge (number of papers) and number of 'review' papers									
# Identify key knowledge gaps  (note importance of these will depend on bayesian network, but perhaps potential to get estimate based on graph theory).									
# Identify potential for meta-analysis for well-informed nodes (consistent, specific use of terms, number of papers, covariates, existing meta-analyses and reviews)									
# Note current doi2text coverage (doi, html access, and check of import quality, for future potential use in extracting terms from methods.									
# Methodological details to experiment with ---
# Group methods: How consistent are the groups between methods, in different themes, and with different year sets? 
# Group similarity can be: for each paper i, how many papers in Gxi are also in Gyi (group x method and group y method).
# Term-node-classifier: this is the hinge in allocating the papers to nodes. I propose by starting only with very specific node terms, even if these are broader themes. Summary could be sum or max. Threshold could be a tuning parameter. Likely more than one mention required, at least in abstract, but mention in title likely suitable by itself. 
# Group methods - search locations: include authors/citations in defining groups or not? Authors likely specialize. Citations should have many more relevant refs, but also open can of worms. 
# Term-lookup - include citations or not? I think initially not. If including will need to determine a higher thresholds, likely these depend on topics, cultures
# ** many of these might require good documentation on where terms come from. 
# Summaries can be as a network (all paths between nodes) or as pathways (allowing for all pairwise arcs along the path).

# Develop as an R package. --- (would be ideal, but for now, series of functions)
# Tools to create
# Tools to query
# Database (including term lists)

# How to have user interaction? How to have users contribute papers, and suggested nodes, and changes? ---
# Suggest main version should be 'starting data' and 'updates' available for each of the terms
# Updates can then use version number lists + updates to replicate data.


# setup -----------------------------
library(tidyverse)   # programming
library(scimeetr)    # grouping, identification of reading lists
library(tm)          # text mining (term document)
library(tidytext)    # tidy methods for term document matrices
library(revtools)   # grouping (optional)
#library(synthesisr) 
#library(litsearchr) # automated development of new search terms (optional)
#library(bnlearn)    # mapping nodes

# initiate database -----------------------------

# import data. 

# Data can come as a bibliography in many formats. We can use revtools to import the data from either .bib or .txt (.bib I think works better) 

# dd <- "./data/raw_WoS_20201105/as_bib/"
# rev.df <- revtools::read_bibliography(list.files(dd, full.names = TRUE))

# this currently fails 
#dd <- "./data/raw_WoS_20201105/as_txt/"
#revtools::read_bibliography(list.files(dd, full.names = TRUE))

# data are given revtools names:
# https://github.com/mjwestgate/revtools/blob/master/R/tag_lookup.R

# Alternatively, scimeetr imports to a specialised object, from wos ris text files ONLY, automatically removing duplicates. 
dd <- "./data/raw_WoS_20201105/as_txt/"
sci.df <- scimeetr::import_wos_files(dd)$com1$dfsci

# we can use the revtools tags to convert between formats
# ! this process can remove 'duplicate' columns

convert_names <- function(df, from, to){
  # format types = (rev)bib, ris, ris_write, wos, medline
  if(from == "rev") from <- "bib"
  if(to == "rev") to <- "bib"
  if(to == "ris") to <- "ris_write"
  
  # first, convert to bib format if not already in it (fill in duplicates with missing)
  if(!from == "bib"){
    lookup <- revtools::tag_lookup(type = from)
    names(df) <- tibble(ris = names(df)) %>% 
      left_join(lookup, by = 'ris') %>% 
      mutate(bib = map2_chr(bib, ris, ~if_else(is.na(.x), .y, .x))) %>% 
      pull(bib)
    df <- df %>% select(any_of(unname(unlist(lookup['bib']))))
  }
  
  # then, convert to other format, if not bib  
  if(!to == "bib"){
    lookup <- revtools::tag_lookup(type = to)
    names(df) <- tibble(bib = names(df)) %>% 
      left_join(lookup, by = 'bib') %>% 
      mutate(ris = map2_chr(ris, bib, ~if_else(is.na(.x), .y, .x))) %>% 
      group_by(bib) %>% 
      slice_head(n=1) %>% 
      ungroup() %>% 
      pull(ris)
    df <- df %>% select(any_of(unname(unlist(lookup['ris']))))
  }
  return(df)
}

# apply to the scimeetr import, to convert to 'revtools' format
scib.df <- sci.df %>% convert_names(from = "wos", to = "bib")

# for our DB, we can alter into a simpler document by removing all non-required columns ----------------
# note there may be some differences here depending on which source the data is being extracted from, what is available
DB0 <- scib.df %>% 
  tibble::tibble() %>% 
  dplyr::select(doi, 
                author, 
                title, 
                keywords = author_keywords, 
                keywords_plus, 
                abstract, 
                cited_references,  
                language, 
                document_type, 
                n_cited_allwos,  
                year, 
                early_access_date,
                #journal, journal_iso,
                journal_iso = source_abbreviation_29char, 
                volume, 
                pages, 
                article_number,
                date_generated
  ) %>% 
  tibble::rownames_to_column(var = 'PID') %>% 
  dplyr::mutate(PID = paste0('PID', stringr::str_pad(PID, 5, 'left', 0))) %>% 
  mutate(across(everything(), na_if, y=""))

summary(DB0)

# quick check function

quick_check <- function(df){
  tibble(
    names = df %>% names(),
    class = df %>% 
    summarise(across(everything(), class)) %>% 
    as.character(),
  na = df %>% 
    summarise(across(everything(), is.na)) %>% 
    colSums(),
  hidden_na = df %>% 
    mutate(across(everything(), na_if, y="")) %>% 
    mutate(across(everything(), na_if, y="NA")) %>% 
    summarise(across(everything(), is.na)) %>% 
    colSums()
  )
}

quick_check(DB0)

# Clean and prepare text columns for grouping -----------------------------------------------------------

# text cleaning includes removing numbers, punctuation, whitespace, meaningless words, stemming words, 
# and filling in missing keywords and keywords plus

# we want to clean the Title, Abstract, Keywords, and keywords_plus

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


# Cleaning will be applied per column. First clean, then stem.

clean_text <- function(x, sep = NULL){
  if(!is.null(sep)){
    x <- x %>% str_split(., pattern = sep) %>% .[[1]]
  }
  x <- x %>% 
    tolower() %>% 
    removeNumbers() %>% 
    removePunctuation() %>% 
    removeWords(words = meaningless_words) %>% 
    stripWhitespace() %>% 
    trimws()
  if(!is.null(sep)){
    x <- x %>% 
      paste(collapse = sep) %>% 
      na_if(y="NA")
  }
  return(x)
}

# clean_text(DB0$abstract[[1]])
# clean_text(DB0$keywords[[2]], sep = "; ")

# clean our data into a new version of the database

DB1 <- DB0 %>% 
  mutate(across(c(keywords, keywords_plus), ~map_chr(.x, clean_text, sep = "; "))) %>%    ### here we need to keep keyword bigrams together!
  mutate(across(c(title, abstract), clean_text))

quick_check(DB1)

# stemming with SnowballC::wordStem
# to be more interpretable, we would like to complete the words with the stem.
# to do that, we need to create a stem dictionary from the entire columns
# I've used ngram here, but could have used tidytext.

stem_dictionary <- function(df, ...){
  .cols = enquos(...)
  df %>% 
    select(!!!.cols) %>% 
    unite("txt", c(!!!.cols), sep = " ", na.rm = TRUE) %>% 
    pull(txt) %>% 
    toString() %>% 
    removePunctuation() %>% 
    ngram::ngram(n = 1) %>% 
    ngram::get.phrasetable() %>% 
    tibble() %>% 
    select(word = ngrams) %>% 
    mutate(word = trimws(word)) %>% 
    mutate(stem = SnowballC::wordStem(word))
}

stemDict <- stem_dictionary(DB1, title, abstract, keywords, keywords_plus)

# functions to stem, them complete using the dictionary

stemR <- function(x, token_split = NULL){
  if(!is.null(token_split)){
    x <- x %>% str_split(., pattern = token_split) %>% .[[1]]
  }
  x <- map_chr(x, ~str_split(.x, pattern = " ") %>% 
    map(SnowballC::wordStem) %>% 
    unlist() %>% 
    paste(collapse = " ")) %>% 
    na_if(y="NA")
  if(!is.null(token_split)){
    x <- x %>% 
      paste(collapse = token_split) %>% 
      na_if(y="NA")
  }
  return(x)
}

stemCompletR <- function(x, .dict, token_split = NULL) {
  if(!is.null(token_split)){
    x <- x %>% str_split(., pattern = token_split) %>% .[[1]]
  }
  x <- map(x, str_split, pattern = " ") %>% unlist(recursive = FALSE)
  x <- map_chr(x, function(.x){
    .dict[ match(.x, .dict$stem) , "word"] %>% 
    pull(word) %>% 
    paste(collapse = " ") %>% 
    unlist()  
  }) %>% 
    na_if(y="NA")
  if(!is.null(token_split)){
    x <- x %>% 
      paste(collapse = token_split) %>% 
      na_if(y="NA")
  }
  return(x)
}


# apply to our database
DB1 <- DB1 %>% 
  mutate(across(c(title, abstract), stemR))  %>% 
  mutate(across(c(keywords, keywords_plus), ~map_chr(.x, stemR, token_split = "; "))) %>%  
  mutate(across(c(title, abstract), stemCompletR, .dict = stemDict)) %>% 
  mutate(across(c(keywords, keywords_plus), ~map_chr(.x, stemCompletR, .dict = stemDict, token_split = "; ")))  

quick_check(DB1)

# for those entries without author supplied keywords, we now want to use the title and abstract to fill these in, with the 10 most frequent words used.
# I've used ngram here, but could have used tidytext.

keywords_from_columns <- function(df, .cols, .sep = "; "){
  .cols <- enquos(.cols)
  txt <- df %>% 
    select(!!!.cols) %>% 
    unite("txt", everything(), sep = " ", na.rm = TRUE) %>% 
    pull(txt) %>% 
    toString() 
  
  n0 <- ngram::ngram(txt, n = 1) %>% ngram::get.phrasetable()
  n1 <- ngram::ngram(txt, n = 1) %>% ngram::get.phrasetable() %>% filter(freq > 2)
  n2 <- ngram::ngram(txt, n = 2) %>% ngram::get.phrasetable() %>% filter(freq > 1)
  n3 <- ngram::ngram(txt, n = 3) %>% ngram::get.phrasetable() %>% filter(freq > 1)

  # first, remove the n1 rows contained in n2, then remove the n2 rows contained within n2 
  if(nrow(n2) > 0){
      n1 <- n1 %>% 
        mutate(n2 = map(n2$ngrams, str_detect, pattern = n1$ngrams) %>% 
                 as_tibble(.name_repair="unique") %>% 
                 rowSums() ) %>% 
        filter(n2==0) %>% 
      bind_rows(n2)
  }
  if(nrow(n3) > 0){
    n1 <- n1 %>% 
      mutate(n3 = map(n3$ngrams, str_detect, pattern = n1$ngrams) %>% 
               as_tibble(.name_repair="unique") %>% 
               rowSums() ) %>% 
      filter(n3==0)  %>% 
      bind_rows(n3) 
  }
  # then take to top 10 phrases from n1, filling the rest with n0
  bind_rows(
    arrange(n1, desc(freq), desc(prop)),
    arrange(n0, desc(freq), desc(prop)) 
  ) %>% 
    slice_head(n=10) %>% 
    pull(ngrams) %>% 
    trimws() %>% 
    paste(collapse = .sep) %>% 
    na_if(y="NA") %>% 
    na_if(y="")
}

# apply in a loop
for (i in which(is.na(DB1$keywords))){
  DB1$keywords[i] <- DB1[i, ] %>% keywords_from_columns(.cols = c(title, abstract))
}

quick_check(DB1)

DB1[is.na(DB0$keywords),]$keywords %>% head()

# KeyWords Plus are words or phrases that frequently appear in the titles of an article's references, but do not appear in the title of the article itself... KeyWords Plus enhances the power of cited-reference searching by searching across disciplines for all the articles that have cited references in common. (https://support.clarivate.com/ScientificandAcademicResearch/s/article/KeyWords-Plus-generation-creation-and-changes?language=en_US)
# for these, we can use the references that are both cited and in our database to generate these, as the cited ref does not contain the title.
# they will differ from the original keywords plus, but approximate them.

# match  the cited references that are in our database ---------------------------------

# the cited_references are the FirstAuthor, Year, Journal.J9, volume, pagestart, doi: for example 
# Winfree R, 2007, ECOL LETT, V10, P1105, DOI 10.1111/j.1461-0248.2007.01110.x
# this includes also papers NOT in our list. 
# we want to know which of these are in our database, 
# but there are issues here in that there are many odd cases and typos in the cited_references column
# Here we do the basics of cleaning to address this, more effort could be expended on this but guessing with diminishing returns.

# extract the cited references into a list, creating a 'CID' identifier, and retaining PID information -----

cited_references_tolist <- function(df, padCID = 6){
  x <- df %>% 
    select(fromPID = PID, cited_references) %>% 
    mutate(cited_references = map(cited_references, ~str_split(.x, pattern = "; "))) %>% 
    unnest(cited_references) %>% 
    unnest(cited_references) 

  xu <- x %>% 
    select(cited_references) %>% 
    distinct(cited_references) %>% 
    tibble::rownames_to_column(var = 'CID') %>% 
    dplyr::mutate(CID = paste0('CID', stringr::str_pad(CID, padCID, 'left', 0)))
  
  left_join(x, xu, by = 'cited_references') %>% 
    select(CID, fromPID, cited_references)
  
}

cited_references_cleanlist <- function(df){
  df %>% 
    select(-fromPID) %>% 
    distinct() %>% 
    # remove all grey reference sources starting with "(" or "*" that definitely wont be in our database
    filter(!str_starts(cited_references, pattern = "\\(")) %>% 
    filter(!str_starts(cited_references, pattern = "\\*")) %>% 
    
    # split up into name, year, Journal, Volume, Page and DOI ## option to continue improving this section in comments
    # mutate(cited_references = str_remove_all(cited_references, pattern = "\\[,") %>% str_remove_all(", \\]" )) %>% 
    separate(cited_references, into = c("author", "year", "journal_iso", "volume", "pages", "doi"), 
             sep = ", ", extra = "drop", fill = "right") %>% 
    
    # Remove everything that didn't play nicely and clean up some of the columns
    mutate(across(everything(), na_if, y="")) %>% 
    mutate(author = map_chr(author, ~str_split(.x, pattern = " ")[[1]][[1]]) %>% toupper) %>% 
    filter(!is.na(author)) %>% 
    filter(!author == "[ANONYMOUS]")  %>% 
    mutate(doi = str_replace(doi, pattern = "DOI \\[", replacement = "DOI ")) %>% 
    mutate(doi = str_replace(doi, pattern = "DOI DOI ", replacement = "DOI ")) %>% 
    filter(str_starts(volume, "V") | str_starts(doi, "DOI")) %>% 
    filter(str_starts(pages, "P")| str_starts(doi, "DOI")) 
}

CitedRefList <-  cited_references_tolist(DB1)  # Starts at almost 141K lines
CitedRefList$CID %>% unique() %>% length()     # of which there are almost 60k unique 
CitedRefClean <-  cited_references_cleanlist(CitedRefList) # reduces to 42,245 lines after cleaning (not excluding duplicates)

# extract the data base references into the same format -----

cited_reference_comparelist <- function(df){
  df %>% 
    select(PID, author, year, journal_iso, volume, pages, doi) %>% 
    mutate(author = map_chr(author, ~str_split(.x, pattern = "; ")[[1]][[1]] %>% str_remove(",")), 
           author = map_chr(author, ~str_split(.x, pattern = " ")[[1]][[1]] %>% toupper()),
           volume = paste0("V", volume),
           pages = paste0("P", pages),
           doi = map_chr(doi, ~ifelse(is.na(.x), .x, paste("DOI", str_split(.x, pattern = " ")[[1]])))) 
}

CompareList <- cited_reference_comparelist(DB1)

# Determine the matches and their acceptability ---
cited_reference_match <- function(crlist, comparelist, ..., filter_missing = TRUE){
  join_by <- enquos(...)
  x <- crlist %>% unite('.match', c(!!!join_by), sep = "_", remove = FALSE) %>% mutate(.match = na_if(.match, "NA"))
  y <- comparelist %>% unite('.match', c(!!!join_by), sep = "_", remove = FALSE) %>% mutate(.match = na_if(.match, "NA"))
  if(filter_missing){
    x <- x %>% filter(!is.na(.match))
    y <- y %>% filter(!is.na(.match))                 
  }
  inner_join(x,y, by = '.match')
}

# matches 
mDOI <- cited_reference_match(CitedRefClean, CompareList, doi)
mALL <- cited_reference_match(CitedRefClean, CompareList, author, year, journal_iso, volume, pages)
mAYVP <- cited_reference_match(CitedRefClean, CompareList, author, year, volume, pages)

# we can check some of these at this stage, for example
# confirm that the mAYVP missing from mALL are due to bad Source names
mAYVP %>% 
  filter(!(CID %in% mALL$CID)) %>% 
  select(author.x, author.y, year.x, year.y, 
         journal_iso.x, journal_iso.y, volume.x, volume.y, pages.x, pages.y)
# confirm those with mismatched DOI but matched AYVP are correctly matched  
mAYVP %>% 
  filter(!(CID %in% mDOI$CID)) %>% 
  filter(!is.na(doi.x)) %>% 
  filter(!is.na(doi.y)) %>% 
  select(doi.x, doi.y, author.x, author.y, year.x, year.y, 
         journal_iso.x, journal_iso.y)
# Brown OK - doi y from jstor
# Do OK - doi y from bioONE complete
# Godfree OK - doi x from bioONE complete
# Holzschuh OK - doi x capitalised
# Padilla OK - doi y capitalised
# Reynaldi OK - doi y has a typo (extra 0)
# Steffan-Dewenter OK - doi y has a typo (missing i)

CitedRefMatch <- bind_rows(mDOI, mALL, mAYVP) %>% 
  select(CID,  PID) %>% 
  distinct()

# this results in 1393 cited references that can be tied to the records in our database. 
# This could be refined further, and maybe get a few more, but would likely be a low benefit:cost.

# note, we can also use this as an opporutnity to fill in some of the missing DOIs

newDOI <- mAYVP %>% 
  filter(!is.na(doi.x) & is.na(doi.y)) %>% 
  select(PID, doi = doi.x) %>% 
  arrange(PID)

for (i in seq(newDOI$PID)){
  DB1$doi[which(DB1$PID == newDOI$PID[i])] <- newDOI$doi[i]
}

# To fill in the missing keywordsPlus

# apply in a loop
for (i in which(is.na(DB1$keywords_plus))){
  pid <- DB1[i, ] %>% pull(PID)
  
  # extract the cited reference list
  cr_list <- CitedRefList %>% filter(fromPID == pid)
  
  if(!is.na(cr_list$cited_references[1])){
    
      # find any that match our DB
      cr_matched <- CitedRefMatch %>% filter(CID %in% cr_list$CID)
      
      if(nrow(cr_matched)>0){
        
        # extract the details of these
        cr_details <- DB1 %>% filter(PID %in% cr_matched$PID)
        
        # create new keywords from this
        new_keyword_plus <- keywords_from_columns(cr_details, .cols = c(title))
          
        # insert into the database
         DB1$keywords_plus[i] <- new_keyword_plus
      }
    }
}

quick_check(DB1)

DB1[is.na(DB0$keywords_plus),]$keywords_plus %>% head()

# add a column that determines the number of cited references that are found in our database
DB1 <- DB1 %>% 
  mutate(number_cited_matches = map_dbl(PID, function(x){
  cr_list <- CitedRefList %>% filter(fromPID == x)
  if(!is.na(cr_list$cited_references[1])){
    cr_matched <- CitedRefMatch %>% filter(CID %in% cr_list$CID)
    y <- nrow(cr_matched)
  } else {y <- 0}
  return(y)
  })
)

DB1$number_cited_matches %>% hist()

# check wheher we could have extracted keywords_plus
DB1[is.na(DB1$keywords_plus),] %>% select(title:cited_references, number_cited_matches) %>% print(n=50)

quick_check(DB1) 

# finally, cleaning the author list. Authors are commonly given various initials, so we want to clean this to AUTHOR I, AUTHOR2 I, ...
# capitalised, surname, single initial.

clean_author_list <- function(x, .sep = "; "){
  x %>% 
    map( ~str_split(.x, pattern = .sep)[[1]]) %>% 
    map( ~toupper(.x) %>% 
         removePunctuation() %>% 
         stripWhitespace() %>% 
         trimws() %>% 
         str_extract("[A-Z]+[ ]{1}[A-Z]{1}"))
  # extract the first word (any letter any number of times) a space, and the first initial
}

DB1 <- DB1 %>% 
  mutate(author_list = clean_author_list(author))

# and we also want to split the references similarly
# here we only do the most basic of cleaning (toupper) 
# because there are many different formats
# consider cleaning through revtools matches??
DB1 <- DB1 %>% 
  mutate(cited_references_list = map(cited_references, ~str_split(.x, pattern = "; ")[[1]])) %>% 
  mutate(cited_references_list = map(cited_references_list, toupper))

quick_check(DB1) %>% print(n=21)

# end of cleaning and preparing the data
# we now have cleaned data:
# all titles, keywords, keywords_plus, and abstract are cleaned and stemmed
# and gaps of doi, keywords, keywords_plus filled in where possible
# authors and cited references split into constituent lists

saveRDS(DB1, "flashLIT_DB1.RDS")
DB1 <- readRDS("flashLIT_DB1.RDS")

# Grouping the data ------------------------------------------------------------------------------

DB1_copy <- DB1

# we want to retain one copy of the data, create a new group ID, and progressively cluster into smaller groups. 

# data going in could consist of 
# title, abstract, keywords, cited_references (i.e. reference lists), internal citations (known citation matches within our data), and authors

# There are so many different types of clustering algorithms we could use. 

# Scimeetr uses a graph clustering method (louvain or fast_greedy). It takes a graph where the nodes are papers, and the arcs are shared links (cited_references, title, abstract, keywords, authors). Titles and abstracts are composited into dtm before coupling. This coupling can weight by the type of the link, as well as the strength of the link.

# revtools uses an LDA topic model approach (from a dtm), also including bigrams where appropriate.

# the benefits of lourvain are that it determines when to stop.
# the LDA topic models need a number of topics input.

# To create a dtm using a dictionary, we:
# y <- c("text", "this")
# x <- c("this is a text", "this another one")  # this should be a clean text!
# x <- tm::VCorpus(tm::VectorSource(x))
# dtm <- tm::DocumentTermMatrix(x, control = list(dictionary = y)) 
# inspect(dtm)

# we can bring this back to a tidy format using 
# tidytext::tidy(dtm)

# DB1 %>% 
#  select(author, keywords) %>%   
#  tidytext::unnest_tokens(keyword, keywords) %>%  # unnest by "word", "ngram" (K = 2), etc
#  dplyr::count(keyword, sort = TRUE)


# basic graph method
# make edge matrices nodes x nodes showing the proportion of shared links for each element
# combine these together
# convert to igraph and cluster
# select out the cluster

# for cited references, we want the proportion intersecting, as a proportion of the union (as we want non-directional links in this case)

# Essentially we can make a dtm, then for each doc pair, take the sum of the rows, sum(min>=1)/sum(max>=1)
# this could be sped up or parrallellised https://blog.jdblischak.com/posts/pairwise-overlaps/ 

make_PID_dtm <- function(df, term_col, term_type = "list", splitpattern = "; ",  ngramlength = 2){
  term_col <- enquo(term_col)
  dtm <- df %>% 
    select(PID, term = !!term_col)
  
  if(term_type == "list"){
    dtm <- dtm %>% unnest(term) 
  } else if (term_type == "text"){
    dtm <- dtm %>% tidytext::unnest_tokens(term, term)
  } else if (term_type == "keyword_string"){
    dtm <- dtm %>% tidytext::unnest_tokens(term, term, token = stringr::str_split, pattern = splitpattern) 
  } else if (term_type == "text_ngrams"){
    dtm <- dtm %>% tidytext::unnest_tokens(term, term, token = "ngrams", n = ngramlength)
  } else {
    stop("term_type not recognised. Options are 'text', 'keyword_string', or 'text_ngrams'")
  }
  
  dtm %>% 
    group_by(PID) %>%  
    dplyr::count(term, sort = TRUE) %>% 
    pivot_wider(names_from = PID, values_from = n, values_fill = 0)
}

make_edgelist <- function(ids){
  ids %>% 
    combn(., 2) %>% 
    t() %>% 
    as.data.frame(stringsAsFactors = FALSE) %>% 
    setNames(., c("from", "to")) %>% 
    as_tibble()
}  

overlap_byEdgelist <- function(dtm, edgelist){
    edgelist %>% 
      mutate(overlap = map2_dbl(from, to, function(from.pid, to.pid){
        xy <- dtm %>% 
          select(x = matches(from.pid), y = matches(to.pid)) %>% 
          mutate(pair.overlap = map2_dbl(x, y, ~sum(.x>=1, .y>=1)))
        sum(xy$pair.overlap==2)/sum(xy$pair.overlap==1)
      })) %>% 
        pull(overlap)
}

overlap_byEdgelist_parallel <- function(dtm, edgelist){
  foreach::foreach(i=iter(edgelist, by = 'row'),
                   .packages = c("dplyr", "purrr"),
                   .combine = c) %dopar% {
  i %>% 
   map2_dbl(from, to, function(from.pid, to.pid){
      xy <- dtm %>% 
        select(x = matches(from.pid), y = matches(to.pid)) %>% 
        mutate(pair.overlap = map2_dbl(x, y, ~sum(.x>=1, .y>=1)))
      sum(xy$pair.overlap==2)/sum(xy$pair.overlap==1)
   })
                     }
}

library(iterators)
library(foreach)

title_dtm <- make_PID_dtm(DB1, title, "text") 

edgelist <- make_edgelist(DB1$PID[1:100]) 

edgelist <- edgelist %>% 
  add_column(title_overlap = overlap_byEdgelist(title_dtm, edgelist))

# this is SLOW 
# Scimeetr uses https://jangorecki.gitlab.io/-/data.table/-/jobs/640724/artifacts/public/html/data.table.html 


overlap_fromPID(DB1[1:4,], PID00001,PID00002, cited_references_list, "list")  # elements already in a list
overlap_fromPID                 # a single block of text (sep words with " ")
overlap_fromPID(DB1[1:4,], PID00001,PID00002, abstract, "text_ngrams")        # a single block of text, with ngrams (sep words with " ", into ngrams)
overlap_fromPID(DB1[1:4,], PID00001,PID00002, keywords, "keyword_string")     # a single string, elements seperated with splitpattern.


overlap_PID(DB1, PID00001,PID00002, cited_references_list, "list") # about 5 second per pair.

# OK TAKE 2 - 
# --- modify scimeetr coupling ---

scimeetr_coupling <- function(df, 
                              coupling_by = "abstract, title, keywords, keywordsplus, bibliography, authors, cr_journal", 
                              title_ngrams = 2, 
                              abstract_ngrams = 2){
  cby_abstract <- grepl(coupling_by, "abstract")
  cby_title <- grepl(coupling_by, "title")
  cby_keywords <- grepl(coupling_by, "keywords")
  cby_keywordsplus <- grepl(coupling_by, "keywordsplus")
  cby_bibliography <- grepl(coupling_by, "bibliography")
  cby_authors <- grepl(coupling_by, "authors")
  cby_cr_journal <- grepl(coupling_by, "cr_journal")
  
  df <- df %>% 
    mutate(NR = map_dbl(length(cited_references_list)))
  
  # Bibliographic coupling -------------------- here each cr is an element. we used the cleaned list.
  cr_list <- df$cited_references_list
  names(cr_list) <- df$PID
  cr_df <- data.frame('PID' = rep(names(cr_list), sapply(cr_list, length)),
                      'CR' = unlist(cr_list),
                       stringsAsFactors=F)
  couple_df <- inner_join(cr_df, cr_df, by = 'CR') %>%
    filter(PID.x > PID.y) %>%
    group_by(PID.x, PID.y) %>%
    summarise(count = n()) %>%
    left_join(df[,c('PID', 'NR')], by = c('PID.x' = 'PID'))%>%
    left_join(df[,c('PID', 'NR')], by = c('PID.y' = 'PID'))%>%
    mutate(w_ij = count/sqrt(NR.x * NR.y)) %>%
    select(PID.x, PID.y, w_ij)
  rm(cr_df)
  couple_df$w_ij[couple_df$w_ij == Inf] <- 0
  couple_df_bic <- filter(couple_df, w_ij != 0)
  
  # keyword coupling -------------------  using the author supplied keywords
  kw_list <- strsplit(df$keywords, "[;][ ]")
  names(kw_list) <- df$PID
  kw_df <- data.frame('PID'= rep(names(kw_list), sapply(kw_list, length)),
                       'KW' =  unlist(kw_list),
                       stringsAsFactors=F)
  kw_length <- group_by(kw_df, PID) %>%
    summarize(NK = n())
  couple_df <- inner_join(kw_df, kw_df, by = 'KW') %>%
    filter(PID.x > PID.y) %>%
    group_by(PID.x, PID.y) %>%
    summarise(count = n()) %>%
    left_join(kw_length, by = c('PID.x' = 'PID'))%>%
    left_join(kw_length, by = c('PID.y' = 'PID'))%>%
    mutate(w_ij = count/sqrt(NK.x * NK.y)) %>%
    select(PID.x, PID.y, w_ij)
  couple_df$w_ij[couple_df$w_ij == Inf] <- 0
  couple_df_kw <- filter(couple_df, w_ij != 0)
  
  # keyword_plus coupling ------------------- keywords drawn from the titles of the cited references
  kwp_list <- strsplit(df$keywords_plus, "[;][ ]")
  names(kwp_list) <- df$PID
  kwp_df <- data.frame('PID'= rep(names(kwp_list), sapply(kwp_list, length)),
                       'KW' =  unlist(kwp_list),
                       stringsAsFactors=F)
  kwp_length <- group_by(kwp_df, PID) %>%
    summarize(NK = n())
  couple_df <- inner_join(kwp_df, kwp_df, by = 'KW') %>%
    filter(PID.x > PID.y) %>%
    group_by(PID.x, PID.y) %>%
    summarise(count = n()) %>%
    left_join(kwp_length, by = c('PID.x' = 'PID'))%>%
    left_join(kwp_length, by = c('PID.y' = 'PID'))%>%
    mutate(w_ij = count/sqrt(NK.x * NK.y)) %>%
    select(PID.x, PID.y, w_ij)
  couple_df$w_ij[couple_df$w_ij == Inf] <- 0
  couple_df_kwp <- filter(couple_df, w_ij != 0)
  
  # Title coupling ---------------------------- this uses the pre-cleaned title. this and below needs to be modified still. need to include ngrams in this process
  documents <- tm::Corpus(tm::VectorSource(df$title))
  myTdm <- tm::DocumentTermMatrix(documents)
  myTdm2 <- tm::removeSparseTerms(myTdm, sparse = 0.99)
  dtm2list <- apply(myTdm2, 1, function(x) {
    paste(rep(names(x), x), collapse=" ")
  })
  ti_list <- strsplit(dtm2list, "[ ]")
  names(ti_list) <- df$PID
  ti_df <- data.frame('PID'= rep(names(ti_list), sapply(ti_list, length)),
                       'TI' =  unlist(ti_list),
                       stringsAsFactors=F)
  couple_df <- inner_join(ti_df, ti_df, by = 'TI') %>%
    filter(PID.x > PID.y) %>%
    group_by(PID.x, PID.y) %>%
    summarise(count = n()) %>%
    left_join(df[,c('PID', 'title')], by = c('PID.x' = 'PID'))%>%
    left_join(df[,c('PID', 'title')], by = c('PID.y' = 'PID'))%>%   ### up to here
    mutate(NR.x = str_count(TI.x, " "),
           NR.y = str_count(TI.y, " "),
           w_ij = count/sqrt(NR.x * NR.y)) %>%
    select(UT.x, UT.y, w_ij)
  couple_df$w_ij[couple_df$w_ij == Inf] <- 0
  couple_df_tic <- filter(couple_df, w_ij != 0)
  
  # Abstract coupling ---------------------------- this and below needs to be modified still. need to include ngrams in this process
  documents <- tolower(dfsci$TI)
  documents <- tm::VCorpus(tm::VectorSource(documents))
  documents <- tm::tm_map(documents, tm::removePunctuation)
  documents <- tm::tm_map(documents, tm::removeNumbers)
  meaningless_word <- c(tm::stopwords("english"), 'use', 'used', 'using', 'uses',
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
                        'increase')
  documents <- tm::tm_map(documents, tm::removeWords, meaningless_word)
  myTdm <- tm::DocumentTermMatrix(documents)
  myTdm2 <- tm::removeSparseTerms(myTdm, sparse = 0.99)
  dtm2list <- apply(myTdm2, 1, function(x) {
    paste(rep(names(x), x), collapse=" ")
  })
  TI_list <- strsplit(dtm2list, "[ ]")
  names(TI_list) <- dfsci$UT
  tiutdf <- data.frame('UT'= rep(names(TI_list), sapply(TI_list, length)),
                       'TI' =  unlist(TI_list),
                       stringsAsFactors=F)
  couple_df <- inner_join(tiutdf, tiutdf, by = 'TI') %>%
    filter(UT.x > UT.y) %>%
    group_by(UT.x, UT.y) %>%
    summarise(count = n()) %>%
    left_join(dfsci[,c('UT', 'TI')], by = c('UT.x' = 'UT'))%>%
    left_join(dfsci[,c('UT', 'TI')], by = c('UT.y' = 'UT'))%>%
    mutate(NR.x = str_count(TI.x, " "),
           NR.y = str_count(TI.y, " "),
           w_ij = count/sqrt(NR.x * NR.y)) %>%
    select(UT.x, UT.y, w_ij)
  couple_df$w_ij[couple_df$w_ij == Inf] <- 0
  couple_df_tic <- filter(couple_df, w_ij != 0)
  
  # JOC ------------------------------------------------------ here the journals came from the cited references. 
  cr_list <- strsplit(dfsci$CR, split="; ")
  names(cr_list) <- dfsci$UT
  crutdf <- data.frame('UT'= rep(names(cr_list), sapply(cr_list, length)),
                       'CR' = toupper(unlist(cr_list)),
                       stringsAsFactors=F)
  rm(cr_list)
  tmp <- left_join(crutdf, splitted_cr, by = c('CR' = 'record')) %>%
    select(UT, journal) %>%
    group_by(UT, journal) %>%
    summarise(jo_freq = n()) %>%
    filter(!is.na(journal))
  couple_df <- inner_join(tmp, tmp, by = 'journal') %>%
    filter(UT.x > UT.y) %>%
    mutate(min_jo = min(jo_freq.x, jo_freq.y)) %>%
    select(UT.x, UT.y, min_jo) %>%
    group_by(UT.x, UT.y) %>%
    summarise(count = sum(min_jo)) %>%
    left_join(dfsci[,c('UT', 'NR')], by = c('UT.x' = 'UT'))%>%
    left_join(dfsci[,c('UT', 'NR')], by = c('UT.y' = 'UT'))%>%
    mutate(w_ij = count/sqrt(NR.x * NR.y)) %>%
    select(UT.x, UT.y, w_ij)
  rm(tmp)
  gc()  # It can be useful to call gc after a large object has been removed, as this may prompt R to return memory to the operating system.
  couple_df$w_ij[couple_df$w_ij == Inf] <- 0
  couple_df_joc <- filter(couple_df, w_ij != 0)
  
  # author coupling --------------------------------------------------- need to modify this 
  cr_list <- strsplit(dfsci$CR, split="; ")
  names(cr_list) <- dfsci$UT
  crutdf <- data.frame('UT'= rep(names(cr_list), sapply(cr_list, length)),
                       'CR' = toupper(unlist(cr_list)),
                       stringsAsFactors=F)
  rm(cr_list)
  tmp <- left_join(crutdf, splitted_cr, by = c('CR' = 'record')) %>%
    select(UT, journal) %>%
    group_by(UT, journal) %>%
    summarise(jo_freq = n()) %>%
    filter(!is.na(journal))
  couple_df <- inner_join(tmp, tmp, by = 'journal') %>%
    filter(UT.x > UT.y) %>%
    mutate(min_jo = min(jo_freq.x, jo_freq.y)) %>%
    select(UT.x, UT.y, min_jo) %>%
    group_by(UT.x, UT.y) %>%
    summarise(count = sum(min_jo)) %>%
    left_join(dfsci[,c('UT', 'NR')], by = c('UT.x' = 'UT'))%>%
    left_join(dfsci[,c('UT', 'NR')], by = c('UT.y' = 'UT'))%>%
    mutate(w_ij = count/sqrt(NR.x * NR.y)) %>%
    select(UT.x, UT.y, w_ij)
  rm(tmp)
  gc()  # It can be useful to call gc after a large object has been removed, as this may prompt R to return memory to the operating system.
  couple_df$w_ij[couple_df$w_ij == Inf] <- 0
  couple_df_joc <- filter(couple_df, w_ij != 0)
  
  # Join kec, tic and bic, joc ---------------------------------------------
  couple_df <- full_join(ungroup(couple_df_bic), ungroup(couple_df_kec), by = c('UT.x', 'UT.y'), suffix = c(".bic", ".kec")) %>%
    full_join(ungroup(couple_df_tic), by = c('UT.x', 'UT.y'), suffix = c('', '.tic')) %>%
    full_join(ungroup(couple_df_joc), by = c('UT.x', 'UT.y'), suffix = c('', '.joc'))
  couple_df$w_ij.bic[is.na(couple_df$w_ij.bic)] <- 0
  couple_df$w_ij.kec[is.na(couple_df$w_ij.kec)] <- 0
  couple_df$w_ij[is.na(couple_df$w_ij)] <- 0
  couple_df$w_ij.joc[is.na(couple_df$w_ij.joc)] <- 0
  couple_df <- mutate(couple_df,
                      w_ij = 2 * w_ij.bic + w_ij.kec + w_ij + 0.5 * w_ij.joc) %>%
    select(UT.x, UT.y, w_ij)
  # Keep going normally
  m <- sum(couple_df$w_ij)
  coup2 <- data.frame(wos_id = c(as.vector(couple_df$UT.x), as.vector(couple_df$UT.y)),
                      w_ij = rep(couple_df$w_ij,2))
  k_i <- coup2 %>%
    group_by(wos_id) %>%
    summarise(k_i = sum(w_ij))
  rm(coup2)
  nam <- as.character(k_i$wos_id)
  k_i <- k_i$k_i
  names(k_i) <- nam
  couple_df <- couple_df %>%
    mutate(asso_stre = (2* w_ij * m)/ (k_i[UT.x] * k_i[UT.y]))
  rm(k_i)
  names(couple_df) <- c("rec1",
                        "rec2",
                        "bc",
                        "weight"
  )
  couple_df <- ungroup(couple_df)
  tmp <- dfsci$UT[which(!(dfsci$UT %in% unique(c(couple_df$rec2, couple_df$rec1))))]
  if(length(tmp) >= 1){
    missing_df <- data.frame('rec1' = dfsci$UT[which(!(dfsci$UT %in% unique(c(couple_df$rec2, couple_df$rec1))))],
                             'rec2' = dfsci$UT[1],
                             'bc' = 0,
                             'weight' = 0,
                             stringsAsFactors = F)
    couple_df <- rbind(couple_df, missing_df)
  }
  graph <- igraph::graph_from_data_frame(d=couple_df, directed= F)
  
}
return(graph)
}










## OLD CODE ---------------------------------------------------------------------------
# Make the groups -----------------------------
# we do this with the native scimeetr object 
# coupling by 'bickecticauc' might be better, but here we go fast
# if we can split manually, then we might want to consider parallel

### with scimeetr, the output includes nested groups, named with a "_" for each split level. 
### we can just keep splitting until the number of groups with the maximum "_" in their names (ie the last split level) is 
### equal to our group target. for other topic model options and revtools, this might need more processing. 

#' group summary of a (scimeetr) object 
#'
#' @param .x a scimeetr object
#'
#' @return a tibble with group names, level (of splitting), and size of each group.
#' @export
#'
#' @examples
grSummary <- function(.x, kw = FALSE){
  df <- tibble(
    grNames = names(.x),
    grLevel = names(.x) %>% stringr::str_count(pattern = "_"),
    grSize = map_dbl(.x, ~nrow(.x$dfsci))
  )
  if(kw) df <- df %>% 
      add_column(grTag = map(.x, ~.x$tag),
                 grKeywords = map(.x, ~.x$kw),
                 grTitle = map(.x, ~.x$ti),
                 grAbstract = map(.x, ~.x$ab))
  return(df)
}

#' iterative split of a scimeetr object into groups
#' 
#' will iteratively apply scimeetr::scimap to a scimeetr object until the stopping rules are satisfied.
#' rules are specified as both groupsize AND number of groups need to reach targets, or the end is triggered by 
#' a timeout.
#'
#' @param x a scimeetr object (can already be split)
#' @param coupling style of coupling. see ?scimeetr::scimap.
#' @param grouptarget number of groups desirable
#' @param mingroupsize minimum size of the groups
#' @param maxgroupsize maximum desired size of the groups
#' @param timeout maximum time (minutes) to iterate for
#'
#' @return
#' @export
#'
#' @examples
scimap_iterative  <- function(x, coupling = 'bickec', 
                              grouptarget = 20, mingroupsize = 5, maxgroupsize = 30,
                              timeout = 5){
  timeStart <- Sys.time()
  timeElapsed <- 0
  grSum <- grSummary(x)
  ngroups <- grSum %>% filter(grLevel == max(grLevel)) %>% nrow()
  maxGrSize <- grSum %>% filter(grLevel == max(grLevel)) %>% pull(grSize) %>% max
  
  while((ngroups  < grouptarget | maxGrSize > maxgroupsize) & !(timeElapsed > timeout)){
    
    x <- scimeetr::scimap(x, coupling_by = coupling, min_com_size = mingroupsize)
    
    grSum <- grSummary(x)
    ngroups <- grSum %>% filter(grLevel == max(grLevel)) %>% nrow()
    maxGrSize <- grSum %>% filter(grLevel == max(grLevel)) %>% pull(grSize) %>% max
    timeElapsed <- difftime(Sys.time(), timeStart, units = "mins")
  }
  attr(x,'STATUS') <- c(ngroups = ngroups, maxGrSize = maxGrSize, timeElapsed = timeElapsed)
  return(x)
}  

sciGroups <- scimap_iterative(scimeetr_list)  ## takes ~4 minuites for this data
attr(sciGroups,'STATUS') 

grSummary(sciGroups)
grSummary(sciGroups) %>% pull(grLevel) %>% max() ## goes up to 4 levels deep
grSummary(sciGroups) %>% filter(grLevel==4) %>% pull(grSize) %>% hist()

# From here, we may want to look at the keywords that define the groups (tags), by group level to help inform some of the terms
# note the full keyword, title, and abstract are the frequency of terms in these. We will use these next, in the terms

sciSum <- grSummary(sciGroups, kw = TRUE)

sciSum %>% 
  select(grNames:grTag) %>% 
  filter(grLevel <3) %>% 
  unnest_wider(grTag)

sciSum %>% 
  filter(grLevel >3) %>% 
  unnest(grTag) %>% 
  pull(grTag) %>% 
  table()

# Term document matrices  -----------------------------

# the second component we need is a list of themes, nodes, and associated terms
# here we develop this really simply, first, for terms that should split the 4 pathways, i.e. terms UNIQUE to these pathways
# we also identify several covariates - many as 'flags' that this information occurs.
# these terms will be refined, so this is just the start.

beeDictList <- list(
  terms_floralcompetition = c('floral resources', 'flower', 'nectar',
                              'floral visits', 'foraging', 'visitation rate', 'foraging period',
                              'nectar availability', 'nectar volume', 'nectar robbing',
                              'pollen availability', 'pollen collection', 
                              'niche overlap', 'floral preference', 'diet breadth',
                              'behaviour'),
  terms_nestcompetition = c('nest'),
  terms_pesticides = c('pesticide', 'LD50', 'insecticide', 'neonicotinoids', 'ecotoxicology', 'risk assessment'),
  terms_pathogen = c('mite', 'varroa', 'small hive beetle', 'bacteria', 'virus', 'transmission', 'pathogen load', 'nosema', 'apicystis', 'infection'),
  terms_predator = c('predation'),
  terms_introgression = c('introgression', 'genetic', 'mitochondrial', 'hybridization'),
  terms_plants = c('pollination', 'seed', 'fruit', 'yield', 'self-incompatability',  'self-compatability','reproductive assurance'),
  terms_es = c('ecosystem services')
)

beeDict <- tibble( type = rep(names(beeDictList), map_dbl(beeDictList, length)),
                   term = unlist(beeDictList) %>% unname(),
                   sterm = tm::stemDocument(term))

# covariates - let's ignore these for now
# terms_species_native <- c('solitary', 'semi-social', 'ground nest', 'stem nest')
# terms_species_exotic <- c('Apis mellifera', 'Bombus terrestris', 'Bombus ruderatus') 
# terms_habitat <- c('forest', 'savanna', 'savannah', 'shrubland', 'tundra', 'grassland', 'wetland', 'desert', 'arable', 'pasture', 'plantation', 'crop', 'rural', 'urban', 'residential', 'degraded forest', 'grazing', 'karst', 'rocky', 'coastal')
# terms_climate <- c('boreal', 'subarctic', 'subantarctic', 'temperate', 'subtropical', 'tropical', 'dry', 'wet', 'lowland', 'montane', 'mediterranean', 'continental', 'precipitation', 'temperature')
# terms_outcomes <- c('abundance', 'weight', 'width', 'mass', 'size', 'fitness', 'reproduc', 'diversity', 'density')
# terms_researchcontext <- c('experiment', 'lab', 'field', 'observation', 'transect', 'pan trap', 'net', )
# terms_countries <- raster::ccodes()[,"NAME"]
# terms_regions <- raster::ccodes()[,"UNREGION1"]
# terms_continents <- terms_regions <- raster::ccodes()[,"continent"]

## to create these, we need to use tm to create a corpus, and them apply the term document matrix using a supplied dictionary of terms. 
## to create the corpus ready for this we need to clean it (lower case, remove whitespace, and stem words)  then the dictionary we supply also needs to have stemmed words.
## we can convert back to tidy format with broom::tidy and broom::cast

tm::stemDocument()
tm::stemCompletion()
inspect(DocumentTermMatrix(reuters, list(dictionary = c("prices", "crude", "oil"))))


# PID - TID	n	(term document matrix)							
# PID - AID	n	(additional term document matrix to potentially identify e.g. lab, field, experiment, review, metaanalysis, etc.)							

# Group characteristics ---									
# GID	n	keywords	some sort of metric indicating overlap between grouping methods						
# GID	theme	RID	method	(groups identified to theme - if done)					
# GID	node	RID	method	(groups identified to node - if done)		
# GID readinglist type, PIDs allocated as this (readinglists e.g. of review, highly cited)

# Term-node identifiers	---
# theme node
# category list (fro Aterms)
# TID	term	node	node specificity	RID	method	source	date	
# AID	Aterm	category	these can be identifying e.g. review papers, or species, locations/habitats etc. Can initiate with e.g. country lists, common habitat types, common bee groups and species. 						
# 
# Node and theme characteristics	--- (created from the above, used for queries)								
# NID	number of papers, reading lists (scimeetr), number of likely review papers.								
# Theme	number of papers, reading lists (scimeetr), number of likely review papers.		