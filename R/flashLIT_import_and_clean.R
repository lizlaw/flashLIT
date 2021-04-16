# flashLIT import and clean

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

# Search from WoS - multiple different versions searched on 15.04.2021, Searching 1987-current (all available)

#M: Managed	(African* XXX bee) OR Apis OR Bombus OR "bumble bee" OR bumblebee* OR "honey bee" OR honeybee OR ((introduc* OR invas* OR non-native OR nonnative OR commercial OR exotic OR feral OR managed) XXX (bee OR pollin*))
#N: Native	(((cavity OR ground OR stem) XXX nesting) OR (native OR solitary OR wild)) XXX (bee OR pollin*)
#I: Interaction	pollinat* OR network* OR niche OR partition* OR interact* OR competit* OR facilitat* OR mutualis* OR “resource limitation” OR hybridization OR introgression OR dependen* OR assemblag* OR overlap OR spillover OR impact* OR community OR diversity

# where XXX is AND, NEAR (defaults to NEAR/15), NEAR/3, and NEAR/1

# also, we search either TOPIC, or TITLE + ABSTRACT + KEYWORDS (TAK)

# we summarise with AND, for M+N and M+N+I

# we are not interested in all these groups, some are slightly larger

list.files("data/raw_WoS_20210415")

# here, we develop the grouping using the groups initially using the smallest group
# MNI_n1_TAK


# ==================================================================

# import data. 
# Data can come as a bibliography in many formats, and can be imported via revtools or scimeetr.
# Revtools imports from either .bib or .txt (.bib I think works better). data are given revtools names:
# https://github.com/mjwestgate/revtools/blob/master/R/tag_lookup.R
# Alternatively, scimeetr imports to a specialised object, from wos or scopus ris text files ONLY, automatically removing duplicates. 

# revtools bib
# dd <- "./data/raw_WoS_20201105/as_bib/"
# rev.df <- revtools::read_bibliography(list.files(dd, full.names = TRUE))

# revtools txt - currently fails 
# dd <- "./data/raw_WoS_20201105/as_txt/"
# revtools::read_bibliography(list.files(dd, full.names = TRUE))

# scimeetr import 
# dd <- "./data/raw_WoS_20201105/as_txt/"
# sci.df <- scimeetr::import_wos_files(dd)$com1$dfsci

# customised scimeetr import (straight to df with file and folder information)
dd <- "./data/raw_WoS_20210415/MNI_n1_TAK/"
sci.df <- import_wos_files_DFonly(dd)

# convert the scimeetr import to 'revtools' format
scib.df <- sci.df %>% convert_names(from = "wos", to = "bib")

# choose columns that we will use (note choice may be dictated by import type) and create PID
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
                date_generated,
                folder,
                file
  ) %>% 
  dplyr::mutate(
    title_orig = title,
    abstract_orig = abstract,
    keywords_orig = keywords
    ) %>% 
  tibble::rownames_to_column(var = 'PID') %>% 
  dplyr::mutate(PID = paste0('PID', stringr::str_pad(PID, 5, 'left', 0))) %>% 
  dplyr::mutate(across(everything(), na_if, y=""))

# check current status
summary(DB0)
quick_check(DB0) %>% print(n=nrow(.))

# clean, stem, fill dataframe ===================================================================

# text cleaning includes removing numbers, punctuation, whitespace, meaningless words 
# we want to clean the Title, Abstract, Keywords, and keywords_plus - into a new version of the database

DB1 <- DB0 %>% 
  mutate(across(c(keywords, keywords_plus), ~map_chr(.x, clean_text, sep = "; "))) %>%    ### here we need to keep keyword phrases together!
  mutate(across(c(title, abstract), clean_text))

quick_check(DB1) %>% print(n=nrow(.))

# stemming means that similar terms will be made comparable. here, we re-complete words (via the created dictionary) in order for the words to be understandable by humans.

stemDict <- stem_dictionary(DB1, title, abstract, keywords, keywords_plus)

DB1 <- DB1 %>% 
  mutate(across(c(title, abstract), stemR))  %>% 
  mutate(across(c(keywords, keywords_plus), ~map_chr(.x, stemR, token_split = "; "))) %>%  
  mutate(across(c(title, abstract), stemCompletR, .dict = stemDict)) %>% 
  mutate(across(c(keywords, keywords_plus), ~map_chr(.x, stemCompletR, .dict = stemDict, token_split = "; ")))  

quick_check(DB1) %>% print(n=nrow(.))

# for those entries without author supplied keywords, we now want to use the title and abstract to fill these in, with the 10 most frequent words used.
# apply in a loop to fill in the missing keywords

for (i in which(is.na(DB1$keywords))){
  DB1$keywords[i] <- DB1[i, ] %>% keywords_from_columns(.cols = c(title, abstract))
}

quick_check(DB1) %>% print(n=nrow(.))

DB1[is.na(DB0$keywords),]$keywords %>% head()

# KeyWords Plus are words or phrases that frequently appear in the titles of an article's references, but do not appear in the title of the article itself... KeyWords Plus enhances the power of cited-reference searching by searching across disciplines for all the articles that have cited references in common. (https://support.clarivate.com/ScientificandAcademicResearch/s/article/KeyWords-Plus-generation-creation-and-changes?language=en_US)
# for these, we can use the references that are both cited and in our database to generate these, as the cited ref does not contain the title.
# they will differ from the original keywords plus, but approximate them.

# so, we first need to identify all those cited references that match

# extract the cited references into a list, creating a 'CID' identifier, and retaining PID information 
CitedRefList <-  cited_references_tolist(DB1)  # Starts at almost 80K lines
CitedRefList$CID %>% unique() %>% length()     # of which there are almost 37k unique 
CitedRefClean <-  cited_references_cleanlist(CitedRefList) # reduces to 25,854 lines after cleaning

# extract a comparison list from the documents
CompareList <- cited_reference_comparelist(DB1) 

# matches based on DOI, all columns, and all but DOI
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

# Do OK - doi y from bioONE complete
# Holzschuh OK - doi x capitalised
# Steffan-Dewenter OK - doi y has a typo (missing i)

CitedRefMatch <- bind_rows(mDOI, mALL, mAYVP) %>% 
  select(CID, PID) %>% 
  distinct()

# this results in 727 cited references that can be tied to the records in our database, approx 54% of our database, and only 3% of all the unique (cleaned) bibliography entries.
# This could be refined further, and maybe get a few more, but would likely be a low benefit:cost.

# Add a column which includes the cited references that can be matched (i.e. PID's of matched cited references) in each PID
# join fromPID to the matched ones
CitedMatched <- CitedRefMatch %>% 
  left_join(CitedRefList, by = "CID") %>% 
  select(citedPID = PID, PID = fromPID) %>% 
  group_by(PID) %>% 
  summarise(citedPID_list = list(citedPID))

DB1 <- DB1 %>% 
  left_join(CitedMatched, by = "PID")

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

quick_check(DB1) %>% print(n=nrow(.))

DB1[is.na(DB0$keywords_plus),]$keywords_plus %>% head()

# add a column that determines the number of cited references that are found in our database
# this could now be done by counting the list entries in the citedPID_list
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

quick_check(DB1)  %>% print(n=nrow(.))

# cleaning the author list. Authors are commonly given various initials, so we want to clean this to capitalised, surname, single initial.
DB1 <- DB1 %>% 
  mutate(author_list = clean_author_list(author))

# and we also want to split the references similarly
# here we only do the most basic of cleaning (toupper) 
# because there are many different formats
# consider cleaning through revtools matches??
DB1 <- DB1 %>% 
  mutate(cited_references_list = map(cited_references, ~str_split(.x, pattern = "; ")[[1]])) %>% 
  mutate(cited_references_list = map(cited_references_list, toupper))

quick_check(DB1)  %>% print(n=nrow(.))

# finally, fill in age from year or early access date
DB1 <- DB1 %>% 
  mutate(age = as.integer(format(Sys.time(), "%Y")) - as.integer(ifelse(is.na(year), str_extract(string = early_access_date, pattern = "[0-9]{4}"), year)))

