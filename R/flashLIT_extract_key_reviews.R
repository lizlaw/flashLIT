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
  filter(grepl(pattern ='review|metaanalysis|meta analysis', x = abstract)) %>% 
  add_ref()
# 69 papers

compreviews <- allreviews %>% 
  filter(grepl(pattern ='competition', x = abstract))
# 13 papers

floralcompreviews <- compreviews %>% 
  filter(grepl(pattern ='floral|flower|food', x = abstract)) %>% 
  pull(ref)
  #pull(abstract_orig)
# 10 papers

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
"CHANDLER D, et al., 2019, Are there risks to wild European bumble bees from using commercial stocks of domesticated Bombus terrestris for crop pollination?, J APICULT RES"
"STOUT J, et al., 2009, Ecological impacts of invasive alien species on bees, APIDOLOGIE"  

# for each of these papers, we can download both their bibliographies (or the lists within a paper), and the papers that cite them. 
# for those with specific and seperate categories (e.g. presented in tables in the papers), we can futher identify the papers allocated to each category
# these can be used to test the categorisation.

# I did this from WoS, but matt says:
# We have updated citationchaser (citationchaser (shinyapps.io)) so you dont need an API from the Shinyserver (still need one in the R package). That might help with back/forward citation.

list.files("data/raw_WoS_20210415/FromKeyReviews")
list.files("data/KeyReferences") # includes ones we can further categorise into lists, and the potential lists that we can use from e.g. tables

# Importing the key review papers -----------------------------

KR0_cites <- import_wos_files_DFonly("./data/raw_WoS_20210415/FromKeyReviews/cites/")
KR0_refs <- import_wos_ref_files("./data/raw_WoS_20210415/FromKeyReviews/refs/")

# select and clean
KR1_cites <- KR0_cites %>% 
  convert_names(from = "wos", to = "bib") %>% 
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
  tibble::rownames_to_column(var = 'KRC_ID') %>% 
  dplyr::mutate(KRC_ID = paste0('KRC_ID', stringr::str_pad(KRC_ID, 5, 'left', 0))) %>% 
  dplyr::mutate(across(everything(), na_if, y="")) %>% 
  mutate(across(c(title, abstract), clean_text)) %>% 
  mutate(across(c(title, abstract), stemR))  %>% 
  mutate(across(c(title, abstract), stemCompletR, .dict = stemDict))

KR1_refs <- KR0_refs %>% 
  convert_names(from = "wos", to = "bib") %>% 
  tibble::tibble() %>% 
  dplyr::select(doi, 
                author, 
                title, 
                abstract, 
                n_cited_allwos,  
                year, 
                early_access_date,
                #journal, journal_iso,
                journal_iso = source, 
                volume, 
                pages, 
                article_number,
                folder,
                file
  ) %>% 
  dplyr::mutate(
    title_orig = title,
    abstract_orig = abstract
  ) %>% 
  tibble::rownames_to_column(var = 'KRR_ID') %>% 
  dplyr::mutate(KRR_ID = paste0('KRR_ID', stringr::str_pad(KRR_ID, 5, 'left', 0))) %>% 
  dplyr::mutate(across(everything(), na_if, y="")) %>% 
  mutate(across(c(title, abstract), clean_text)) %>% 
  mutate(across(c(title, abstract), stemR))  %>% 
  mutate(across(c(title, abstract), stemCompletR, .dict = stemDict))

# ==================================================================
# clean the KR1_cites and KR1_refs to unique entries before proceeding

# KR1_cites ---

dup_title <- KR1_cites %>% filter(duplicated(title)) %>% pull(title)
dup_doi <- KR1_cites %>% filter(!is.na(doi)) %>% filter(duplicated(doi)) %>% pull(doi)

# check manually - OK happy all duplicates are matches based on the title/author
KR1_cites %>% filter(doi %in% dup_doi | title %in% dup_title) %>% arrange(title) %>% select(title_orig, year, journal_iso, volume, pages, author) %>% print(n=nrow(.))

# check doi matches are the same as the title matches
KR1_cites %>% filter(!is.na(doi)) %>% filter(duplicated(doi)) %>% pull(title) %>% .[]== dup_title

# check others for fuzzy matches
fuzzymatch_potential_duplicates(KR1_cites$title_orig)

# happy no other fuzzy matches are duplicates

# Now, we pivot wider the title and file, then join back onto the de-duplicated data

kr1_files <- KR1_cites %>% 
  select(title, file) %>% 
  mutate(val = 1) %>% 
  pivot_wider(names_from = file, values_from = val, values_fill = 0)

kr1_dedup <- KR1_cites %>% 
  filter(!duplicated(title)) %>% 
  select(-file)

KR2_cites <- full_join(kr1_dedup, kr1_files)

# KR2_refs --- 

dup_title <- KR1_refs %>% filter(duplicated(title)) %>% pull(title)
dup_doi <- KR1_refs %>% filter(!is.na(doi)) %>% filter(duplicated(doi)) %>% pull(doi)

KR1_refs %>% filter(doi %in% dup_doi | title %in% dup_title) %>% arrange(title) %>% select(title_orig) %>% summarise(n=n()) %>% filter(n<2)
KR1_refs %>% filter(doi %in% dup_doi | title %in% dup_title) %>% arrange(title) %>% select(author) %>% summarise(n=n()) %>% filter(n<2)

# check manually - OK happy all duplicates are matches based on the title/author
KR1_refs %>% filter(doi %in% dup_doi | title %in% dup_title) %>% arrange(title) %>% select(title_orig, year, journal_iso, volume, pages, author) %>% print(n=nrow(.))

# doi matches are different to title matches this time, but all doi are within title matches
KR1_refs %>% filter(!is.na(doi)) %>% filter(doi %in% dup_doi) %>% pull(title) %in% dup_title

# fuzzy matches
fuzzymatch_potential_duplicates(KR1_refs$title_orig) %>% print(n=nrow(.))
fuzzymatch_potential_duplicates(KR1_refs$title) %>% print(n=nrow(.))

# manual check of potential matches
checkmatch_potential_duplicates("Honey bee colony losses", list(KR1_refs), title_orig, year, journal_iso, author)
checkmatch_potential_duplicates("OCCURRENCE OF CHALKBROOD CAUSED BY ASCOSPHAERA-AGGREGATA", list(KR1_refs), title_orig, year, journal_iso, author)

# the first is duplicated, but this is captured in the true duplicates, and the fuzzy match is not a duplicate.

# Now, we pivot wider the title and file, then join back onto the de-duplicated data

kr1_files <- KR1_refs %>% 
  select(title, file) %>% 
  mutate(val = 1) %>% 
  pivot_wider(names_from = file, values_from = val, values_fill = 0)

kr1_dedup <- KR1_refs %>% 
  filter(!duplicated(title)) %>% 
  select(-file)

KR2_refs <- full_join(kr1_dedup, kr1_files)

# ====================================================

# Matches with our database

# first check the fuzzy matches with the KR2cites:
fuzzymatch_potential_duplicates(V1 = DB1$title_orig, V2 = KR2_cites$title_orig)
# no fuzzy matches

# match exact matches using PID
KR2_cites <- KR2_cites %>% 
  left_join(DB1 %>% select(PID, title_orig))

# then check the fuzzy matches with the KR2refs:
fuzzymatch_potential_duplicates(V1 = DB1$title_orig, V2 = KR2_refs$title_orig)

# happy no fuzzy matches to contend with

# match the exact matches using PID
KR2_refs <- KR2_refs %>% 
  left_join(DB1 %>% select(PID, title_orig)) 

# examine output
KR2_refs %>% filter(!is.na(PID)) %>% summarise(across(Arena_refs.txt:Wood_refs.txt, sum)) %>% t()
KR2_cites %>% filter(!is.na(PID)) %>% summarise(across(cites_Arena.txt:cites_Wood.txt, sum)) %>% t()
