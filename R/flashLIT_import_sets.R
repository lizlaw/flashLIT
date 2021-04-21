# flashLIT import sets

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

# 1. download all sets into seperate dfs
# 2. strip the folder/file names from a copy, and join together
# 3. deduplicate
# 4. re-join set information
# 5. fuzzy-deduplicate

# download all sets ---

ddd <- list.files("data/raw_WoS_20210415", full.names = TRUE)[-1]
DL0 <- list()

for (di in seq(ddd)){
  DL0[[di]] <- import_wos_files_DFonly(paste0(ddd[di], "/"))
}

DL0 <- DL0 %>% 
  bind_rows() %>% 
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
  ) 

# remove duplicates ---

# remove exact duplicates of all columns
DL1 <- DL0 %>% 
  select(-folder, -file) %>% 
  filter(!duplicated(.))
  
# check fuzzy title duplicates
m <- fuzzymatch_potential_duplicates(DL1$title) 
m %>% print(n=nrow(.))

checkmatch_potential_duplicates("Variation in gut microbial communities and its association with pathogen", list(DL1), title, year, author, journal_iso) 
checkmatch_potential_duplicates("Wild Pollinators Enhance Fruit Set of Crops Regardless of Honey Bee", list(DL1), title, year, author, journal_iso) 
checkmatch_potential_duplicates("Evidence for multiple introductions of an invasive wild bee species", list(DL1), title, year, author, journal_iso) 

# match lookup
match_lookup <- tibble(
  change = c(
    "Variation in gut microbial communities and its association with pathogen infection in wild bumble bees (Bombus) (vol 8, pg 2369, 2014)",
    "Wild Pollinators Enhance Fruit Set of Crops Regardless of Honey Bee Abundance(vol 344, 2014)",
    "Evidence for multiple introductions of an invasive wild bee species currently under rapid range expansion in Europe (vol 21, 17, 2021)"
  ),
  keep = c(
    "Variation in gut microbial communities and its association with pathogen infection in wild bumble bees (Bombus)",
    "Wild Pollinators Enhance Fruit Set of Crops Regardless of Honey Bee Abundance",
    "Evidence for multiple introductions of an invasive wild bee species currently under rapid range expansion in Europe"
  )
)

DL1 <- DL1 %>% 
  mutate(alt_title = map_chr(title, ~ifelse(.x %in% match_lookup$keep,
                                        match_lookup$change[which(match_lookup$keep == .x)],
                                        .x))) %>% 
  mutate(title = map_chr(title, ~ifelse(.x %in% match_lookup$change,
                                          match_lookup$keep[which(match_lookup$change == .x)],
                                          .x)))

# check exact title duplicates (that the other details also match)
mt <- DL1 %>% 
  filter(duplicated(DL1$title)) %>% 
  pull(title)

m <- DL1 %>% filter(title %in% mt) %>% arrange(title) %>% group_by(title, author, year, journal_iso, volume, pages) %>% summarise(n=n()) 

m %>% filter(n<2) 

# correct version for "Peer review of the pesticide risk assessment for the active substance sulfoxaflor in light of confirmatory data submitted" is the 2020 volume 18 version
DL1 <- DL1 %>% filter(!(title == "Peer review of the pesticide risk assessment for the active substance sulfoxaflor in light of confirmatory data submitted" & volume == "17"))


# check for duplicated titles within this set
m %>% filter(duplicated(title))
m %>% filter(grepl("Variation in gut microbial communities and its association with pathogen infection in", title))

## corrected version is the pages 2550 version (corregerium)
DL1 <- DL1 %>% filter(!(grepl("Variation in gut microbial communities and its association with pathogen infection in", title) & pages == "2369"))

# now is ok to remove all duplicated title versions
DL1 <- DL1 %>% filter(!duplicated(title))

# match the groups on based on (corrected) title, note if in the folder ---------------------------

DL0w <- DL0 %>% 
  select(title, folder) %>% 
  mutate(title = map_chr(title, ~ifelse(.x %in% match_lookup$change,
                                        match_lookup$keep[which(match_lookup$change == .x)],
                                        .x))) %>% 
  filter(!duplicated(.)) %>% 
  mutate(folder = str_remove(folder, "data/raw_WoS_20210415/") %>% str_remove("/")) %>% 
  mutate(set = 1) %>% 
  pivot_wider(names_from = folder, values_from = set, values_fill = 0)

# DL0w %>% filter_at(vars(MN_n1_TAK:MNI_n3_TOPIC), any_vars(. > 1) ) %>% select(title, year, volume, pages, MN_n1_TAK:MNI_n3_TOPIC)

# join
DL2 <- DL1 %>% 
  right_join(DL0w)

# and clean, adding a PID
DL2 <- DL2 %>% 
  tibble::rownames_to_column(var = 'PID') %>% 
  dplyr::mutate(PID = paste0('PID', stringr::str_pad(PID, 5, 'left', 0))) %>% 
  dplyr::mutate(across(everything(), na_if, y=""))


  
# ===================== 
  # older import of one set only
  
  # # import data. 
  # # Data can come as a bibliography in many formats, and can be imported via revtools or scimeetr.
  # # Revtools imports from either .bib or .txt (.bib I think works better). data are given revtools names:
  # # https://github.com/mjwestgate/revtools/blob/master/R/tag_lookup.R
  # # Alternatively, scimeetr imports to a specialised object, from wos or scopus ris text files ONLY, automatically removing duplicates. 
  # 
  # # revtools bib
  # # dd <- "./data/raw_WoS_20201105/as_bib/"
  # # rev.df <- revtools::read_bibliography(list.files(dd, full.names = TRUE))
  # 
  # # revtools txt - currently fails 
  # # dd <- "./data/raw_WoS_20201105/as_txt/"
  # # revtools::read_bibliography(list.files(dd, full.names = TRUE))
  # 
  # # scimeetr import 
  # # dd <- "./data/raw_WoS_20201105/as_txt/"
  # # sci.df <- scimeetr::import_wos_files(dd)$com1$dfsci
  # 
  # # customised scimeetr import (straight to df with file and folder information)
  # dd <- "./data/raw_WoS_20210415/MNI_n1_TAK/"
  # sci.df <- import_wos_files_DFonly(dd)
  # 
  # # convert the scimeetr import to 'revtools' format
  # scib.df <- sci.df %>% convert_names(from = "wos", to = "bib")
  # 
  # # check and remove possible duplicates
  # scib.df %>% filter(duplicated(title)) %>% select(title, year, author, source)
  # fuzzymatch_potential_duplicates(scib.df$title) 
  # checkmatch_potential_duplicates("Native bee diversity and abundance in an urban green space in Bengaluru, India", list(scib.df), title, year, author, source) 
  # checkmatch_potential_duplicates("The native bee fauna of Carlinville, Illinois, revisited after 75 years: a case for persistence", list(scib.df), title, year, author, source) 
  # checkmatch_potential_duplicates("Wild Pollinators Enhance Fruit Set of Crops Regardless of Honey Bee Abundance", list(scib.df), title, year, author, source) 
  # 
  # scib.df <- scib.df %>% filter(!duplicated(title)) %>% filter(!title == "Wild Pollinators Enhance Fruit Set of Crops Regardless of Honey Bee Abundance(vol 344, 2014)")
  # 
  # # choose columns that we will use (note choice may be dictated by import type) and create PID
  # DB0 <- scib.df %>% 
  #   tibble::tibble() %>% 
  #   dplyr::select(doi, 
  #                 author, 
  #                 title, 
  #                 keywords = author_keywords, 
  #                 keywords_plus, 
  #                 abstract, 
  #                 cited_references,  
  #                 language, 
  #                 document_type, 
  #                 n_cited_allwos,  
  #                 year, 
  #                 early_access_date,
  #                 #journal, journal_iso,
  #                 journal_iso = source_abbreviation_29char, 
  #                 volume, 
  #                 pages, 
  #                 article_number,
  #                 date_generated,
  #                 folder,
  #                 file
  #   ) %>% 
  #   dplyr::mutate(
  #     title_orig = title,
  #     abstract_orig = abstract,
  #     keywords_orig = keywords
  #   ) %>% 
  #   tibble::rownames_to_column(var = 'PID') %>% 
  #   dplyr::mutate(PID = paste0('PID', stringr::str_pad(PID, 5, 'left', 0))) %>% 
  #   dplyr::mutate(across(everything(), na_if, y=""))