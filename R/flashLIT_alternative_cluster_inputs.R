# flashLIT alternative cluster inputs

# ===================================

## cluster options:
# 1) TAK-raw, 2) TAK-clean 3)TAK-stem
# 4) TAKKP-stem, 5) Bib, 6) TAKKP-stem-Au and 6) TAKKP-stem-Bib.
# we will be comparing across 1:3 and 3:7 separately.

## distance method options: euclidean or bray-curtis

# ==================================

# these may need to be modified to retain all the PIDs possible.

# text-based only ---

TAK_raw <- flashLIT_dtm(df = DB1, 
                        textcols = vars(title_orig, abstract_orig), 
                        sepcols = vars(keywords_orig),
                        listcols = NULL,
                        sepchr = "; ",
                        smlwordlimit = 3,
                        minDocFreq = 0.01,
                        maxDocFreq = 0.85)

TAK_clean <- flashLIT_dtm(df = DB1, 
                          textcols = vars(title_clean, abstract_clean), 
                          sepcols = vars(keywords_clean),
                          listcols = NULL,
                          sepchr = "; ",
                          smlwordlimit = 3,
                          minDocFreq = 0.01,
                          maxDocFreq = 0.85)

TAK_stem <- flashLIT_dtm(df = DB1, 
                         textcols = vars(title_stem, abstract_stem), 
                         sepcols = vars(keywords_stem),
                         listcols = NULL,
                         sepchr = "; ",
                         smlwordlimit = 3,
                         minDocFreq = 0.01,
                         maxDocFreq = 0.85)

TAKKP_stem <- flashLIT_dtm(df = DB1, 
                           textcols = vars(title_stem, abstract_stem), 
                           sepcols = vars(keywords_stem, keywordsPlus_stem),
                           listcols = NULL,
                           sepchr = "; ",
                           smlwordlimit = 3,
                           minDocFreq = 0.01,
                           maxDocFreq = 0.85)

# text and other information based ---
# note change in smlwordlimit and min-max docfreq
Au <- flashLIT_dtm(df = DB1, 
                   textcols = NULL, 
                   sepcols = NULL,
                   listcols = vars(author_list),
                   sepchr = "; ",
                   smlwordlimit = 0,
                   minDocFreq = 1/nrow(DB1),
                   maxDocFreq = 0.99)

Bib <- flashLIT_dtm(df = DB1, 
                   textcols = NULL, 
                   sepcols = NULL,
                   listcols = vars(cited_references_list),
                   sepchr = "; ",
                   smlwordlimit = 0,
                   minDocFreq = 1/nrow(DB1),
                   maxDocFreq = 0.99)

## combine to distance matrices ---
TAK_raw_E <- TAK_raw %>% flashLIT_distance(doclist = DB1$PID, distMethod = "euclidean")
TAK_raw_BC <- TAK_raw %>% flashLIT_distance(doclist = DB1$PID, distMethod = "bray-curtis")
TAK_clean_E <- TAK_clean %>% flashLIT_distance(doclist = DB1$PID, distMethod = "euclidean")
TAK_clean_BC <- TAK_clean %>% flashLIT_distance(doclist = DB1$PID, distMethod = "bray-curtis")
TAK_stem_E <- TAK_stem %>% flashLIT_distance(doclist = DB1$PID, distMethod = "euclidean")
TAK_stem_BC <- TAK_stem %>% flashLIT_distance(doclist = DB1$PID, distMethod = "bray-curtis")

TAKKP_stem_E <- TAKKP_stem %>% flashLIT_distance(doclist = DB1$PID, distMethod = "euclidean")
TAKKP_stem_BC <- TAKKP_stem %>% flashLIT_distance(doclist = DB1$PID, distMethod = "bray-curtis")

Au_E <- Au %>% flashLIT_distance(doclist = DB1$PID, distMethod = "euclidean") %>% rename(au_distance = text_distance)
Au_BC <- Au %>% flashLIT_distance(doclist = DB1$PID, distMethod = "bray-curtis") %>% rename(au_distance = text_distance)

Bib_E <- Bib %>% flashLIT_distance(doclist = DB1$PID, distMethod = "euclidean") %>% rename(bib_distance = text_distance)
Bib_BC <- Bib %>% flashLIT_distance(doclist = DB1$PID, distMethod = "bray-curtis") %>% rename(bib_distance = text_distance)

# 6) TAKKP-stem-Au and 6) TAKKP-stem-Bib.
TAKKP_stem_Au_E <- full_join(TAKKP_stem_E, Au_E)
TAKKP_stem_Au_BC <- full_join(TAKKP_stem_BC, Au_BC)
TAKKP_stem_Bib_E <- full_join(TAKKP_stem_E, Bib_E)
TAKKP_stem_Bib_BC <- full_join(TAKKP_stem_BC, Bib_BC)