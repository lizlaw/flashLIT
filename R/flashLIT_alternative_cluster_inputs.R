# flashLIT alternative cluster inputs

# ===================================

## cluster options:
# 1) TAK-raw, 2) TAK-clean 3)TAK-stem
# 4) TAKKP-stem, 5) Bib, 6) TAKKP-stem-Au and 6) TAKKP-stem-Bib.
# the latter two can be normalised and then joined to a single distance with mean or euclidean. 
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

Au_clean_E <- Au %>% flashLIT_distance(doclist = DB1$PID, distMethod = "euclidean") %>% rename(au_distance = text_distance)
Au_clean_BC <- Au %>% flashLIT_distance(doclist = DB1$PID, distMethod = "bray-curtis") %>% rename(au_distance = text_distance)

Bib_clean_E <- Bib %>% flashLIT_distance(doclist = DB1$PID, distMethod = "euclidean") %>% rename(bib_distance = text_distance)
Bib_clean_BC <- Bib %>% flashLIT_distance(doclist = DB1$PID, distMethod = "bray-curtis") %>% rename(bib_distance = text_distance)

# 6) TAKKP-stem-Au and 6) TAKKP-stem-Bib. ---
# First join the information 
TAKKPAu_stem_E <- full_join(TAKKP_stem_E, Au_clean_E)
TAKKPAu_stem_BC <- full_join(TAKKP_stem_BC, Au_clean_BC)
TAKKPBib_stem_E <- full_join(TAKKP_stem_E, Bib_clean_E)
TAKKPBib_stem_BC <- full_join(TAKKP_stem_BC, Bib_clean_BC)

# then normalise and average, or euclidean distance to one metric
TAKKPAu_stem_Em <- dist_rescaler(TAKKPAu_stem_E, text_distance, au_distance, mean)
TAKKPAu_stem_BCm <- dist_rescaler(TAKKPAu_stem_BC, text_distance, au_distance, mean)
TAKKPBib_stem_Em <- dist_rescaler(TAKKPBib_stem_E, text_distance, bib_distance, mean)
TAKKPBib_stem_BCm <- dist_rescaler(TAKKPBib_stem_BC, text_distance, bib_distance, mean)

TAKKPAu_stem_Ee <- dist_rescaler(TAKKPAu_stem_E, text_distance, au_distance, euclidean)
TAKKPAu_stem_BCe <- dist_rescaler(TAKKPAu_stem_BC, text_distance, au_distance, euclidean)
TAKKPBib_stem_Ee <- dist_rescaler(TAKKPBib_stem_E, text_distance, bib_distance, euclidean)
TAKKPBib_stem_BCe <- dist_rescaler(TAKKPBib_stem_BC, text_distance, bib_distance, euclidean)