# flashLIT functions natural clusters

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com


# ==================================================================


#' flashLIT_dtm
#' 
#' calculates a tidy document-term table. Inputs can be textcols (strings of text seperated by spaces), sepcols (strings of terms, seperated by a sepchr, e.g. "; ") or listcols (each row containing a list of entries). This function will process all into a tidy document-term-matrix.
#' 
#' textcols are, for example, title and abstract. These are united, and then split into single words, bigrams and trigrams. Words with fewer than smlwordlimit characters are removed before processing further.
#' 
#' sepcols are, for example, keywords and keywords_plus. These are united, then split into tokens using the sepchr. 
#' 
#' listcols are, for example, already cleaned author lists. These are unnested into tokens.
#' 
#' Then, textcols, sepcols, and splitcols are bound into a single dataframe, and frequencies are calculated using tidytext::bind_tf_idf. Sparse terms (i.e. terms found in a < minDocFreq proportion of documents) are removed, as are those that are too frequently found (i.e. terms found in a > maxDocFreq proportion of documents).
#' 
#' NOTE1: columns are not further cleaned before processing. all cleaning should be done prior to using this function.
#' NOTE2: no standardisation or weighting is applied. this should be done after function is applied.
#' NOTE3: it is recommended to seperate word-based columns from author-based columns from bibliography-based columns. these do not make sense to combine at this stage.
#'
#' @param df the tibble dataframe containing PID and desired columns
#' @param textcols specify the text columns as var(col1, col2) else NULL
#' @param sepcols specify the sep columns as var(col1, col2) else NULL
#' @param listcols specify the split columns as var(col1, col2) else NULL
#' @param smlwordlimit 
#' @param minDocFreq 
#' @param maxDocFreq  
#'
#' @return
#' @export
#'
#' @examples
flashLIT_dtm <- function(df, 
                        textcols = vars(title, abstract), 
                        sepcols = vars(keywords, keywords_plus),
                        listcols = NULL,
                        sepchr = "; ",
                        smlwordlimit = 3,
                        minDocFreq = 0.01,
                        maxDocFreq = 0.85
                        ){
  # initiate empty sets
  dtx <- tibble(PID = character(), token = character())
  dty <- tibble(PID = character(), token = character())
  dtz <- tibble(PID = character(), token = character())
  
  # merge text columns, and tokenize, to words, and ngrams, then filter out small words
  # text columns are sets of words seperated by a space
  if(!is.null(textcols)){
    x <- df %>% 
      select(PID, !!!textcols) %>% 
      unite("token", !!!textcols, sep = " ")
    dtx <- bind_rows(
      x %>% unnest_tokens(word, token) %>% rename(token = word),
      x %>% unnest_tokens(ngram, token, token = "ngrams", n = 2) %>% rename(token = ngram),
      x %>% unnest_tokens(ngram, token, token = "ngrams", n = 3) %>% rename(token = ngram)
    ) %>% 
    filter(!(str_length(token) < smlwordlimit))   
  }
  
  # merge sep columns, and tokenize
  # sep columns are sets of tokens seperated by a character
  if(!is.null(sepcols)){
    y <- df %>% 
      select(PID, !!!sepcols) %>% 
      unite("token", !!!sepcols, sep = sepchr)
    dty <-  y %>% unnest_tokens(word, token, token = stringr::str_split, pattern = sepchr) %>% rename(token = word)
  }  
  
  # merge split columns, and tokenize
  # split columns are already lists
  # we need to unnest these 
  if(!is.null(listcols)){
    dty <- df %>% 
      select(PID, !!!listcols) %>% 
      pivot_longer(cols = c(!!!listcols), values_to = "token") %>% 
      select(-name) %>% 
      unnest(token) 
  }  
   
  # combine all 
    dt <- bind_rows(dtx, dty, dtz) 

   # calculate frequencies  (tf = n/totalinPID, idf = ln(nDocs/nDocsWithTerm), tf_idf = tf*idf)
    dt <- dt %>% 
      count(PID, token, sort = TRUE) %>% 
      bind_tf_idf(token, PID, n)
    
   # remove the sparse terms
    dt <- dt %>% 
      filter(idf < log(nrow(df)/((minDocFreq)*nrow(df))) ) %>% 
      filter(idf > log(nrow(df)/((maxDocFreq)*nrow(df))) )
  return(dt)
}

#' flashLIT_distance
#' 
#' The document-term tidy_dtm is cast into a sparse matrix, and into a distance metric using ecodist (and the method specified by the distMethod). Euclidean distances are recommended for relatively dense matrices, whereas bray-curtis may be better for counts with many zeros. see other options in ?ecodist::distance
#'
#' @param dt document-term tidy dtm 
#' @param distMethod character indicating ecodist::distance method
#'
#' @return a distance object (includes all PIDs, no diagonals, unidirectional) converted into tidy format.
#' @export
#'
#' @examples
flashLIT_distance <- function(dt, doclist, distMethod = "euclidean"){  
  
  # cast to matrix, complete with missing PIDs
  tokens <- dt$token %>% unique() %>% sort()
  dm <- Matrix::sparseMatrix(
    i = match(dt$PID, doclist), 
    j = match(dt$token, tokens), 
    x = dt$n,
    dims = c(length(doclist), length(tokens)), 
    dimnames = list(doclist, tokens)
  )

  # calculate distance
  if(distMethod == "euclidean"){
    dst <- dist(dm, method = distMethod)
  } else if(distMethod == "bray-curtis"){
    dst <- ecodist::bcdist(dm)
  } else {
    dst <- dist(dm, method = distMethod)
  }
  
  # tidy and return  (includes all PIDs, no diagonals, unidirectional)
  tdst <- tidy(dst) %>% 
    rename(PID.x = item1, PID.y = item2, text_distance = distance)
  return(tdst)
}

# ==================================================================

#' scimeetr_coupling
#' 
#' a modified scimeetr process for 'coupling' i.e. creating a graph object by quantifying the similarity between document pairs. this creates a graph weighted by the sum similarity across the different components for the pairs of documents. Scimeetr uses https://jangorecki.gitlab.io/-/data.table/-/jobs/640724/artifacts/public/html/data.table.html whereas here we opt for the more transparent (but slower) dplyr process.
#' 
#' scores are first scaled to 0-1 (with NA replaced by 0) before being weighted, and summed. 
#' the graph is an (undirect, weighted) igraph object
#'
#' @param df a cleaned and processed dataframe
#' @param coupling_by text string, any combination of "abstract, title, keywords, keywordsplus, bibliography, authors, cr_journal"
#' @param w.tic weight for title component
#' @param w.kw weight for keyword (author supplied) component
#' @param w.kwp weight for keyword-plus component
#' @param w.abc weight for abstract component
#' @param w.auc weight for author matches component
#' @param w.joc weight for journal matches component
#' @param w.bic weight for bibliography matches component
#' @param ... placeholder for other subcalls, e.g. to include ngrams within the abstract and title overlaps
#'
#' @return an (undirect, weighted) igraph graph object
#' @export
#'
#' @examples
scimeetr_coupling <- function(df, 
                              coupling_by = "abstract, title, keywords, keywordsplus, bibliography, authors, cr_journal", 
                              w.tic = 1, w.kw = 1, w.kwp = 1, w.abc = 1, w.auc = 1, w.joc = 1, w.bic = 1, 
                              ...){
  
  # logical of which vars to couple by
  cby_title <- grepl(x = coupling_by, pattern = "title")
  cby_keywords <- grepl(x = coupling_by, pattern = "keywords")
  cby_keywordsplus <- grepl(x = coupling_by, pattern = "keywordsplus")
  cby_abstract <- grepl(x = coupling_by, pattern = "abstract")
  cby_authors <- grepl(x = coupling_by, pattern = "authors")
  cby_cr_journal <- grepl(coupling_by, pattern = "cr_journal")  
  cby_bibliography <- grepl(x = coupling_by, pattern = "bibliography")
  
  # initiate blanks of all of these
  couple_df_tic <- tibble(PID.x=character(), PID.y=character(), w_ij=numeric())
  couple_df_kw <- tibble(PID.x=character(), PID.y=character(), w_ij=numeric())
  couple_df_kwp <- tibble(PID.x=character(), PID.y=character(), w_ij=numeric())
  couple_df_abc <- tibble(PID.x=character(), PID.y=character(), w_ij=numeric())
  couple_df_auc <- tibble(PID.x=character(), PID.y=character(), w_ij=numeric())
  couple_df_joc <- tibble(PID.x=character(), PID.y=character(), w_ij=numeric())
  couple_df_bic <- tibble(PID.x=character(), PID.y=character(), w_ij=numeric())
  
  # add number of references column
  df <- df %>% 
    mutate(NR = map_dbl(cited_references_list, length))
  
  # then extract if required
  
  # Title coupling ---------------------------- 
  # this uses the pre-cleaned title. 
  # here, we use words from the titles to determine overlap - so draw these out using a dtm
  # still need to include ngrams in this process
  if(cby_title) {
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
      left_join(df[,c('PID', 'title')], by = c('PID.y' = 'PID'))%>%   
      mutate(NR.x = str_count(title.x, " "),
             NR.y = str_count(title.y, " "),
             w_ij = count/sqrt(NR.x * NR.y)) %>%
      select(PID.x, PID.y, w_ij)
    couple_df$w_ij[couple_df$w_ij == Inf] <- 0
    couple_df_tic <- filter(couple_df, w_ij != 0) %>% 
      ungroup()
  }
  
  # keyword coupling -------------------  
  # here each author supplied keyword or keyword phrase is an element. 
  # author supplied keywords are words authors have determined are relevant to their work.
  # these should be pre-cleaned for best results
  if(cby_keywords){
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
    couple_df_kw <- filter(couple_df, w_ij != 0) %>% 
      ungroup()
  }
  
  # keyword_plus coupling ------------------- 
  # keywords_plus are a WoS thing, aimed at extracting information from the entire cited references
  # and therefore useful for overlap in papers.
  # again, this should be a pre-cleaned list.
  if(cby_keywordsplus){
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
    couple_df_kwp <- filter(couple_df, w_ij != 0) %>% 
      ungroup()
  } 
  
  # Abstract coupling ---------------------------- 
  # here, we use words from the abstracts to determine overlap - so draw these out using a dtm
  # ideally using pre-cleaned abstracts
  # still need to include ngrams in this process
  if(cby_abstract){
    documents <- tm::Corpus(tm::VectorSource(df$abstract))
    myTdm <- tm::DocumentTermMatrix(documents)
    myTdm2 <- tm::removeSparseTerms(myTdm, sparse = 0.99)
    dtm2list <- apply(myTdm2, 1, function(x) {
      paste(rep(names(x), x), collapse=" ")
    })
    ab_list <- strsplit(dtm2list, "[ ]")
    names(ab_list) <- df$PID
    ab_df <- data.frame('PID'= rep(names(ab_list), sapply(ab_list, length)),
                        'AB' =  unlist(ab_list),
                        stringsAsFactors=F)
    couple_df <- inner_join(ab_df, ab_df, by = 'AB') %>%
      filter(PID.x > PID.y) %>%
      group_by(PID.x, PID.y) %>%
      summarise(count = n()) %>%
      left_join(df[,c('PID', 'abstract')], by = c('PID.x' = 'PID'))%>%
      left_join(df[,c('PID', 'abstract')], by = c('PID.y' = 'PID'))%>%   
      mutate(NR.x = str_count(abstract.x, " "),
             NR.y = str_count(abstract.y, " "),
             w_ij = count/sqrt(NR.x * NR.y)) %>%
      select(PID.x, PID.y, w_ij)
    couple_df$w_ij[couple_df$w_ij == Inf] <- 0
    couple_df_abc <- filter(couple_df, w_ij != 0) %>% 
      ungroup()
  } 
  
  # author coupling --------------------------------------------------- 
  # similar to keywords, here we use a pre-cleaned author list
  if(cby_authors){
    au_list <- df$author_list
    names(au_list) <- df$PID
    au_df <- data.frame('PID'= rep(names(au_list), sapply(au_list, length)),
                        'AU' =  unlist(au_list),
                        stringsAsFactors=F)
    au_length <- group_by(au_df, PID) %>%
      summarize(NAU = n())
    couple_df <- inner_join(au_df, au_df, by = 'AU') %>%
      filter(PID.x > PID.y) %>%
      group_by(PID.x, PID.y) %>%
      summarise(count = n()) %>%
      left_join(au_length, by = c('PID.x' = 'PID'))%>%
      left_join(au_length, by = c('PID.y' = 'PID'))%>%
      mutate(w_ij = count/sqrt(NAU.x * NAU.y)) %>%
      select(PID.x, PID.y, w_ij)
    couple_df$w_ij[couple_df$w_ij == Inf] <- 0
    couple_df_auc <- filter(couple_df, w_ij != 0) %>% 
      ungroup()
  }
  
  # Journal coupling ------------------------------------------------------ 
  # here the journals came from the cited references (not the journals of the papers themselves)
  # we havent yet pre-cleaned this, aside from pre-cleaning the cited references
  if(cby_cr_journal){
    cr_list <- df$cited_references_list
    names(cr_list) <- df$PID
    cr_df <- data.frame('PID' = rep(names(cr_list), sapply(cr_list, length)),
                        'CR' = unlist(cr_list),
                        stringsAsFactors=F)
    # extract the entry after the year and before the next comma
    cr_df <- cr_df %>% 
      mutate(CRJ = map_chr(CR, ~str_extract(.x, pattern = "[0-9]{4}[,][ ][^,]+") %>% str_remove(pattern = "[0-9]{4}[,][ ]"))) %>% 
      select(-CR)
    rm(cr_list)
    tmp <- cr_df %>% 
      group_by(PID, CRJ) %>%
      summarise(jo_freq = n()) %>%
      filter(!is.na(CRJ)) %>% 
      filter(!CRJ == "NA")
    couple_df <- inner_join(tmp, tmp, by = 'CRJ') %>%
      filter(PID.x > PID.y) %>%
      mutate(min_jo = min(jo_freq.x, jo_freq.y)) %>%
      select(PID.x, PID.y, min_jo) %>%
      group_by(PID.x, PID.y) %>%
      summarise(count = sum(min_jo)) %>%
      left_join(df[,c('PID', 'NR')], by = c('PID.x' = 'PID'))%>%
      left_join(df[,c('PID', 'NR')], by = c('PID.y' = 'PID'))%>%
      mutate(w_ij = count/sqrt(NR.x * NR.y)) %>%
      select(PID.x, PID.y, w_ij)
    rm(tmp)
    gc()  # It can be useful to call gc after a large object has been removed, as this may prompt R to return memory to the operating system.
    couple_df$w_ij[couple_df$w_ij == Inf] <- 0
    couple_df_joc <- filter(couple_df, w_ij != 0) %>% 
      ungroup()
  }
  
  # Bibliographic coupling -------------------- 
  # here each cr is an element. 
  # this should be a pre-cleaned list named cr_list
  if(cby_bibliography){
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
    couple_df_bic <- filter(couple_df, w_ij != 0) %>% 
      ungroup()
  } 
  
  
  # Join all the coupled data groups ---------------------------------------------
  # interestingly this misses a suffix on the second join, no matter what it is 
  couple_df <- couple_df_tic %>% 
    full_join(couple_df_kw, by = c('PID.x', 'PID.y'), suffix = c(".tic", ".kw")) %>%
    full_join(couple_df_kwp, by = c('PID.x', 'PID.y'), suffix = c('', '.kwp')) %>%
    full_join(couple_df_abc, by = c('PID.x', 'PID.y'), suffix = c('', '.abc')) %>%
    full_join(couple_df_auc, by = c('PID.x', 'PID.y'), suffix = c('', '.auc')) %>%
    full_join(couple_df_joc, by = c('PID.x', 'PID.y'), suffix = c('', '.joc')) %>%
    full_join(couple_df_bic, by = c('PID.x', 'PID.y'), suffix = c('', '.bic')) %>% 
    rename(w_ij.kwp = w_ij)
  
  # replace all the na with zero, and scale to 1
  couple_df <- couple_df %>% 
    mutate(across(everything(), replace_na, replace = 0)) %>% 
    mutate(across(w_ij.tic:w_ij.bic, scales::rescale)) 
  
  # weight 
  couple_df <- couple_df %>% 
    mutate(w_ij = 
             w_ij.tic * w.tic +
             w_ij.kw  * w.kw +
             w_ij.kwp * w.kwp +
             w_ij.abc * w.abc +
             w_ij.auc * w.auc +
             w_ij.joc * w.joc +
             w_ij.bic * w.bic 
    ) %>%
    select(PID.x, PID.y, w_ij)
  
  # and create the graph
  m <- sum(couple_df$w_ij)
  coup2 <- data.frame(PID = c(as.vector(couple_df$PID.x), as.vector(couple_df$PID.y)),
                      w_ij = rep(couple_df$w_ij,2))
  k_i <- coup2 %>%
    group_by(PID) %>%
    summarise(k_i = sum(w_ij))
  rm(coup2)
  nam <- as.character(k_i$PID)
  k_i <- k_i$k_i
  names(k_i) <- nam
  couple_df <- couple_df %>%
    mutate(asso_stre = (2* w_ij * m)/ (k_i[PID.x] * k_i[PID.y]))
  rm(k_i)
  names(couple_df) <- c("PID1",
                        "PID2",
                        "bc",
                        "weight"
  )
  couple_df <- ungroup(couple_df)
  tmp <- df$PID[which(!(df$PID %in% unique(c(couple_df$PID2, couple_df$PID1))))]
  if(length(tmp) >= 1){
    missing_df <- data.frame('PID1' = df$PID[which(!(df$PID %in% unique(c(couple_df$PID2, couple_df$PID1))))],
                             'PID2' = df$PID[1],
                             'bc' = 0,
                             'weight' = 0,
                             stringsAsFactors = F)
    couple_df <- rbind(couple_df, missing_df)
  }
  graph <- igraph::graph_from_data_frame(d=couple_df, directed= F)
  
  return(graph)
}

# ===============================================================

#' clusterize
#' 
#' graph clusterisation wrapper, directly copied from scimeetr. Two options are offered here, although there are many other possible cluster algorithms.
#'
#' @param graph an igrpah object
#' @param community_algorithm character string, either 'louvain' (default) or 'fast greedy'
#'
#' @return returns the igraph community object output
#' @export
#'
#' @examples
clusterize <- function(graph = graph, community_algorithm = 'louvain'){
  if(community_algorithm == 'louvain'){
    community <- igraph::cluster_louvain(graph)
  } else if(community_algorithm == 'fast greedy'){
    community <- igraph::cluster_fast_greedy(graph)
  }
  return(community)
}

# ====================================================================

#' iteratively_cluster
#' 
#' apply a clustering algorithm iteratively to the data, until all groups are smaller than the specified maxgroupsize
#'
#' @param graph the igraph object to cluster (made with scimeetr_coupling())
#' @param maxgroupsize if group size is larger than this number, the group will be further split
#' @param community_algorithm character string, either 'louvain' (default) or 'fast greedy' passed to the cluster algorithm
#'
#' @return a tibble with PID of all the documents, and columns with each iterative split level. 
#' @export
#'
#' @examples
iteratively_cluster <- function(graph,            
                                maxgroupsize = igraph::gorder(graph)/20,  
                                community_algorithm = 'louvain'           
){
  
  # initiate the first cluster and determine current groups
  clg <- clusterize(graph = graph) %>% igraph::membership()
  cl <- tibble(PID = names(clg),
               clg = as.character(clg))
  clbase <- cl %>% add_column(cl_0 = "0", .before = 2)
  
  # determine how many in the groups, and therefore whether to cluster them further (if still above maxgroupsize). 
  gtb <- table(cl$clg)
  to_split <- gtb[gtb > maxgroupsize] %>% names()
  no_split <- gtb[gtb <= maxgroupsize] %>% names()
  if(length(to_split) == 0){stop("All clusters already smaller than maxgroupsize.")}
  
  while(length(to_split) > 0){
    
    # map over these, subset the graph, calculate the clusters, and append to the cluster object
    clsplit <- map_dfr(to_split, function(x){
      vpids <- cl %>% filter(clg == x) %>% pull(PID)
      sg <- igraph::induced_subgraph(graph = graph, vids = vpids, impl = "auto")
      clg <- clusterize(graph = sg) %>% igraph::membership() 
      clsplit <- tibble(
        PID = names(clg),
        clg = paste(x, clg, sep = "_"))
      clsplit
    })
    
    # for the non-split groups, add these directly with a zero group
    clns <- cl %>% 
      filter(clg %in% no_split) %>% 
      mutate(clg = paste(clg, 0, sep = "_"))
    
    # join together, add to original, and count the number of groups to determine if the split needs to be continued
    clnew <- bind_rows(clsplit, clns) 
    
    clbase <- clbase %>% full_join(clnew, by = "PID")
    
    cl <- clnew
    gtb <- table(cl$clg)
    to_split <- gtb[gtb > maxgroupsize] %>% names()
    no_split <- gtb[gtb <= maxgroupsize] %>% names()
  }  
  
  # then rename all the clbase columns
  names(clbase)[-1] <- paste("cl", 0:(length(clbase)-2), sep = "_")
  return(clbase)
}
