# flashLIT group summaries

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# !! documentation needs to be completed for these

# ==================================================================

#' age_based_residuals
#' 
#' using mgcv::gam, this is a poisson model with a default k=10 relating e.g the number of citations against age. flexible so other models can be specified with the same formatting. Output is the residuals (with included NA). This default k is not always appropriate, and it will fail with an error, so recommend adjusting the k or wrapping with tryCatch if this might be problematic (e.g. with smaller group sizes).
#'
#' @param df a bibliographic dataframe
#' @param .ycol the y variable column for the age model
#' @param .xcol the x variable colomn, ie the age, for the age model
#' @param .k the knot level of the gam model
#'
#' @return the residuals from the age based model, with NA inserted from original Y vector
#' @export
#'
#' @examples
age_based_residuals <- function(df, .ycol, .xcol = "age", .k = 10){
  y <- df[.ycol] %>% pull()
  x <- df[.xcol] %>% pull()
  z <- mgcv::gam(y ~ s(x, k=.k), family = "poisson")$residuals
  y[!is.na(y)] <- z
  return(y)
}

# ==================================================================

# a function to give a contensed version of the reference for use in the summaries
# Author1, Year, Title, Journal

#' add_ref
#' 
#' a function to give a contensed version of the reference for use in the summaries (Author1, Year, Title, Journal, and optionally DOI)
#'
#' @param df the dataframe to add a ref column to
#' @param doi logical, whether to include DOI
#'
#' @return the original dataframe with appended ref column
#' @export
#'
#' @examples
add_ref <- function(df, .doi = FALSE){
  df %>% 
    mutate(ref = pmap_chr(list(df$author_list, df$year, df$title_orig, df$journal_iso, df$doi),
                          function(au, yr, ti, jo, di){
                            if(length(au)>1) au <- paste(au[[1]], "et al.", sep = ", ") 
                            x <- paste(au[[1]], yr, ti, jo, sep = ", ")
                            if(.doi) x <- paste(x, di, sep = ", ")
                            x
                          }))
}

# ==================================================================

# functions to add groupspecific frequencies to the data 
# most are modified from scimeetr
# input is the group dataset, ie. only the group.

# ------------------------------------------------------------------
# most cited within group (i.e. by the group)
add_group_most_cited <- function(df,...){  
  ncites <- df %>% 
    pull(citedPID_list) %>% 
    unlist() %>% 
    as_tibble() %>% 
    rename(PID = value) %>% 
    group_by(PID) %>% 
    summarise(cites_within_group=n())
  df %>% 
    left_join(ncites, by = "PID")
}

# ------------------------------------------------------------------
# cites highest number of group papers 
# need to filter out the within group pids from the list of cited pids and then count them
add_group_cites_most_group <- function(df,...){
  df %>% 
    mutate(cites_of_group = map_dbl(citedPID_list, function(x) sum(x %in% df$PID)))
}

# ------------------------------------------------------------------
# papers by most frequent authors within group - sum of author frequency within group
# NOTE this is NOT scaled by the number of authors on a paper. We tentatively think it is more representative of the key papers this way.
add_most_freq_author_papers <- function(df,...){  
  nauthor <- df %>% 
    pull(author_list) %>% 
    unlist() %>% 
    as_tibble() %>% 
    rename(authori = value) %>% 
    group_by(authori) %>% 
    summarise(authorfreq_within_group=n())
  df %>% 
    mutate(sum_author_freq = map_dbl(author_list, function(x) nauthor %>% filter(authori %in% x) %>% pull(authorfreq_within_group) %>% sum()))
} 


# ------------------------------------------------------------------
# graph methods
# lots of possible measures https://igraph.org/c/doc/igraph-Structural.html#centrality-measures
# we focus on closeness (inverse shortest path to all other verticies in the graph)
# this is most representative in terms of keywords, and most citing/cited by in terms of the references.
# see other interpretations in scimeetr
# these are somewhat slow for large graphs so we restrict our analysis to few of these measures.
add_group_closeness <- function(df, graph, ...){  
  graphi <-igraph::induced_subgraph(graph = graph, vids = df$PID, impl = "auto")
  close_score <- tibble(
    PID = igraph::V(graphi)$name,
    close_score = igraph::closeness(graphi, mode = "all", normalized = T)
  )
  df %>% 
    left_join(close_score, by = "PID")
}


# ==================================================================

# functions to give group summaries (individual functions, then summarise into one)

# ------------------------------------------------------------------

# GRAPH CENTRALITY

# ------------------------------------------------------------------
# papers with high closeness
papers_high_centrality <- function(df, np=5, graph){
  df %>% 
    add_group_closeness(graph = graph) %>% 
    arrange(-close_score, age) %>% 
    slice_head(n = np) %>% 
    add_ref() %>% 
    select(PID, ref, doi)
}

# ------------------------------------------------------------------

# CITATIONS

# ------------------------------------------------------------------
# papers most cited in all wos sources - slice_max(n_cited_allwos, age)
papers_most_cited_all <- function(df, np=5,...){
  df %>% 
    arrange(-n_cited_allwos, age) %>% 
    slice_head(n = np) %>% 
    add_ref() %>% 
    select(PID, ref, doi)
}

# ------------------------------------------------------------------
# papers most cited in all wos sources (age residuals) 
papers_most_cited_all_resid <- function(df, np=5,...){
  df %>% 
    arrange(-n_cited_allwos_resid, age) %>% 
    slice_head(n = np) %>% 
    add_ref() %>% 
    select(PID, ref, doi)
}

# ------------------------------------------------------------------
# papers most cited within group - frequency of occurence in citedPID_list 
papers_most_cited_within_group <- function(df, np=5,...){
  df %>% 
    add_group_most_cited() %>% 
    arrange(-cites_within_group, age) %>% 
    slice_head(n = np) %>% 
    add_ref() %>% 
    select(PID, ref, doi)
}

# ------------------------------------------------------------------
# papers most cited within group (age-residuals)
papers_most_cited_within_group_resid <- function(df, np=5, .k=min(round(nrow(df)/10), 10), ...){
  df2 <- df %>% add_group_most_cited()
  df %>% 
    add_column(n_cited_within_resid = age_based_residuals(df2, .ycol = "cites_within_group", .k = .k)) %>% 
    arrange(-n_cited_within_resid, age) %>% 
    slice_head(n = np) %>% 
    add_ref() %>% 
    select(PID, ref, doi)  
}

# ------------------------------------------------------------------
# papers citing most others within the group - sum of group's occurence in citedPID_list (likely review paper)
papers_citing_most_within_group <- function(df, np=5,...){
  df %>% 
    add_group_cites_most_group() %>% 
    arrange(-cites_of_group, age) %>% 
    slice_head(n = np) %>% 
    add_ref() %>% 
    select(PID, ref, doi)
}

# ------------------------------------------------------------------

# AUTHORS

# ------------------------------------------------------------------
# authors most frequently occuring within group - calculate from unlisted author list
authors_most_frequent <- function(df, np=5, ma=5, mn=3, ...){
  df %>% 
    pull(author_list) %>% 
    unlist() %>% 
    as_tibble() %>% 
    rename(author1 = value) %>% 
    group_by(author1) %>% 
    summarise(freq=n()) %>% 
    slice_max(freq, n=ma) %>% 
    # most frequent authors' top cited papers within groups
    mutate(authors_top_papers = map(author1, function(a1){
      df %>% 
        add_group_most_cited() %>% 
        filter(map_lgl(author_list, ~(a1 %in% .x))) %>% 
        arrange(-cites_within_group, age) %>% 
        slice_head(n = mn) %>% 
        add_ref() %>% 
        select(PID, ref, doi)
    })) %>% 
    unnest(authors_top_papers)
}
# ------------------------------------------------------------------

# KEYTERMS

# ------------------------------------------------------------------
# top keyterms (from (optionally) keywords, keyworads plus, titles and abstracts)
top_keywords <- function(df, nk= 15){
  df %>% 
    pull(keywords) %>% 
    strsplit("[;][ ]") %>% 
    unlist() %>% 
    as_tibble() %>% 
    rename(keywords = value) %>% 
    group_by(keywords) %>% 
    summarise(freq=n()) %>% 
    slice_max(freq, n=nk)
}

# ------------------------------------------------------------------
top_keywords_plus <- function(df, nk= 15){
  df %>% 
    pull(keywords_plus) %>% 
    strsplit("[;][ ]") %>% 
    unlist() %>% 
    as_tibble() %>% 
    rename(keywords_plus = value) %>% 
    group_by(keywords_plus) %>% 
    summarise(freq=n()) %>% 
    slice_max(freq, n=nk)
}

# ------------------------------------------------------------------
top_title_words <- function(df, nk= 15){
  df %>% 
    pull(title) %>% 
    strsplit("[ ]") %>% 
    unlist() %>% 
    as_tibble() %>% 
    rename(title_words = value) %>% 
    group_by(title_words) %>% 
    summarise(freq=n()) %>% 
    slice_max(freq, n=nk)
}
# ------------------------------------------------------------------
top_abstract_words <- function(df, nk = 15){
  df %>% 
    pull(abstract) %>% 
    strsplit("[ ]") %>% 
    unlist() %>% 
    as_tibble() %>% 
    rename(abstract_words = value) %>% 
    group_by(abstract_words) %>% 
    summarise(freq=n()) %>% 
    slice_max(freq, n=nk)
}

# could add bigrams there...

# ------------------------------------------------------------------

# JOURNALS

# ------------------------------------------------------------------
# journals most common
journals_most_common <- function(df, nj = 5){
  df %>%
    select(journal_iso) %>% 
    group_by(journal_iso) %>% 
    summarise(freq=n()) %>% 
    slice_max(freq, n=nj)
}

# ------------------------------------------------------------------
# cited journals most common in group
journals_most_commonly_cited <- function(df, nj = 5){
  df %>% 
    pull(cited_references_list) %>% 
    unlist() %>% 
    strsplit(", ") %>% 
    map_chr(~.x[3]) %>% 
    as_tibble() %>% 
    rename(journal_cited = value) %>% 
    group_by(journal_cited) %>% 
    summarise(freq=n()) %>% 
    slice_max(freq, n=nj)
}
