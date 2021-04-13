# flashLIT functions: term list revision

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# =============================

#' term_clean_stem_complete 
#' 
#' takes a character vector of terms and outputs a tibble with the original, cleaned, stemmed, and stem-completed values. Stem-completion requires a dictionary of terms, and where the term does not exist in the dictionary, NA will be output.
#'
#' @param terms character vector of terms
#' @param .dict stem completer dictionary
#'
#' @return a tibble with the original, cleaned, stemmed, and stem-completed values
#' @export
#'
#' @examples
term_clean_stem_complete <- function(terms, .dict = stemDict){
  tibble(
    original = terms,
    cleaned = clean_text(original),
    stem = stemR(cleaned),
    completed = stemCompletR(stem, .dict = .dict)
  )
}

# =============================

#' lookupDict
#' 
#' create a lookup dictionary from columns in a dataframe. COlumns should be one or more of the pre-cleaned and stemmed ('title', 'abstract', 'keywords', 'keywords_plus') or the original ('title_orig', 'abstract_orig', 'keywords_orig') columns. The latter will be cleaned (but not stemmed), the prior will remain as-is. 
#'
#' @param df the dataframe with the required selected columns
#' @param .cols a string with any combination of 'title', 'abstract', 'keywords', 'keywords_plus', 'title_orig', 'abstract_orig', 'keywords_orig'
#' @param sep.txt seperator string for text columns (default " ")
#' @param sep.kw seperator string for keyword columns (default "; ")
#' @param ngrams integer vector > 0, which ngrams are desired in the output (default 1:2 giving singles and bigrams)
#' @param ... other parameters passed to clean_text()
#'
#' @return a tibble with 
#' @export
#'
#' @examples
lookupDict <- function(df, .cols = 'title, abstract, keywords', sep.txt = " ", sep.kw = "; ", ngrams = 1:2, ...){
  
  # logicals for columns desired
  byti <- grepl(x = .cols, pattern = "title")
  byab <- grepl(x = .cols, pattern = "abstract")
  bykw <- grepl(x = .cols, pattern = "keywords")
  bykwp <- grepl(x = .cols, pattern = "keywords_plus")
  byti_o <- grepl(x = .cols, pattern = "title_orig")
  byab_o <- grepl(x = .cols, pattern = "abstract_orig")
  bykw_o <- grepl(x = .cols, pattern = "keywords_orig")
  
  # if already cleaned columns then take as is
  ti <- if(byti) {df %>% select(title) %>% pull() %>% str_split(pattern = sep.txt) %>% unlist()} else NULL
  abs <- if(byab) {df %>% select(abstract)%>% pull() %>% str_split(pattern = sep.txt) %>% unlist()} else NULL 
  kw <- if(bykw) {df %>% select(keywords) %>% pull() %>% str_split(pattern = sep.kw) %>% unlist()} else NULL 
  kwp <- if(bykwp) {df %>% select(keywords_plus) %>% pull() %>% str_split(pattern = sep.kw) %>% unlist()} else NULL
  
  # if not cleaned, then do basic cleaning
  ti_o <- if(byti_o) {df %>% select(title_orig) %>% clean_text(...) %>% str_split(pattern = sep.txt) %>% unlist()} else NULL
  abs_o <- if(byab_o) {df %>% select(abstract_orig) %>% clean_text(...) %>% str_split(pattern = sep.txt) %>% unlist()} else NULL 
  kw_o <- if(bykw_o) {df %>% select(keywords_orig) %>% pull() %>% str_split(pattern = sep.kw) %>% unlist() %>% clean_text(...) } else NULL 
  
  # then combine 
  map(ngrams, function(nx){
    x <- c(ti, abs, kw, kwp, ti_o, abs_o, kw_o) %>% 
        removePunctuation() %>% 
        .[which(.!="")] %>% 
        .[which(!is.na(.))] %>% 
        .[which(str_count(., pattern = " ") == (nx-1))] 
    
    if(length(x) > 0) {
      x %>% 
        ngram::ngram(n = nx) %>% 
        ngram::get.phrasetable() %>% 
        tibble() %>% 
        rename(word = ngrams) %>% 
        mutate(word = trimws(word))} 
  }) %>% 
    bind_rows() %>% 
    arrange(-freq)
      
}


# =============================

# create a document-term matrix, using the term-list as a dictionary
## to create these, we need to use tm to create a corpus, and them apply the term document matrix using a supplied dictionary of terms. 
## we can convert back to tidy format with broom::tidy 

# To create a dtm using a dictionary, we:
# y <- c("text", "this")
# x <- c("this is a text", "this another one")  # this should be a clean text!
# x <- tm::VCorpus(tm::VectorSource(x))
# dtm <- tm::DocumentTermMatrix(x, control = list(dictionary = y)) 
# inspect(dtm)

#' term_group_proportions
#' 
#' determine how many of the documents within a group contain a particular term (type = "terms") or at least one of a group of terms (type = "documents")
#'
#' @param df the bibliography dataframe
#' @param colnames which columns to search for the terms in
#' @param groupname the name of the group
#' @param communities a dataframe defining the communities
#' @param dict the dictionary terms (character vector)
#' @param type either "terms" or "documents" (see description)
#'
#' @return if type = "documents", a numeric giving the proportion of the documents containing at least one of the group terms, else, if type = "terms", a tibble giving the terms and the proportion of the documents that each term is found in.
#' @export
#'
#' @examples
term_group_proportions <- function(df, 
                                   colnames,
                                   groupname, 
                                   communities, 
                                   dict, 
                                   type = "documents"){
  # df <- DB1
  # groupname <- "2_2_3"
  # communities <- coms_summaries
  # dict <- T_floralcompetition$completed
  # colname <- enquo(colname)
  colnames <- enquos(colnames)
  group_PIDs <- communities %>% filter(group_id == groupname) %>% pull(group_PIDs) %>% .[[1]] %>% pull(PID)
  
  x <- df %>% 
    filter(PID %in% group_PIDs) %>% 
    select(!!!colnames) %>% 
    unite("words", !!!colnames, sep = " ") %>% 
    pull() %>% 
    tm::removePunctuation() 
  
  x <- tm::VCorpus(tm::VectorSource(x))
  dtm <- tm::DocumentTermMatrix(x, control = list(dictionary = dict)) 
  
  # proportion of documents with seperate term
  if(type == "terms") {broom::tidy(dtm) %>% group_by(term) %>% summarise(doc_prop = n()/length(x))}
  # proportion of documents with at least one of the terms
  if(type == "documents") {broom::tidy(dtm) %>% group_by(document) %>% summarise(doc_prop = n()) %>% nrow(.)/length(x)}
}


# =============================

#' group_top_keywords
#' 
#' extract the group top keywords from a community summary object
#'
#' @param community_summaries a community summary object dataframe
#' @param keyword_types a vector of the column names specifiying the keyword types desired (keywords, keywords_plus, title_words, abstract_words)
#'
#' @return a tibble with the keywords and their frequency (i.e. maximum frequency across the keyword columns)
#' @export
#'
#' @examples
group_top_keywords <- function(community_summaries, ...){
  keyword_types <- enquos(...)
  community_summaries %>% 
    select(group_top_keywords) %>% 
    unnest_longer(col = group_top_keywords) %>% 
    unnest(cols = group_top_keywords) %>% 
    select(freq, !!!keyword_types) %>% 
    pivot_longer(cols = c(!!!keyword_types), names_to = "type", values_to = "word") %>%
    remove_missing(na.rm = TRUE) %>% 
    group_by(word) %>% 
    summarise(freq = max(freq)) 
}

#======================================


#' suggest_terms
#' 
#' suggest further terms for consideration in the term lists. The user should then select which ones they want by saving into a new vector and deleting unwanted terms.
#'
#' @param coms_summaries a community summary object 
#' @param existing_words a list of existing words
#' @param ... columns of the group_top_keyword columns desired for keyword suggestion (title_words, abstract_words, keywords, keywords_plus)
#' @param above_freq only return words with a maximum frequency in any of the columns larger than this.
#'
#' @return a formatted string vector (with quotes) of suggested words and the original words for consideration in the term lists. 
#' @export
#'
#' @examples
suggest_terms <- function(coms_summaries, existing_words, ..., above_freq = 1){
  .cols <- enquos(...)
  group_top_keywords(coms_summaries, !!!.cols) %>% 
    filter(freq > 1) %>% 
    arrange(-freq) %>% 
    pull(word) %>% 
    c(., existing_words) %>% 
    sort() %>% 
    shQuote( type = "csh") %>% 
    paste( collapse=", ")
}
