# flashLIT functions: import and clean

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# =============================

# modify https://github.com/MaximeRivest/scimeetr/blob/master/R/import_wos_files.R
# only create the DF
# add columns with filename and foldername

#' import_wos_files_DFonly
#' 
#' A modification of scimeetr::import_wos_files()  https://github.com/MaximeRivest/scimeetr/blob/master/R/import_wos_files.R
#' This version only creates the DF (does not convert to a scimeetr list) and also includes the filename and folder
#'
#' @param files_directory 
#'
#' @return a dataframe 
#' @export
#'
#' @examples
import_wos_files_DFonly <-
  function (files_directory) 
  {
    dfsci <- NULL
    dfsci_temp <- NULL
    
    folder_content <- list.files(files_directory)
    files_quantity <- length(folder_content)
    
    for (files in 1:files_quantity) {
      full_file_path <- paste(files_directory, folder_content[files], 
                              sep = "")
      v_char <- suppressWarnings(readLines(full_file_path, encoding = "UTF-8"))
      v_char <- iconv(v_char, from = "UTF-8", to = "ASCII", 
                      sub = "")
      v_char <- stringr::str_replace(v_char, "^[null]+", "")
      tab_count <- stringr::str_count(v_char[], '\t')
      good_lines <- c(1, which(tab_count == max(tab_count)))
      dfsci_temp <- read.table(text = v_char[good_lines], header = T, quote = "",
                               fileEncoding = 'ASCII',
                               row.names = NULL,
                               comment.char = "",
                               stringsAsFactors = F,
                               sep = "\t")
      dnames <- names(dfsci_temp)[-1]
      dfsci_temp <- dfsci_temp[-length(dfsci_temp)]
      names(dfsci_temp) <- dnames
      
      dfsci_temp <- dfsci_temp %>% 
        add_column(file = folder_content[files]) %>% 
        add_column(folder = files_directory)
      dfsci <- rbind(dfsci, dfsci_temp)
    }
    
    dfsci <- dfsci %>% distinct()
    
    if(sum(is.na(dfsci$CR)) == nrow(dfsci)) {
      warning("The field CR (cited reference) contains no references. \nAlmost everything is scimeetr will break without that field. \nWhen you download records from WOS or Scopus, please carefully follow that steps at: \nhttps://github.com/MaximeRivest/scimeetr#how-to-get-bibliometric-data",call. = F)
    } 
    
    return(dfsci)
  }

# ==============================


#' import_wos_ref_files
#' 
#' A modification of scimeetr::import_wos_files()  https://github.com/MaximeRivest/scimeetr/blob/master/R/import_wos_files.R
#' This version only creates the DF (does not convert to a scimeetr list) and also includes the filename and folder
#' AND is specific to importing the WoS output for the reference list of a paper(s).
#'
#' @param files_directory 
#'
#' @return a dataframe 
#' @export
#'
#' @examples
import_wos_ref_files <-
  function (files_directory) 
  {
    dfsci <- NULL
    dfsci_temp <- NULL
    
    folder_content <- list.files(files_directory)
    files_quantity <- length(folder_content)
    
    for (files in 1:files_quantity) {
      full_file_path <- paste(files_directory, folder_content[files], 
                              sep = "")
      v_char <- suppressWarnings(readLines(full_file_path, encoding = "UTF-8"))
      v_char <- iconv(v_char, from = "UTF-8", to = "ASCII", 
                      sub = "")
      v_char <- stringr::str_replace(v_char, "^[null]+", "")
      tab_count <- stringr::str_count(v_char[], '\t')
      good_lines <- which(tab_count == max(tab_count))
      dfsci_temp <- read.table(text = v_char[good_lines], header = T, quote = "",
                               fileEncoding = 'ASCII',
                               row.names = NULL,
                               comment.char = "",
                               stringsAsFactors = F,
                               sep = "\t")
      dfsci_temp <- dfsci_temp %>% 
        add_column(file = folder_content[files]) %>% 
        add_column(folder = files_directory)
      dfsci <- rbind(dfsci, dfsci_temp)
    }
    
    dfsci <- dfsci %>% distinct()
    
    if(sum(is.na(dfsci$CR)) == nrow(dfsci)) {
      warning("The field CR (cited reference) contains no references. \nAlmost everything is scimeetr will break without that field. \nWhen you download records from WOS or Scopus, please carefully follow that steps at: \nhttps://github.com/MaximeRivest/scimeetr#how-to-get-bibliometric-data",call. = F)
    } 
    
    return(dfsci)
  }

# ==============================

#' convert_names
#' 
#' uses the revtools tags to convert between formats. NOTE this process can remove columns that are duplicated.
#'
#' @param df bibliography database e.g. import table from revtools or scimeetr
#' @param from format type, i.e. (rev)bib, ris, ris_write, wos, medline
#' @param to format type, see 'from'
#'
#' @return a dataframe with revised column names
#' @export
#'
#' @examples
convert_names <- function(df, from, to){
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
    df <- df %>% select(any_of(c(unname(unlist(lookup['bib'])), "folder", "file")))
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
    df <- df %>% select(any_of(c(unname(unlist(lookup['ris'])), "folder", "file")))
  }
  return(df)
}

# =============================

# quick check function

#' quick_check
#' 
#' rapid check of dataframe to identify how many NA entries in the different columns, inlcuding ones that may be hidden as text or missing data without being explicitly NA. 
#'
#' @param df dataframe to check
#'
#' @return tibble inlcuding column names, class, na (is.na), and 'hidden_na' (identifies all no data instances, not just the NA)
#' @export
#'
#' @examples
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


# =============================

# checking for potential duplicates 

#' fuzzymatch_potential_duplicates
#' takes in two string vectors, and fuzzymatches them using tidystringdist (wrapper for stringdist). Two methods are selected, "jw", and "osa". It will return matched entries where 0 < jw < 0.2, or  0 < osa < 25 & jw < 0.3. This is reasonably generous and should find all the potential matches.
#'
#' @param V1 A character vector
#' @param V2 Another character vector (or the same as V1)
#'
#' @return a tibble with V1, V2, osa and jw columns
#' @export
#'
#' @examples
fuzzymatch_potential_duplicates <- function(V1, V2=NULL){
  if(is.null(V2)) {comb <- tidy_comb_all(c(V1))} else {expand_grid(V1 = V1, V2 = V2)}
  tms <- tidystringdist::tidy_stringdist(comb, method = c("osa", "jw"))
  bind_rows(
    tms %>% filter(jw > 0 & jw < 0.2),
    tms %>% filter(osa > 0 & osa < 25 & jw < 0.3) 
  ) %>% filter(!duplicated(.))
}

#' checkmatch_potential_duplicates
#'
#' @param checktext a text string to find in the in_col
#' @param dfs a list of databases to check in
#' @param in_col the name of the input col
#' @param ... additional output columns desired
#'
#' @return a tibble containing the rows with the checktext in the in_col
#' @export
#'
#' @examples
checkmatch_potential_duplicates <- function(checktext, dfs, in_col, ...){
  in_col <- enquo(in_col)
  out_cols <- enquos(...)
  map_dfr(dfs, function(x) x %>% filter(grepl(checktext, !!in_col)) %>% select(!!in_col, !!!out_cols))
}
# checkmatch_potential_duplicates("Large-scale monitoring of effects of clothianidin-dressed oilseed rape seeds", list(DB1, KR2_refs), title_orig, year, author) 



# =============================

# Cleaning will be applied per column. First clean, then stem.

#' clean_text
#' 
#' clean a text column in the dataframe. Converts all to lower case, then removes numbers, punctuation, and (custom) stopwords, then trims whitespace in and around text. If a seperation is specified, text will first be split by this character, then recombined using the same character.
#'
#' @param x a dataframe column
#' @param sep a seperation character, or the default NULL
#'
#' @return c leaned dataframe column
#' @export
#'
#' @examples
clean_text <- function(x, sep = NULL, .stopwords = meaningless_words){
  if(!is.null(sep)){
    x <- x %>% str_split(., pattern = sep) %>% .[[1]]
  }
  x <- x %>% 
    tolower() %>% 
    removeNumbers() %>% 
    removePunctuation(preserve_intra_word_dashes = TRUE) %>% 
    str_replace_all("-", " ") %>% 
    removeWords(words = .stopwords) %>% 
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

# ================================================================================

# functions to stem, them complete using the dictionary.
# stemming with SnowballC::wordStem
# to be more interpretable, we would like to complete the words with the stem.
# to do that, we need to create a stem dictionary from the entire columns
# I've used ngram here, but could have used tidytext.

#' stem_dictionary
#' 
#' creates a term dictionary for the stemming. 
#'
#' @param df the dataframe
#' @param ... the columns to be used in the stem dictionary
#'
#' @return a two column tibble, giving each 'word' and the associated 'stem'
#' @export
#'
#' @examples
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


#' stemR
#' 
#' stems a column. this makes term comparisons more consistent.
#'
#' @param x a text column.
#' @param token_split a seperator character to use to split and then collapse tokens in each entry, if not NULL.
#'
#' @return a column of stemmed words
#' @export
#'
#' @examples
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

#' stemCompletR
#' 
#' complete stemmed words using the specified dictionary. this makes stemmed words more human 'readable'.
#'
#' @param x a column of stemmed words
#' @param .dict a specified dictionary
#' @param token_split a seperator character to use to split and then collapse tokens in each entry, if not NULL.
#'
#' @return a column of completed words.
#' @export
#'
#' @examples
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

# fill in missing keywords ==============================================================

# keywords can be filled in using common words from titles and abstracts (or other specified columns)
# I've used ngram here, but could have used tidytext.

#' keywords_from_columns
#' 
#' keywords can be filled in using common words from titles and abstracts (or other specified columns). this function takes the specified columns, and selects the <10 most frequent words and phrases (up to 3 words) in the data. 
#'
#' @param df the dataframe with keywords to be filled
#' @param .cols columns of text to be used to find keywords (e.g. title and abstract)
#' @param .sep a character to seperate the resulting keywords with.
#' @param .n integer, the number of keywords desired to be returned
#'
#' @return
#' @export
#'
#' @examples
keywords_from_columns <- function(df, .cols, .sep = "; ", .n = 10){
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
    slice_head(n=.n) %>% 
    pull(ngrams) %>% 
    trimws() %>% 
    paste(collapse = .sep) %>% 
    na_if(y="NA") %>% 
    na_if(y="")
}

# note:
# KeyWords Plus are words or phrases that frequently appear in the titles of an article's references, but do not appear in the title of the article itself... KeyWords Plus enhances the power of cited-reference searching by searching across disciplines for all the articles that have cited references in common. (https://support.clarivate.com/ScientificandAcademicResearch/s/article/KeyWords-Plus-generation-creation-and-changes?language=en_US)
# for these, we can use the references that are both cited and in our database to generate these, as the cited ref does not contain the title.
# they will differ from the original keywords plus, but approximate them.

# extract the cited references ==============================================================

#' cited_references_tolist
#' 
#' extract the cited references into a list, creating a 'CID' identifier, and retaining PID information
#' default seperates citation elements with "; "
#'
#' @param df the dataframe to draw the cited_references from
#' @param padCID number of zeros to pad the resulting CID with
#' @param .sep characters for seperating the citation elements
#'
#' @return a tibble, with CID, fromPID and the cited_reference entry
#' @export
#'
#' @examples
cited_references_tolist <- function(df, padCID = 6, .sep = "; "){
  x <- df %>% 
    select(fromPID = PID, cited_references) %>% 
    mutate(cited_references = map(cited_references, ~str_split(.x, pattern = .sep))) %>% 
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


#' cited_references_cleanlist
#' 
#' cleans cited reference list ready to match against a reference dataframe
#' note, this pre-removes all entries unlikely to gain a match.
#'
#' @param df a tibble created by cited_references_tolist()
#'
#' @return a tibble with CID, author, year, journal_iso, volume, pages, doi
#' @export
#'
#' @examples
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

#' cited_reference_comparelist
#'
#' @param df a dataframe to draw the citations from. uses PID, author, year, journal_iso, volume, pages, doi columns, and cleans these to the same format as the bibliography (i.e. strips author to an uppercase single name with single initial, appends V and P to volume and pages, and cleans DOI)
#'
#' @return a dataframe with cleaned PID, author, year, journal_iso, volume, pages, doi columns
#' @export
#'
#' @examples
cited_reference_comparelist <- function(df){
  df %>% 
    select(PID, author, year, journal_iso, volume, pages, doi) %>% 
    mutate(author = map_chr(author, ~str_split(.x, pattern = "; ")[[1]][[1]] %>% str_remove(",")), 
           author = map_chr(author, ~str_split(.x, pattern = " ")[[1]][[1]] %>% toupper()),
           volume = paste0("V", volume),
           pages = paste0("P", pages),
           doi = map_chr(doi, ~ifelse(is.na(.x), .x, paste("DOI", str_split(.x, pattern = " ")[[1]])))) 
}

#' cited_reference_match
#' 
#' match a cleaned cited references list with a cleaned comparison list, by specified columns
#'
#' @param crlist the output from cited_references_cleanlist() 
#' @param comparelist output from cited_reference_comparelist()
#' @param ... the columns that need to be matched (exactly) for the match to be acceptable.
#' @param filter_missing logical defaulting to TRUE, to remove NA entries before matching (as these will match all other NAs)
#'
#' @return a joined dataframe containing the matches
#' @export
#'
#' @examples
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

# finally, cleaning the author list. Authors are commonly given various initials, so we want to clean this to AUTHOR I, AUTHOR2 I, ...
# capitalised, surname, single initial.

#' clean_author_list
#' 
#' cleaning the author list. Authors are commonly given various initials, so we want to clean this to have a surname, and single initial only, and capitalised.
#'
#' @param x the author column of a dataframe
#' @param .sep the character seperating the author names
#'
#' @return a cleaned column of a dataframe
#' @export
#'
#' @examples
clean_author_list <- function(x, .sep = "; ", .sep2 = ", "){
  x %>% 
    map( ~str_split(.x, pattern = .sep)[[1]]) %>% 
    map( ~toupper(.x) %>% 
          stripWhitespace() %>% 
          trimws()) %>% 
    map( function(x){
      y <- str_split(x, pattern = .sep2) 
      map(y, ~paste(.x[1], str_sub(.x[2], 1, 1))) %>% # extract the surname and the first initial
        unlist() %>% 
        removePunctuation(preserve_intra_word_dashes=TRUE)})
}