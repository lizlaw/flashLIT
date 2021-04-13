# flashLIT functions: categorise papers 

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

term_counts <- function(df, dict, colnames = c(title, abstract, keywords)){
  # df <- DB1
  # dict <- T_floralcompetition$completed
  colnames <- enquos(colnames)
  
  x <- df %>% 
    select(!!!colnames) %>% 
    unite("words", !!!colnames, sep = " ") %>% 
    pull() %>% 
    tm::removePunctuation() 
  
  
  x <- tm::VCorpus(tm::VectorSource(x))
  for (i in 1:length(x)) {
    x[[i]]$meta$id <- df$PID[i]
  }
  dtm <- tm::DocumentTermMatrix(x, control = list(dictionary = dict)) 
  
  y <- broom::tidy(dtm) %>% 
    group_by(document) %>% 
    summarise(term_count = n())
  
  tibble(
    PID = df$PID
  ) %>%
    left_join(y, by = c("PID" = "document")) %>% 
    pull(term_count) %>% 
    replace(is.na(.), 0)
}