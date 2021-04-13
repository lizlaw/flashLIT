# functions and summaries of a scimeetr object


#' group summary of a (scimeetr) object 
#'
#' @param .x a scimeetr object
#'
#' @return a tibble with group names, level (of splitting), and size of each group.
#' @export
#'
#' @examples
grSummary <- function(.x, kw = FALSE){
  df <- tibble(
    grNames = names(.x),
    grLevel = names(.x) %>% stringr::str_count(pattern = "_"),
    grSize = map_dbl(.x, ~nrow(.x$dfsci))
  )
  if(kw) df <- df %>% 
      add_column(grTag = map(.x, ~.x$tag),
                 grKeywords = map(.x, ~.x$kw),
                 grTitle = map(.x, ~.x$ti),
                 grAbstract = map(.x, ~.x$ab))
  return(df)
}

#' iterative split of a scimeetr object into groups
#' 
#' will iteratively apply scimeetr::scimap to a scimeetr object until the stopping rules are satisfied.
#' rules are specified as both groupsize AND number of groups need to reach targets, or the end is triggered by 
#' a timeout.
#'
#' @param x a scimeetr object (can already be split)
#' @param coupling style of coupling. see ?scimeetr::scimap.
#' @param grouptarget number of groups desirable
#' @param mingroupsize minimum size of the groups
#' @param maxgroupsize maximum desired size of the groups
#' @param timeout maximum time (minutes) to iterate for
#'
#' @return
#' @export
#'
#' @examples
scimap_iterative  <- function(x, coupling = 'bickec', 
                              grouptarget = 20, mingroupsize = 5, maxgroupsize = 30,
                              timeout = 5){
  timeStart <- Sys.time()
  timeElapsed <- 0
  grSum <- grSummary(x)
  ngroups <- grSum %>% filter(grLevel == max(grLevel)) %>% nrow()
  maxGrSize <- grSum %>% filter(grLevel == max(grLevel)) %>% pull(grSize) %>% max
  
  while((ngroups  < grouptarget | maxGrSize > maxgroupsize) & !(timeElapsed > timeout)){
    
    x <- scimeetr::scimap(x, coupling_by = coupling, min_com_size = mingroupsize)
    
    grSum <- grSummary(x)
    ngroups <- grSum %>% filter(grLevel == max(grLevel)) %>% nrow()
    maxGrSize <- grSum %>% filter(grLevel == max(grLevel)) %>% pull(grSize) %>% max
    timeElapsed <- difftime(Sys.time(), timeStart, units = "mins")
  }
  attr(x,'STATUS') <- c(ngroups = ngroups, maxGrSize = maxGrSize, timeElapsed = timeElapsed)
  return(x)
}  

### APPLY:

# scimeetr import 
dd <- "./data/raw_WoS_20201105/as_txt/"
scimeetr_list <- scimeetr::import_wos_files(dd)

sciGroups <- scimap_iterative(scimeetr_list)  ## takes ~4 minuites for this data
attr(sciGroups,'STATUS') 

grSummary(sciGroups)
grSummary(sciGroups) %>% pull(grLevel) %>% max() ## goes up to 4 levels deep
grSummary(sciGroups) %>% filter(grLevel==4) %>% pull(grSize) %>% hist()

# From here, we may want to look at the keywords that define the groups (tags), by group level to help inform some of the terms
# note the full keyword, title, and abstract are the frequency of terms in these. We will use these next, in the terms

sciSum <- grSummary(sciGroups, kw = TRUE)

sciSum %>% 
  select(grNames:grTag) %>% 
  filter(grLevel <3) %>% 
  unnest_wider(grTag)

sciSum %>% 
  filter(grLevel >3) %>% 
  unnest(grTag) %>% 
  pull(grTag) %>% 
  table()