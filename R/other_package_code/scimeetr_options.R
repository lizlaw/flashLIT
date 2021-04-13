# scimeetr options example

# Example: Mapping the evidence for impacts of managed bees on wild bees

# Search terms: 

# *Managed*: (African* NEAR bee) OR Apis OR Bombus OR "bumble bee"
# OR bumblebee* OR "honey bee" OR honeybee OR ((introduc* OR inva* OR non-native
# OR nonnative OR commercial OR exotic OR feral OR managed) NEAR (bee OR
# pollin*)) 

# *Native*: (((cavity OR ground) NEAR nesting) OR (native OR solitary
# OR wild)) NEAR (bee OR pollin*) 

# *Interaction*: pollinat* OR network* OR “niche
# overlap” OR “partitioning” OR interact* OR competit* OR facilitat* OR
# mutualis* OR “resource limitation” OR hybridization OR introgression OR
# dependence OR assemblag* OR overlap OR spillover OR impact* 

# Search Web of Science: TOPIC searches Title, Abstract, Keywords (supplied by the author),
# and Keywords Plus (algorithm extraction of expanded terms stemming from the
# cited references)
#
# TOPIC: ((African*  NEAR bee)  OR Apis  OR Bombus  OR "bumble bee"  OR
# bumblebee*  OR "honey bee"  OR honeybee  OR ((introduc*  OR inva*  OR
# non-native  OR nonnative  OR commercial  OR exotic  OR feral  OR managed)
# NEAR (bee  OR pollin* ))) AND TOPIC: ((((cavity  OR ground)  NEAR nesting)  OR
# (native  OR solitary  OR wild))  NEAR (bee  OR pollin* )) AND TOPIC:
# (pollinat*  OR network*  OR “niche overlap”  OR “partitioning”  OR interact*
# OR competit*  OR facilitat*  OR mutualis*  OR “resource limitation”  OR
# hybridization  OR introgression  OR dependence  OR assemblag*  OR overlap  OR
# spillover  OR impact*) 

# Timespan: All years. Indexes: SCI-EXPANDED, SSCI,A&HCI, ESCI. 
# Search date: 05.Nov.2020 
# Results: 2,400 Export full record and cited references in Tab-delimited (Win, UTF-8) format (for scimeetr)

# --------------------------------------------------------------------------------------------------------------------------------------------------------

# # install required packages
# remotes::install_github("MaximeRivest/scimeetr")

# load packages
library(tidyverse) # programming
library(scimeetr) # bibiometric analysis and determination of sub-communities

# import to scimeetr object -----------------------------------------------
dd <- "./data/raw_WoS_20201105/as_txt/"
scimeetr_list <- import_wos_files(dd)
scimeetr_list
summary(scimeetr_list)

# In the dfsci, the 68 variables imported from WoS and their frequency of non NA data:
fieldtags <- read_csv("wos_fieldtags.csv")
indf <- summarise_all(scimeetr_list$com1$dfsci, function(x) sum(!is.na(x))) %>% 
  t() %>% 
  as_tibble(rownames = "TAG") %>% 
  select(TAG, "is_not_NA" = V1) %>% 
  right_join(., fieldtags) %>% 
  filter(is_not_NA > 0)

indf %>% print(n = nrow(indf))

# # Use scimeetr to define sub-communities -----------------------------------------------

# create nested sub-communities - we'll do this twice to start, using the recommended bickecticjoc method.

sci_bktj30 <- scimap(scimeetr_list, coupling_by = 'bickecticjoc', community_algorithm = 'louvain', min_com_size = 30)
scisub_bktj30 <- scimap(sci_bktj30, coupling_by = 'bickecticjoc', community_algorithm = 'louvain', min_com_size = 30)

saveRDS(scisub_bktj30, "data/raw_WoS_20201105/scisub_bktj30.RDS")
scisub_bktj30 <- readRDS("data/raw_WoS_20201105/scisub_bktj30.RDS")

summary(scisub_bktj30)
plot(summary(scisub_bktj30, com_size = 30))

# define a custom 'deepdive' function ---------------------------------------------------------
deep_dive <- function(object, kw = 10, kr = 10, mr = 3, ...){
  
  out <- list()
  
  out$Overview <- list(
    str_c("Number of papers: ", nrow(object[[1]]$dfsci)),
    str_c("Number of communities: ", length(object)),
    object %>% map_dbl(function(x) nrow(x$dfsci)) %>% .[order(names(.))]
  )
  
  out$Word_frequencies <- list(
    tags = object[[1]]$tag,
    keywords = data.frame(
      "key_words" = object[[1]]$kw$ID[1:kw],
      "key_words_de" = object[[1]]$de$ID[1:kw],
      "title_words" = object[[1]]$ti$ID[1:kw],
      "abstract_words" = object[[1]]$ab$ID[1:kw], 
      "authors" = object[[1]]$au$ID[1:kw],
      "cited_refs" = object[[1]]$cr$ID[1:kw], 
      stringsAsFactors = F)
  )
  
  # prepare function for ranking papers
  get_pubs <- function(rlist){
    .kr <- ifelse(rlist %in% "by_expert_LC", ceiling(kr/mr), kr)
    scilist(object, k = .kr, m = mr, reading_list = rlist) %>% 
      .[[1]] %>% 
      select(publication) %>% 
      add_column(type = rlist) 
  }
  
  title_tab <- object[[1]]$dfsci %>% 
    select(RECID, Title = TI) %>% 
    mutate(publication = map_chr(RECID, function(x) x %>%  
                                   stringr::str_replace_all(', DOI.*', '') %>% 
                                   stringr::str_replace_all('V(?=[0-9]{1,6})', '') %>% 
                                   stringr::str_replace_all('P(?=[0-9]{1,6})', ''))) %>% 
    select(-RECID)
  
  cited_tab <- object[[1]]$dfsci %>% 
    select(RECID, TimesCited = Z9, Title = TI, PY, EA) %>% 
    mutate(Year = pmap_dbl(list(PY, EA), function(py,ea) ifelse(!is.na(py), py, ea %>% str_extract('[0-9]{4}') %>% as.numeric()))) %>% 
    select(-PY, -EA)
  
  cited_tab <-  cited_tab %>% 
    add_column(residuals = mgcv::gam(TimesCited ~ Year, data = cited_tab, family = "poisson")$residuals)
  
  cited_most <- bind_rows(
    top_n(cited_tab, kr, TimesCited) %>% add_column("cited_most" = 1, "cited_resid" = 0),
    top_n(cited_tab, kr, TimesCited) %>% add_column("cited_most" = 0, "cited_resid" = 1)
  ) %>%
    group_by(RECID, Title, Year, TimesCited, residuals) %>% 
    summarise_at(vars(cited_most:cited_resid), sum) %>% 
    arrange(desc(residuals), desc(TimesCited))
  
  out$Articles_highcited_in <- list(
    description = str_c("Highest cited papers in the group"),
    publist = cited_most
  )
  
  out$Articles_highcited_by <- list(
    description = str_c("Highest cited papers by the group"),
    publist = map(c("core_papers", "core_residual"), get_pubs) %>% 
      bind_rows() %>% 
      with(., table(publication, type)) %>% 
      as_tibble() %>% 
      pivot_wider(names_from = "type", values_from ="n") %>% 
      arrange(desc(core_papers), desc(core_residual)) %>% 
      left_join(title_tab, by = "publication")
  )
  
  out$Articles_experts <- list(
    description = str_c("Highest cited papers by highly cited authors"),
    publist = map(c("by_expert_LC"), get_pubs) %>% 
      bind_rows() %>% 
      with(., table(publication, type)) %>% 
      as_tibble() %>% 
      pivot_wider(names_from = "type", values_from ="n") %>% 
      arrange(desc(by_expert_LC)) %>% 
      left_join(title_tab, by = "publication")
  )
  
  out$Aricles_review <- list(
    description = str_c("Papers likely to provide a good overview of the category"),
    publist = map(c("cite_most_others", "connectness"), get_pubs) %>% 
      bind_rows() %>% 
      with(., table(publication, type)) %>% 
      as_tibble() %>% 
      pivot_wider(names_from = "type", values_from ="n") %>% 
      arrange(desc(cite_most_others), desc(connectness)) %>% 
      left_join(title_tab, by = "publication")
  )
  
  return(out)
}

# apply function -----------------------------------------------
comlist <- sort(names(scisub_bktj30))

deep_dive_results <- tibble(
  community = comlist,
  dd_results = map(comlist, function(x) deep_dive(dive_to(scisub_bktj30, aim_at = x)))
)

deep_dive_results$dd_results[[1]]
