# flashLIT group summaries

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

# group summaries 

# given a group and a dataframe (and for some summaries, a graph object)
# most of these can come as 'within group' (i.e. in terms of the papers within the group) or as externally figured (e.g. number of citations in WoS+)
# in our case, we really want to characterise the groups so we focus on the within group measures

# CITATIONS
# papers most cited in all wos sources - slice_max(n_cited_allwos, year2)
# papers most cited in all wos sources (age residuals) - slice_max(n_cited_allwos_resid, year2)
# papers most cited within group - frequency of occurence in citedPID_list * also version of this age-residuals
# papers citing most others within the group - sum of group's occurence in citedPID_list (likely review paper)
# AUTHORS
# authors most frequently occuring within group - calculate from unlisted author list
# most frequent authors' top cited papers within groups
# KEYTERMS
# top keyterms (from (optionally) keywords, keywords plus, titles and abstracts)
# JOURNALS
# journals most common
# cited journals most common in group

# we possibly also want to generate some group metrics
# are the groups well-connected or?
# these might be useful in determining our confidence in allocating rogue papers. 
# but we might not do that here just yet.

# non group-specific frequencies required 
# include the age based residuals for total citations
DB1 <- DB1 %>% 
  add_column(n_cited_allwos_resid = age_based_residuals(DB1, "n_cited_allwos"))

# visual check of what this is doing
# ggplot(DB1) + 
#   geom_point(aes(x=age, y = n_cited_allwos, color = n_cited_allwos_resid), pch = 16) +
#   geom_line(aes(x=age, y = mgcv::gam(n_cited_allwos ~ s(age, k=10), data = DB1, family = "poisson")$fitted.values), color = 'red') +
#   scale_color_viridis_c()

# now, for each community group level, including the whole DB, we need to map over the group summaries

coms_summaries <- coms_DB1 %>% 
  pivot_longer(cols = c(names(coms_DB1)[-1]), names_to = "group_level", values_to = "group_id") %>% 
  group_by(group_id) %>% 
  nest(group_PIDs = c(PID)) %>% 
  mutate(n = map_int(group_PIDs, nrow))

# apply group summary functions - this can take some time if including group centrality  
# we do this in sections, to aid debugging.

coms_summaries <- coms_summaries %>% 
  mutate(group_top_keywords = map(group_PIDs, function(x){
    .df <- DB1 %>% filter(PID %in% pull(x)) 
    list(
      top_keywords = top_keywords(.df),
      top_keywords_plus = top_keywords_plus(.df),
      top_title_words = top_title_words(.df),
      top_abstract_words = top_abstract_words(.df)
    )
  }
  )) 

coms_summaries <- coms_summaries %>% 
  mutate(group_top_journals = map(group_PIDs, function(x){
    .df <- DB1 %>% filter(PID %in% pull(x)) 
    list(
      journals_most_common = journals_most_common(.df),
      journals_most_commonly_cited = journals_most_commonly_cited(.df)
    )
  }
  )) 

coms_summaries <- coms_summaries %>% 
  mutate(group_top_papers = map(group_PIDs, function(x){
    .df <- DB1 %>% filter(PID %in% pull(x)) 
    list(
      papers_high_centrality = tryCatch(papers_high_centrality(.df, graph = graph_DB1), error = function(e) e),
      papers_most_cited_all = tryCatch(papers_most_cited_all(.df), error = function(e) e),
      papers_most_cited_within_group = tryCatch(papers_most_cited_within_group(.df), error = function(e) e),
      papers_citing_most_within_group = tryCatch(papers_citing_most_within_group(.df), error = function(e) e),
      authors_most_frequent = tryCatch(authors_most_frequent(.df), error = function(e) e)
    )
  }
  )) 


# these models sometimes fail, so we do them seperately with cautions for errors (warnings will likely also occur)
coms_summaries <- coms_summaries %>% 
  mutate(group_top_papers_resid = map(group_PIDs, function(x){
    if (nrow(x)>30){
      .df <- DB1 %>% filter(PID %in% pull(x)) 
      list(
        papers_most_cited_all_resid = tryCatch(ifelse(nrow(x)>30, papers_most_cited_all_resid(.df), NA), error = function(e) e),
        papers_most_cited_within_group_resid = tryCatch(ifelse(nrow(x)>30, papers_most_cited_within_group_resid(.df), NA), error = function(e) e)
      )
    } else {
      list("Fewer than 30 documents, summary model not completed")
    }
  }
  )) 

# endscript