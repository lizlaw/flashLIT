# flashLIT stratified sample

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

# strategically sample
# sample at least one from all groups, and fill the remainder with up to 10% of the nrows 
s1 <- DB1 %>% 
  left_join(coms_DB1, by = "PID") %>% 
  select(PID, cl_5)  %>% 
  group_by(cl_5) %>% 
  slice_sample(n=1)
# gives 106 items

s2 <- DB1 %>% 
  left_join(coms_DB1, by = "PID") %>% 
  filter(! PID %in% s1$PID) %>% 
  select(PID, cl_5)  %>% 
  group_by(cl_5) %>% 
  slice_sample(prop=(nrow(DB1)*0.1 - nrow(s1))/nrow(DB1))
# gives an additional 82

stratsamp <- DB1 %>% 
  left_join(coms_DB1, by = "PID") %>% 
  mutate(select10 = ifelse(PID %in% c(s1$PID, s2$PID), 1, 0)) %>% 
  select(PID, journal_iso, year, author, keywords, keywords_orig, keywords_plus, title_orig, abstract_orig, select10) 
