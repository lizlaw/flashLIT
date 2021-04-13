# flashLIT categorise papers 

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

# here we aim to categorise papers into the 3 primary streams (foraging, nesting, and PPP) (with a 'other' or undefined)

# inclusion into a stream ---
# papers should be 'in' a stream if they include several instances of the stream keywords
# if they include high instances of several streams' keywords, then they can be potentially classified as both.
# if they do not mention many of any streams' keywords, then they can be allocated 'other'
# We can identify some papers as 'definitely' in the group (high instances of a single stream keywords, and verified)
# others we can use e.g. random forest models OR some other distance type methods to assign likelihoods to the other observations.

# group characterisation ---
# group characterisation should be at level 4 groups only (can be informed by higher groups, but these can often be multi-channeled)

# =====================================================================

# first, we need to create a dtm for each stream's words, again, here we only use the (cleaned) title, abstract, keywords of each document
# and cound how many times any of the terms are used in any of the components
DB1_dtms <- tibble(
  PID = DB1$PID,
  floralcompetition_in_TAK = term_counts(df = DB1, colnames = c(title, abstract, keywords), dict = T2_floralcompetition),
  nestcompetition_in_TAK = term_counts(df = DB1, colnames = c(title, abstract, keywords), dict = T2_nestcompetition),
  ppp_in_TAK = term_counts(df = DB1, colnames = c(title, abstract, keywords), dict = T2_ppp)
)  

dtmplot <- DB1_dtms %>% 
  left_join(coms_DB1, by = "PID")

# examining over the smallest group possible -----------------------
# suggests a few groups might be easily defined this way, but not many

ggplot(dtmplot) +
  geom_point(aes(x=floralcompetition_in_TAK, y = nestcompetition_in_TAK, color = ppp_in_TAK), alpha = 0.2) +
  facet_wrap(vars(cl_4)) +
  scale_color_viridis_c() +
  geom_hline(aes(yintercept = 4)) +
  geom_vline(aes(xintercept = 4))

ggplot(dtmplot) +
  geom_histogram(aes(x=floralcompetition_in_TAK))
ggplot(dtmplot) +
  geom_histogram(aes(x=nestcompetition_in_TAK))
ggplot(dtmplot) +
  geom_histogram(aes(x=ppp_in_TAK))


# return some of the papers with the highest mentions of stream terms --------------
DB1 <- DB1 %>% add_ref()

# floral competition, ok but not ideal
DB1 %>% 
  select(ref) %>% 
  filter(dtmplot$floralcompetition_in_TAK > 5) %>% 
  print(n=nrow(.))

# nest competition, most about nesting, but not always comparisons of native AND wild bees
DB1 %>% 
  select(ref) %>% 
  filter(dtmplot$nestcompetition_in_TAK > 3) %>% 
  print(n=nrow(.))

# ppp, better in getting ppp, but not always comparisons of native AND wild bees
DB1 %>% 
  select(ref) %>% 
  filter(dtmplot$ppp_in_TAK > 9) %>% 
  print(n=nrow(.))

# How about the small groups themselves, how good are they? ------------------------
DB1 %>% 
  select(ref) %>% 
  filter(dtmplot$cl_4 == dtmplot$cl_4[[1]]) %>%  # select any number between 1:2400
  print(n=nrow(.))

coms_summaries %>% 
  filter(group_id == dtmplot$cl_4[[1]]) %>%  # select any number between 1:2400
  select(group_top_keywords) %>% 
  unnest(cols = c(group_top_keywords)) %>% 
  print(n=nrow(.))

