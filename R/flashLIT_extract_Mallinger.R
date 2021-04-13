# flashLIT extract Mallinger

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# Mallinger ==================================================================

list.files("data/KeyReferences/Mallinger_2017")

# Table 1. Studies published from 1900±2016 examining potential competitive effects of managed bees on wild bees.
MT1 <- read_csv("data/KeyReferences/Mallinger_2017/Mallinger_2017_table1.csv")
# Table 2. Studies published from 1900±2016 examining the potential effect of managed bees on wild bees through changes in plant communities, including the spread of exotic plants.
MT2 <- read_csv("data/KeyReferences/Mallinger_2017/Mallinger_2017_table2.csv")
# Table 3. Studies published from 1900±2016 examining the potential transmission of pathogens from managed to wild bees.
MT3 <- read_csv("data/KeyReferences/Mallinger_2017/Mallinger_2017_table3.csv")
# Citations
MS2 <- read_csv("data/KeyReferences/Mallinger_2017/Mallinger_2017_S2_citations.csv")

# MT1 could be split into floral or nest or undetermined - this can be done by the 6, 7 and 8 column terms (measured variables)
t1 <- MT1[,c(6:8)] %>% pull() %>% paste(collapse = ", ") %>% 
  str_remove_all("\\([^)]*\\)") %>% str_split(", ") %>% .[[1]] %>% 
  trimws() %>% unique()

# we could categorise these into:

# presence (inc distribution, occupancy) or abundance (of a species) or relative abundance/densities
# composition & diversity (of a community)
# fitness (inc number of nests, cells, body sizes, colony mass, fecundity, survival, biomass)
# foraging (pollen & honey, visitaiton rate, foraging behaviour) & (pollen/nectar availability, concentration, tongue length, visitaiton rates)
# the displacement and niche overlap ones should be checked, but likely also floral. 

# However short of time for now so keep them in the original categories

# Assign these to a PID if possible -------------------------------------------------

# within the tables, these have an Author .. Year format e.g "Dick 2001",  "Aizen & Feinsinger 1994", "Abe et al. 2010","Kaiser- Bunbury & Müller 2009", "Sun et al. 2013a" 

# lets create a list with just these and the groups
Mallinger <- tibble(
  Reference = c(MT1$Reference, MT2$Reference, MT3$Reference),
  Group = c(rep("Competition_FN", length(MT1$Reference)),
            rep("Plants", length(MT2$Reference)),
            rep("Pathogen", length(MT3$Reference))
            )
)

# citations have the all authors, year, Title, full journal, volume, pages (no DOI) # note MS2 %>% filter(grepl("Kaiser", Citation)) ## is "Kaiser-Bunbury"

# to match abbreviated with full, we:
# Get first author name and year: take the 1st word, and the last word (words defined by a space or hyphen) 
Mallinger <- Mallinger %>% 
  mutate(Author1 = str_extract(Reference, "^[^ -]+")) %>% 
  mutate(Year = str_extract(Reference, "[a-z0-9]+$")) %>% 
  select(Author1, Year, Group)

# then we search the citations for strings that include both of these, and start with the first author
Mallinger <- Mallinger %>% 
  mutate(fullref = map2(Author1, Year, function(a,y,cv){
    cv[grepl(paste0("^", a), cv) & grepl(y, cv)]
  },
  cv = MS2$Citation))%>% 
  mutate(nm = map_dbl(fullref, length))

# there are a number of entries with typos (e.g. Abe 2011/2010 and Connor/Conner) or multiple entries in the list due to failues to distinguish with a letter
# the 0 group needs to be manually found, the doubles group first take distinct rows, as they come from the same groups.
Mallinger0 <- Mallinger %>% filter(nm == 0) %>% select(-fullref, -nm)
Mallinger1 <- Mallinger %>% filter(nm == 1) %>% select(-nm) %>% unnest(fullref)
Mallinger2 <- Mallinger %>% filter(nm == 2) %>% select(-nm) %>% unnest(fullref) %>% distinct()

Mallinger0
# check which is correct and correct the search, then we can match them
MS2$Citation[grepl("^Abe", MS2$Citation)]
Mallinger0[1,2] <- "2011"
MS2$Citation[grepl("^Conn", MS2$Citation)]
Mallinger0[2,1] <- "Conner"     
MS2$Citation[grepl("^Nishikawa", MS2$Citation)]
Mallinger0[3,2] <- "2016"
MS2$Citation[grepl("^F", MS2$Citation)]
Mallinger0[4,1] <- "Fürst"
MS2$Citation[grepl("^Otter", MS2$Citation)]
Mallinger0[5,1] <- "Otterstatter"

Mallinger0 <- Mallinger0 %>% 
  mutate(fullref = map2(Author1, Year, function(a,y,cv){
    cv[grepl(paste0("^", a), cv) & grepl(y, cv)]
  },
  cv = MS2$Citation))%>% 
  mutate(nm = map_dbl(fullref, length))

Mallinger0
Mallinger0 <- Mallinger0 %>% select(-nm) %>% unnest(fullref)

Mallinger <- bind_rows(Mallinger0, Mallinger1, Mallinger2) %>%
  arrange(Author1, Year)

# now we can also pull out the title (or at least some words from it), and clean the letters from the year
# and to make it easier for matching on title, lets take only the first 5 words, and convert to lowercase, remove punctuation
Mallinger <- Mallinger %>% 
  mutate(Title = map2_chr(Year, fullref, function(y,x){
    str_extract(x, paste0(y,"[a-z.]+[^.]+")) %>% 
      str_remove("^[^A-Z]+")
  })) %>% 
  mutate(Year = Year %>% str_remove("[a-z]") %>% as.integer()) %>% 
  mutate(title_5w = str_extract(Title, "[^ ]+([ ]+[^ ]+){4}") %>% tolower() %>% tm::removePunctuation()) %>% 
  mutate(Author1 = tolower(Author1))

# and use these to match it to our PID
 
test <- DB1 %>% 
  select(PID, author, Year = year, title_orig) %>% 
  mutate(Author1 = str_extract(author, "^[^ -,]+") %>% tolower()) %>% 
  mutate(title_5w = str_extract(title_orig, "[^ ]+([ ]+[^ ]+){4}") %>% tolower() %>% tm::removePunctuation())

joined <- left_join(Mallinger, test)

isjoined <- joined %>% filter(!is.na(PID)) #n=96
notjoined <- joined %>% filter(is.na(PID)) #n=57

# how many of the not joined can be joined on just the title?? 
tj <- notjoined %>% 
  select(-PID, -author, -title_orig) %>% 
  left_join(test, by = c("title_5w" = "title_5w")) 

tjj <- tj %>% 
  filter(!is.na(PID))

tjj[5,] %>% t()

tjj <- tjj[-5, ] %>% 
  select(Author1 = Author1.x, Year = Year.x, Group, fullref, Title, title_5w, PID, author, title_orig)

isjoined <- isjoined %>% bind_rows(tjj) #104

notjoined <- notjoined %>% filter(!fullref %in% isjoined$fullref ) #n=49

# which of these can we join by author-year? 
# notjoined$Author1 in test$Author1, and Year within 1 of the year 

auj <- notjoined %>% 
  select(-PID, -author, -title_orig) %>% 
  mutate(match_au_yrish = map2(Author1, Year, function(a,y){
    test %>% 
      filter(grepl(a, Author1)) %>% 
      filter(abs(y - Year)<2)
  })) %>% 
  mutate(nm = map_dbl(match_au_yrish, nrow)) %>% 
  filter(nm > 0)
  
auj1 <- auj %>% filter(nm == 1) %>% unnest(match_au_yrish, names_repair = 'unique')
auj1 %>% select(Title, title_orig) ## none of these are matches

auj2 <- auj %>% filter(nm > 1) %>% unnest(match_au_yrish, names_repair = 'unique')
auj2 %>% select(Title, title_orig)

aujj <- auj2[15,] %>% # is the only match, title has extra spaces
 select(Author1 = Author1...1, Year = Year...2, Group, fullref, Title, title_5w = title_5w...6, PID, author, title_orig)

isjoined <- isjoined %>% bind_rows(aujj) #105

notjoined <- notjoined %>% filter(!fullref %in% isjoined$fullref ) # n=48