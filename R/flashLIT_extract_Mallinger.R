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


