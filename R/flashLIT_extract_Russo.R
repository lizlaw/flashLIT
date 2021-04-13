# flashLIT extract Russo

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# Russo ==================================================================
list.files("data/KeyReferences/Russo_2016")
# table 2, table S2

# Table 2  (non APis/bombus)
RussoNest1 <- c(15, 19, 22, 79, 80, 83, 87, 88)  # Nest sites
RussoFloral1 <- c(14, 22, 52, 56, 73, 76, 78, 84)  # Floral resources
RussoPath1 <- c(19, 60, 63, 73, 82, 86, 96, 97)  # Pathogens
RussoPlant1 <- c(19, 28, 39, 52, 73, 77, 78, 81, 84, 89, 90, 91, 92, 93, 94) # Plants
RussoPollNet1 <- c(22, 26, 90) # Networks
RussoPollin1 <- c(22, 62, 95) # change pollination
RussoEmpirical1 <- c(56, 63, 80, 81, 82, 84, 86, 91, 92, 93, 94, 96)

# Table S2 (APis, bombus)
RussoNest2 <- c(1,6,7,60, 62, 74:76)
RussoFloral2 <- c(2, 5, 6, 8:19, 20:22, 51, 52:54, 60, 62, 63, 67, 68, 70,74,75, 77:84, 85, 86)
RussoPath2 <- c(1:3, 5, 9, 23:26, 61, 26,64:66, 69,70, 76, 87:90)
RussoPlant2 <- c(6, 27:35, 58,59, 69,71,72,91:95, 96)
RussoPollNet2 <- c(36:38, 39:43, 44, 55:57, 68,81,97,98)
RussoIntrogression <- c(4, 62,76,88, 99,100)
RussoPollin2 <- c(37, 38,39, 45,46,47,48,49,50, 73,76,49,86,97,101:103)
RussoEmpirical2 <- c(8,9,10,11,12,13,14,15,16,17,18,19,23,24,25,26,27,28,29,30,31,32,33,34,35,37,38,39,40,41,42,43,44,45,48,49,50,
                     51,52,53,54,55,56,57,
                     58,59,
                     60,61,
                     28,58,59,68,69,71,72,73,
                     59,
                     4,25,26,29,38,39,47,49,58,67,68,70,72,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,97,98,99,100,101,102,103) %>% 
  sort() %>% unique()

# convert numbers to references
RussoRef1 <- read_csv("data/KeyReferences/Russo_2016/Russo_2016_references.csv")
RussoRef2 <- read_csv("data/KeyReferences/Russo_2016/Russo_2016_SIreferences.csv")

# match these to our dataset
RussoRef1 <- RussoRef1 %>% 
  mutate(
    Nest = RefID %in% RussoNest1,
    Floral = RefID %in% RussoFloral1,
    Path = RefID %in% RussoPath1,
    Plant = RefID %in% RussoPlant1,
    PollNet = RefID %in% RussoPollNet1,
    Pollin = RefID %in% RussoPollin1,
    Empirical = RefID %in% RussoEmpirical1
  )

RussoRef2 <- RussoRef2 %>% 
  mutate(
    Nest = RefID %in% RussoNest2,
    Floral = RefID %in% RussoFloral2,
    Path = RefID %in% RussoPath2,
    Plant = RefID %in% RussoPlant2,
    PollNet = RefID %in% RussoPollNet2,
    Introg = RefID %in% RussoIntrogression,
    Pollin = RefID %in% RussoPollin2,
    Empirical = RefID %in% RussoEmpirical2
  )

# Match these with DB1
# we can match the first author and the title_orig words from DB1 to the Citation
