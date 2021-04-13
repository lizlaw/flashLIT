#revtools example - shiny

library(revtools)
library(tidyverse)

file_names <- list.files("./data/raw_WoS_20201105/as_bib/", full = TRUE)
data <- revtools::read_bibliography(file_names)

data_unique <- screen_duplicates(data)

data_ti_screened <- screen_titles(data %>% sample_n(10))

data_ab_screened <- screen_abstracts(data %>% sample_n(5))

data_topic_screened <- screen_topics(data %>% sample_n(100))
