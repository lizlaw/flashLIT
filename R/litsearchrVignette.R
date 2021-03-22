# litsearchr vignette

library(litsearchr)

search_directory <- system.file("extdata", package="litsearchr")

naiveimport <-
  litsearchr::import_results(directory = search_directory, verbose = TRUE)
naiveresults <-
  litsearchr::remove_duplicates(naiveimport, field = "title", method = "string_osa")

rakedkeywords <-
  litsearchr::extract_terms(
    text = paste(naiveresults$title, naiveresults$abstract),
    method = "fakerake",
    min_freq = 2,
    ngrams = TRUE,
    min_n = 2,
    language = "English"
  )

taggedkeywords <-
  litsearchr::extract_terms(
    keywords = naiveresults$keywords,
    method = "tagged",
    min_freq = 2,
    ngrams = TRUE,
    min_n = 2,
    language = "English"
  )

all_keywords <- unique(append(taggedkeywords, rakedkeywords))
naivedfm <-
  litsearchr::create_dfm(
    elements = paste(naiveresults$title, naiveresults$abstract),
    features = all_keywords
  )
naivegraph <-
  litsearchr::create_network(
    search_dfm = naivedfm,
    min_studies = 2,
    min_occ = 2
  )

cutoff <-
  litsearchr::find_cutoff(
    naivegraph,
    method = "cumulative",
    percent = .80,
    imp_method = "strength"
  )
reducedgraph <-
  litsearchr::reduce_graph(naivegraph, cutoff_strength = cutoff[1])
searchterms <- litsearchr::get_keywords(reducedgraph)
head(searchterms, 20)

mysearchterms <-
  list(
    c(
      "picoides arcticus",
      "black-backed woodpecker",
      "cavity-nesting birds",
      "picoides tridactylus",
      "three-toed woodpecker"),
    c(
      "wildfire",
      "burned forest",
      "post-fire",
      "postfire salvage logging",
      "fire severity",
      "recently burned"
    )
  )
my_search <-
  litsearchr::write_search(
    groupdata = mysearchterms,
    languages = "English",
    stemming = TRUE,
    closure = "none",
    exactphrase = TRUE,
    writesearch = FALSE,
    verbose = TRUE
  )

cat(my_search)


## from 156 input
## Results: 103
# (from Web of Science Core Collection) Date = 19.11.2020
# You searched for: TOPIC: ((("picoid* arcticus*"  OR "black-back* woodpeck*"  OR "cavity-nest* bird*"  OR "picoid* tridactylus*"  OR "three-to* woodpeck*")  AND (wildfir*  OR "burn* forest*"  OR post-fir*  OR "postfir* salvag* logging"  OR "fire* sever*"  OR "recent* burn*")))
# Timespan: All years. Indexes: SCI-EXPANDED, SSCI, A&HCI, ESCI.
