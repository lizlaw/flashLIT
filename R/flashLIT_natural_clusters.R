# flashLIT clustering documents (natural clusters)

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

# clustering the data into natural clusters 

# data going in could consist of 
# title, abstract, keywords, cited_references (i.e. reference lists), internal citations (known citation matches within our data), and authors

# There are many different types of clustering algorithms we could use. 

# Scimeetr uses a graph clustering method (louvain or fast_greedy). It takes a graph where the nodes are papers, and the arcs are shared links (cited_references, title, abstract, keywords, authors). Titles and abstracts are composited into dtm before coupling. This coupling can weight by the type of the link, as well as the strength of the link.

# revtools uses an LDA topic model approach (from a dtm), also including bigrams where appropriate.
# this is slow when working with multiple data types. 
# also would need to work out how to combine the different inputs, or outputs. (no way of weighting them)

# the benefits of louvain are that it determines when to stop automatically. However, there is no garauntee the groups will be the same when updating inputs.
# the LDA topic models need a number of topics input. Similarly, there is no gauranteee the groups will be the same when updating inputs.

# therefore, for the categorisation clustering, we need a more prescriptive method. 
# this can draw from knowledge of the natural clusters though, so we first cluster into natural clusters.

# --- modified scimeetr process ---
# Scimeetr uses https://jangorecki.gitlab.io/-/data.table/-/jobs/640724/artifacts/public/html/data.table.html 
# this creates a graph weighted by the sum similarity across the different components for the pairs of documents.
# this is faster and more flexible clustering (once the graph has been constructed) than offered by the topic models

# create the graph object - !! can be slow
# we want to up-weight the information coming from the title, abstract, and keywords, because these form the basis of categorising the papers
# the abstract might be more important, as titles amd selected keywords can really vary. 
# the other information can help discern the groups, but should not be overweighted.

graph_DB1 <- scimeetr_coupling(DB1,
                               coupling_by = "title, abstract, keywords, keywordsplus, bibliography, authors, cr_journal", 
                               w.tic = 2, w.abc = 3, w.kw = 2, w.kwp = 1, w.bic = 1, w.auc = 1, w.joc = 1)

# saveRDS(graph_DB1, "flashLIT_graph_DB1.RDS")
# graph_DB1 <- readRDS("flashLIT_graph_DB1.RDS")


# clustering once/twice would be:
# # apply the function using the default louvain method (https://igraph.org/r/doc/cluster_louvain.html) and return a vector of group membership
# cl1 <- clusterize(graph = coupled_DB) %>% igraph::membership()
# CL1 <- tibble(
#   PID = names(cl1),
#   cl1 = as.character(cl1)
# )
# 
# # how many in each group?
# table(CL1$cl1)
# 
# #subset the graph according to these memberships, before splitting again
# CL2 <- map_dfr(unique(CL1$cl1), function(x){
#   vpids <- CL1 %>% filter(cl1 == x) %>% pull(PID)
#   sg <- igraph::induced_subgraph(graph = coupled_DB, vids = vpids, impl = "auto")
#   cl2 <- clusterize(graph = sg) %>% igraph::membership() 
#   CL2 <- tibble(
#     PID = names(cl2),
#     cl2 = paste(x, cl2, sep = "_"))
#   CL2
# })

# but we need to iterate this:

coms_DB1 <- iteratively_cluster(graph_DB1, 50)

# endscript