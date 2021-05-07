# flashLIT alternative clustering methods

# here we examine 
# 1) effect of using raw data, cleaned data, or stemmed data
# 2) effect of using TAK, TAKKP, Au, Bib, TAKKPAu or TAKKPBib
# 3) effect of using Euclidean or Bray-Curtis metrics for distance calculations
# 4) effect of combining data types using mean or euclidean distance
# ====================================

# Process input data-----------

# input data:
RCS <- tibble(
  tidyData = indata,
  label = names(indata)) %>% 
  separate(label, into = c("dataType", "cleanLevel", "distMethod"), remove = FALSE)

# we only want the smallest subgroup for this task
setPIDs <- DB1 %>% filter(MNI_n1_TAK == 1) %>% pull(PID)

# cast into a distance object - using just the specified set - and from there into a graph object
RCS <- RCS %>% 
  mutate(distData = map(tidyData, cast_dist, labs = setPIDs)) %>% 
  mutate(graph = map(distData, ~as.matrix(1/.x) %>% 
          igraph::graph_from_adjacency_matrix(., mode = "undirected", weighted = TRUE, diag = FALSE)))
  
# visualise the dissimilarity --- !! SLOW
jpeg("data/plot_distViz_rawCleanStem.jpeg", width = 20, height = 12, units = "in", res = 300)
    cowplot::plot_grid(plotlist = map(RCS$distData, fviz_dist, show_labels = FALSE), labels = RCS$label, nrow = 4)
dev.off()
# these suggest that the grouping may not be so distinct for these comparisons

# heirarchical clustering ----
# alternatives here are stats::hclust, or agnes(), we use fastcluster for speed
# tree cutting alternatives:
# HCct = set cut height, HCcts static cut height with minimum size, HCctd dynamic cuts without and with a PAM assignment stage
# note, dynamic hybrid cuts use a modified function as the package function has a bug which resulted in only the label information being returned.
# the modified version results in all the documented results being returned (labels, cores, etc.)
RCS <- RCS %>% 
  mutate(HCo = map(distData, fastcluster::hclust, method = "ward.D2")) %>% 
  mutate(HCct = map2(HCo, distData, ~dynamicTreeCut::cutreeDynamic(.x,
                              cutHeight = NULL, minClusterSize = 0,
                              deepSplit = FALSE, method = "tree"))) %>% 
  mutate(HCcts = map2(HCo, distData, ~dynamicTreeCut::cutreeDynamic(.x,
                              cutHeight = NULL, minClusterSize = 20,
                              deepSplit = TRUE, method = "tree"))) %>% 
  mutate(HCctd = map2(HCo, distData, ~cutreeHybrid_fullres(.x,
                      cutHeight = NULL, minClusterSize = 20,
                      distM = as.matrix(.y), deepSplit = 4,
                      pamStage = FALSE))) %>% 
  mutate(HCctdp = map2(HCo, distData, ~cutreeHybrid_fullres(.x,
                      cutHeight = NULL, minClusterSize = 20,
                      distM = as.matrix(.y), deepSplit = 4,
                      pamStage = TRUE, pamRespectsDendro = TRUE)))

# graph-based clustering ---
# further alternatives ncut from NCutYX does not yet take bray-curtis, and leiden need additional python modules
RCS <- RCS %>% 
  mutate(fastgreedy = map(graph, igraph::cluster_fast_greedy)) %>% 
  mutate(lourvain = map(graph, igraph::cluster_louvain)) 
 
# Comparing clusters ---

# number of clusters and median cluster size
RCS <- RCS %>% 
  mutate(n_HCct = map_dbl(HCct, n_distinct),
         n_HCcts = map_dbl(HCcts, n_distinct),
         n_HCctd = map_dbl(HCctd, ~.x$labels %>% n_distinct()),
         n_HCctdp = map_dbl(HCctdp, ~.x$labels %>% n_distinct())
  ) %>% 
  mutate(gm_HCct = map_dbl(HCct, median_groupsize),
         gm_HCcts = map_dbl(HCcts, median_groupsize),
         gm_HCctd = map_dbl(HCctd, ~.x$labels %>% median_groupsize()),
         gm_HCctdp = map_dbl(HCctdp, ~.x$labels %>% median_groupsize())
  )

RCS %>% 
  select(label:distMethod, n_HCct:gm_HCctdp)

# cluster quality 
# https://cran.r-project.org/web/packages/clusterCrit/vignettes/clusterCrit.pdf
# intCriteria(traj, part, crit) # internal quality criteria --- see clusterCrit_indices.xlsx 
# extCriteria(part1, part2, crit) # external quality criteria
# bestCriterion(x, crit) # best partitions according to a criteria
# concordance(part1, part2) # concordance matrix between two partitions

# intCriteria
# choose a set of common/modern maximum (3) and minimum (2) based indices
mycrit = c("Silhouette", "Xie",  "PBM", "Davies", "GDI")
# clusterCrit::intCriteria(traj = as.matrix(RCS$distData[[1]]), part = as.integer(RCS$HCct[[1]]), crit = mycrit)

# !! SLOW
# for HCctd and HCctdp need $labels
# can't evaluate the HCcts and HCctd this way because they have the 0 cluster - it is not equivalent
RCS <- RCS %>% 
  mutate(crit_HCct = map2(HCct, distData, ~clusterCrit::intCriteria(traj = as.matrix(.y), part = as.integer(.x), crit = mycrit)),
         #crit_HCcts = map2( HCcts, distData, ~clusterCrit::intCriteria(traj = as.matrix(.y), part = as.integer(.x), crit = mycrit)),
         #crit_HCctd = map2(HCctd, distData, ~clusterCrit::intCriteria(traj = as.matrix(.y), part = as.integer(.x$labels), crit = mycrit)),
         crit_HCctdp = map2(HCctdp, distData, ~clusterCrit::intCriteria(traj = as.matrix(.y), part = as.integer(.x$labels), crit = mycrit))
  )

# summaries of results  ===============================================================
# group size and number of groups  ----
res1 <- RCS %>% 
  select(label:distMethod, n_HCct:gm_HCctdp) %>% 
  pivot_longer(n_HCct:gm_HCctdp) %>% 
  separate(name, into = c("score", "algorithm")) %>% 
  pivot_wider(names_from = score, values_from = value) %>% 
  mutate(cleanLevel = factor(cleanLevel, levels = c("raw", "clean", "stem"))) 

# cluster values ---
clusterCrit::bestCriterion(c(3,5), "davies") # lower is better
clusterCrit::bestCriterion(c(0.1,0.3), "gdi") # higher is better
clusterCrit::bestCriterion(c(50,150), "pbm") # higher is better
clusterCrit::bestCriterion(c(-0.2,0), "silh") # higher is better
clusterCrit::bestCriterion(c(1,20), "xie") # lower is better

res2 <- RCS %>% 
  select(label:distMethod, crit_HCct:crit_HCctdp) %>% 
  pivot_longer(crit_HCct:crit_HCctdp) %>% 
  separate(name, into = c("score", "algorithm")) %>% 
  pivot_wider(names_from = score, values_from = value) %>% 
  unnest_longer(crit) %>% 
  mutate(cleanLevel = factor(cleanLevel, levels = c("raw", "clean", "stem"))) %>% 
  mutate(crit = map2_dbl(crit_id, crit, ~ifelse(.x %in% c("davies_bouldin", "xie_beni"), -.y, .y))) %>% 
  mutate(crit_id = map_chr(crit_id,  ~ifelse(.x %in% c("davies_bouldin", "xie_beni"), paste0("-",.x), .x))) 
            
# effect of cleaning ===============================================================
ggplot(res1 %>% dplyr::filter(dataType == 'TAK')) + 
  geom_path(aes(x=n, y = gm, color = cleanLevel, group = interaction(dataType, distMethod, algorithm)))+
  geom_point(aes(x=n, y = gm, colour = cleanLevel, shape = distMethod), size = 4, alpha = 0.8) +
  #ggrepel::geom_text_repel(aes(x=n, y = gm, colour = cleanLevel, label = cleanLevel), max.overlaps = 20) +
  labs(title = "Effect of cleaning (TAK only)", x = "number of clusters", y = "median cluster size") +
  facet_grid(.~algorithm, scales = "free_x") +
  scale_color_viridis_d()

ggplot(res2 %>% dplyr::filter(dataType == 'TAK')) + 
  geom_line(aes(x=cleanLevel, y = crit, color = cleanLevel, group = interaction(dataType, distMethod)))+
  geom_point(aes(x=cleanLevel, y = crit, colour = cleanLevel, shape = distMethod), size = 4, alpha = 0.8) +
  labs(title = "Effect of cleaning (TAK only)", subtitle = "(higher score is better)", x = "cluster method", y = "Score") +
  facet_grid(crit_id~algorithm, scales = "free_y") +
  scale_color_viridis_d()

# effect of data source ===============================================================

ggplot(res1 %>% dplyr::filter(cleanLevel == 'stem' | (dataType %in% c("Au", "Bib")))) + 
  geom_point(aes(x=n, y = gm, colour = dataType, shape = distMethod), size = 4) +
  labs(title = "Effect of data source (stem only)", x = "number of clusters", y = "median cluster size") +
  facet_grid(.~algorithm, scales = "free_x") +
  scale_color_viridis_d()

ggplot(res2 %>% dplyr::filter(cleanLevel == 'stem' | (dataType %in% c("Au", "Bib")))) + 
  geom_line(aes(x=dataType, y = crit, color = dataType, group = distMethod))+
  geom_point(aes(x=dataType, y = crit, colour = dataType, shape = distMethod), size = 4) +
  labs(title = "Effect of data source (stem only)", subtitle = "(higher score is better)", x = "cluster method", y = "Score") +
  facet_grid(crit_id~algorithm, scales = "free_y") +
  scale_color_viridis_d()