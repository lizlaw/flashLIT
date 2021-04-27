# flashLIT alternative clustering methods

# ====================================

# convert tidy back to distance object
cast_dist <- function(df, labs, distMethod = "euclidean"){
  x <- df
  structure(temp, Size = length(labs), Labels = labs,
                   Diag = FALSE, Upper = FALSE, method = distMethod, 
                   class = "dist")
}

# visualise what custering we might see using factoextra

fviz_pca_ind(prcomp(df), title = "PCA - Iris data", geom = "point", ggtheme = theme_classic())
km.res1 <- kmeans(df, 3)
fviz_cluster(list(data = df, cluster = km.res1$cluster), ellipse.type = "norm", geom = "point", stand = FALSE, palette = "jco", ggtheme = theme_classic())
fviz_dend(hclust(dist(df)), k = 3, k_colors = "jco", as.ggplot = TRUE, show_labels = FALSE)
get_clust_tendency(df, n = nrow(df)-1, graph = FALSE)


library("factoextra")
res.pca <- PCA(df,  graph = FALSE)
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_pca_ind(res.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca, repel = TRUE)

## plotting known groups using habillage
iris.pca <- PCA(iris[,-5], graph = FALSE) 
fviz_pca_ind(iris.pca,
             label = "none", # hide individual labels
             habillage = iris$Species, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE # Concentration ellipses
)

# plotting results of clustering
km.res <- kmeans(scale(USArrests), 4, nstart = 25)
fviz_cluster(km.res, data = df,
             palette = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"),
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot"
)

res <- hcut(USArrests, k = 4, stand = TRUE)
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

my_data <- scale(USArrests)
fviz_nbclust(my_data, kmeans, method = "gap_stat") # optimal k clusters

# ====================================

# Hierarchical clustering with dynamic cutting

# hierarchical clustering - agglomerative with complete/maximal or ward's linkage method. (here based on euclidean distance, but bray-curtis might be better). This method is simple, and should give reasonably good clusters at the fine level (but not at the coarse splits). 
# dynamic cuts can increase ability to derive good clusters, including nested clusters, and reduce impact of outliers (which can be removed)
# this method is likely to be sensitive to different data. 

d <- stats::dist(df, method = "euclidean") 
d <- ecodist::distance(df, method = "bray")  # for abundances, use sørensen for presence-absence
hc1 <- fastcluster::hclust(d, method = "complete")
hc2 <- cluster::agnes(d, method = "complete")
hc2$ac
# the agglomerative coefficient (AC), measures the amount of clustering structure found. Values closer to 0 suggest less well-formed clusters. However, the AC tends to become larger as  n increases, so it should not be used to compare across data sets of very different sizes.

# test multiple methods
m <- c("complete", "ward")
names(m) <- c( "complete", "ward")
ac <- function(x) {
  agnes(ames_scale, method = x)$ac
}
purrr::map_dbl(m, ac)

# Number of members in each hard-cut cluster
sub_grp <- cutree(hc2, k = 8)
table(sub_grp)

# dynamic cutting --- 

# compare 
dendextend::cutree(hc2, h=0.9)
WGCNA::cutreeStatic(hc2, cutHeight = 0.9, minSize = 10)

dynamicTreeCut::cutreeDynamic(hc2, cutHeight = NULL, minClusterSize = 10,
                              method = "hybrid", distM = d, deepSplit = 4,
                              pamStage = TRUE, pamRespectsDendro = TRUE)  ## many other advanced options

# ==================================================
# (simple) fuzzy clustering
#  the simplest being based on K-means (cluster::fanny, or e1071::cmeans, which includes an option to specify initial means). m=2 is recommended

cluster::fanny(d, k, memb.exp = 2)
cm <- e1071::cmeans(df, centres, iter.max, dist = "euclidean", m=2)
library(corrplot)
corrplot(cm$membership, is.corr = FALSE)

library(factoextra)
fviz_cluster(list(data = df, cluster=cm$cluster), 
             ellipse.type = "norm",
             ellipse.level = 0.68,
             palette = "jco",
             ggtheme = theme_minimal())


# ===============================================
# Model based clustering (more advanced fuzzy clustering)
# https://bradleyboehmke.github.io/HOML/model-clustering.html
#The key idea behind model-based clustering is that the data are considered as coming from a mixture of underlying probability distributions. The most popular approach is the Gaussian mixture model (GMM) (Banfield and Raftery 1993) which assume the clusters can be created using  k Gaussian distributions. Can parameterise to have many different options of the shape of the clusters. 

library(mclust)
data(geyser, package = 'MASS')
url <- "https://koalaverse.github.io/homlr/data/my_basket.csv"
my_basket <- readr::read_csv(url)
geyser_mc <- Mclust(geyser, G = 3)

# Plot results
plot(geyser_mc, what = "density")
plot(geyser_mc, what = "uncertainty")

# test multiple methods:
geyser_optimal_mc <- Mclust(geyser, 1:20)
summary(geyser_optimal_mc)
legend_args <- list(x = "bottomright", ncol = 5)
plot(geyser_optimal_mc, what = 'BIC', legendArgs = legend_args)
plot(geyser_optimal_mc, what = 'classification')
plot(geyser_optimal_mc, what = 'uncertainty')

# =================================================
# MuNCut (R package NcutYX) 
# graph-based simulated annealing optimization method which uses multiple layers of different information


# ======================================================
# Comparing clusters 
# https://www.datanovia.com/en/lessons/comparing-cluster-dendrograms-in-r/
# https://stats.stackexchange.com/questions/63546/comparing-hierarchical-clustering-dendrograms-obtained-by-different-distances
# http://talgalili.github.io/dendextend/reference/cor.dendlist.html

# cophenetic ---------
# From cophenetic: The cophenetic distance between two observations that have been clustered is defined to be the intergroup dissimilarity at which the two observations are first combined into a single cluster. Note that this distance has many ties and restrictions.cor_cophenetic calculates the correlation between two cophenetic distance matrices of the two trees.The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar. For exact p-value one should result to a permutation test. One such option will be to permute over the labels of one tree many times, and calculating the distriubtion under the null hypothesis (keeping the trees topologies constant).Notice that this measure IS affected by the height of a branch.

cor_cophenetic(dend1, dend2)
cor.dendlist(dend_list, method = "cophenetic")

# Bakers -----
# Calculate Baker's Gamma correlation coefficient for two trees (also known as Goodman-Kruskal-gamma index). Baker's Gamma (see reference) is a measure of accosiation (similarity) between two trees of heirarchical clustering (dendrograms).It is calculated by taking two items, and see what is the heighst possible level of k (number of cluster groups created when cutting the tree) for which the two item still belongs to the same tree. That k is returned, and the same is done for these two items for the second tree. There are n over 2 combinations of such pairs of items from the items in the tree, and all of these numbers are calculated for each of the two trees. Then, these two sets of numbers (a set for the items in each tree) are paired according to the pairs of items compared, and a spearman correlation is calculated. The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar. For exact p-value one should result to a permutation test. One such option will be to permute over the labels of one tree many times, and calculating the distriubtion under the null hypothesis (keeping the trees topologies constant).Notice that this measure is not affected by the height of a branch but only of its relative position compared with other branches

cor_bakers_gamma(dend1, dend2)
cor.dendlist(dend_list, method = "baker")

dend1 <- df %>% dist %>% hclust("complete") %>% as.dendrogram
dend2 <- df %>% dist %>% hclust("single") %>% as.dendrogram
dend3 <- df %>% dist %>% hclust("average") %>% as.dendrogram
dend4 <- df %>% dist %>% hclust("centroid") %>% as.dendrogram


# plot corrplot
dend_list <- dendlist("Complete" = dend1, "Single" = dend2,
                      "Average" = dend3, "Centroid" = dend4)
cors <- cor.dendlist(dend_list)

library(corrplot)
corrplot(cors, "pie", "lower")

# entanglement plots
dendlist(dend1, dend2) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram()  
dendlist(dend1, dend2) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  entanglement()


# Modularity ----
# Graph methods typically are evaluated with Newman and Girvan’s Modularity method

# Silhouette / agglomerative coefficient
# agglomerative coefficient (ac) is similar to the silhouette coefficient (sc). It measures the dissimilarity of an object to the first cluster it joins, divided by the dissimilarity of the final merger in the cluster analysis, averaged across all samples. Low values reflect tight clustering of objects, larger values indicate less well-formed clusters. The agglomerative coefficient increases with sample sizes, making comparisons among data sets difficult, but in our case we have the same groups for the comparisons we wish to make using this method.

wcke<-eclust(pizs, "kmeans", hc_metric="euclidean",k=3)
fviz_cluster(wcke, geom = "point", ellipse.type = "norm", ggtheme = theme_minimal())
sile<-silhouette(wcke$cluster, dist(pizs))
fviz_silhouette(sile)
