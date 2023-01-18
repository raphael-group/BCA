setwd("~/Desktop/Research/structured-priors/experiments/cell-cycle/")
source("utils.R")
# if(!exists("foo", mode="function")) source("util.R")
set.seed(42)

data <- readRDS("data/real-data/HumanPG-guo/germline-human-both_guo.rds")
#geo_dist_true <- dynwrap::calculate_geodesic_distances(data2, data2$waypoint_cells)
#heatmap(geo_dist_true, Colv = NA, Rowv = NA, scale="none")
X <- log1p(data$expression)
# Relabel Cluster labels so that they are in order.; 4, 7, 11, 17, 19
#  10W  11W   4W   7W   8W F17W M19W 
#   41   37   22   55   37   33   47 
# Moreover, merge 10W with 11W as they are so few 10W cells and they are very close to 11W
# data$milestone_percentages$milestone_id[data$milestone_percentages$milestone_id == "4W"] <- "04W"
data$milestone_percentages$milestone_id[data$milestone_percentages$milestone_id == "4W"] <- "04W"
data$milestone_percentages$milestone_id[data$milestone_percentages$milestone_id == "7W"] <- "07W"
data$milestone_percentages$milestone_id[data$milestone_percentages$milestone_id == "8W"] <- "08W"
# data$milestone_percentages$milestone_id[data$milestone_percentages$milestone_id == "8W"] <- "08W"
data$milestone_percentages$milestone_id[data$milestone_percentages$milestone_id == "10W"] <- "11W"


all_clusts <- data$prior_information$groups_id$group_id

all_clusts[all_clusts == "4W"] <- "04W"
all_clusts[all_clusts == "7W"] <- "07W"
all_clusts[all_clusts == "8W"] <- "08W"
all_clusts[all_clusts == "10W"] <- "11W"

clusts <- all_clusts

# clusts[all_clusts == "07W"] <- "04W"
# clusts[all_clusts == "11W"] <- "04W"

#clusts[all_clusts == "07W"] <- "11W"

clust_order <- sort(unique(all_clusts))
int_clusts <- clusts_to_ints(all_clusts, sort(unique(all_clusts)))

# assign cells to closest milestone
fwrite(as.list(all_clusts), "./data/real-data/HumanPG-guo/clusters-PG.csv")

############################## PCA ###################################
X_svd <- svd(scale(X, center=TRUE, scale=TRUE))
Y_pca <- X_svd$u %*% diag(X_svd$d)

fwrite(as.matrix(Y_pca[,1:50]), "data/real-data/HumanPG-guo/PCA-50-PG.csv")

sds_pca <- slingshot(Y_pca[,1:50], clusterLabels = clusts) 
ps_pca <- sds_pca@assays@data@listData[["pseudotime"]]
taus <- all_branches_score(ps_pca, int_clusts)

plot_sling_curve(Y_pca[,1:50], clusts, sds_pca, "PC", "", clust_order, "Cell Type", 5, FALSE)
ggsave("plots/pca-hpgc-slingshot.pdf")

c_mat_pca <- confusion_mat(ps_pca, clusts_to_ints(clusts, clust_order))
cnf_mat_pca <- c_mat_pca[[1]]$table
cat("Classification Accuracy", c_mat_pca[[1]]$overall['Accuracy'], "\n")
pca_title <- "Confusion Matrix for PCA and Slingshot on PGC Data; Lineage 1"
plot_cnf_mat(cnf_mat_pca, clust_order, pca_title)

c_mat_pca <- confusion_mat(ps_pca, clusts_to_ints(clusts, clust_order))
cnf_mat_pca <- c_mat_pca[[2]]$table
cat("Classification Accuracy", c_mat_pca[[2]]$overall['Accuracy'], "\n")
pca_title <- "Confusion Matrix for PCA and Slingshot on PGC Data; Lineage 2"
plot_cnf_mat(cnf_mat_pca, clust_order, pca_title)

#plot_traj(pca_df, "ps", "PCA-sling-ord")
########################################################################
Y_our <- our_method(X, clusts)[,1:4]

sds_our <- slingshot(Y_our, clusterLabels = clusts, start.clus = "04W")
ps_our <- sds_our@assays@data@listData[["pseudotime"]]

taus <- all_branches_score(ps_our, int_clusts)

plot_sling_curve(Y_our, clusts, sds_our, "BCC", 
                 "MBCVE with Slingshot on Human Primordial Germ Cell data (n=272)",
                 "Cell Type", clust_order, 3, FALSE)

###########################################
Y_wour <- wour_method(X, clusts)
fwrite(Y_wour, "data/real-data/HumanPG-guo/BCA-PG.csv")
sds_wour <- slingshot(Y_wour, clusterLabels = all_clusts, start.clus = "04W")
ps_wour <- sds_wour@assays@data@listData[["pseudotime"]]
taus <- all_branches_score(ps_wour, int_clusts)

plot_sling_curve(Y_wour, all_clusts, sds_wour, "BCC", 
                 "MBCVE with Slingshot on Human PGC data (n=272)",
                 "Cell Type", clust_order, 3)

c_mat_wour <- confusion_mat(ps_wour, clusts_to_ints(all_clusts, clust_order))
cnf_mat_wour <- c_mat_wour[[1]]$table
cat("Classification Accuracy", c_mat_wour[[1]]$overall['Accuracy'], "\n")
mbcve_title <- "Confusion Matrix for MBCVE and Slingshot on PGC Data; Lineage 1"
plot_cnf_mat(cnf_mat_wour, clust_order, mbcve_title)

c_mat_wour <- confusion_mat(ps_wour, clusts_to_ints(all_clusts, clust_order))
cnf_mat_wour <- c_mat_wour[[2]]$table
cat("Classification Accuracy", c_mat_wour[[2]]$overall['Accuracy'], "\n")
mbcve_title <- "Confusion Matrix for MBCVE and Slingshot on PGC Data; Lineage 2"
plot_cnf_mat(cnf_mat_wour, clust_order, mbcve_title)

#Kendall's Tau for Lineage 1 is 0.922922174962639. 
#Kendall's Tau for Lineage 2 is 0.747261397179403. 


########################## UMAP #####################################
Y_umap <- umap(X)
fwrite(Y_umap$layout, "data/real-data/HumanPG-guo/UMAP-PG.csv")
sds_umap <- slingshot(Y_umap$layout, clusterLabels = clusts, start.clus = "04W")
ps_umap <- sds_umap@assays@data@listData[["pseudotime"]]
taus <- all_branches_score(ps_umap, int_clusts)

plot_sling_curve(Y_umap$layout, clusts, sds_umap, "UMAP", 
                 "UMAP with Slingshot on Human Primordial Germ Cell data (n=272)",
                 "Cell Type", 3)

############################ mLDA ####################################
Y_lda <- do.lda(Y_pca[,1:50], as.factor(clusts))
#fwrite(Y_lda$Y, "data/real-data/HumanPG-guo/LDA-PG.csv")
sds_lda <- slingshot(Y_lda$Y, clusterLabels = all_clusts, start.clus = "04W")
ps_lda <- sds_lda@assays@data@listData[["pseudotime"]]
taus <- all_branches_score(ps_lda, int_clusts)

plot_sling_curve(Y_lda$Y, all_clusts, sds_lda, "LD", 
                 "", clust_order,
                 "Cell Type", 5, FALSE)
ggsave("plots/lda-hpgc-slingshot.pdf")

c_mat_lda <- confusion_mat(ps_lda, clusts_to_ints(all_clusts, clust_order))
cnf_mat_lda <- c_mat_lda[[1]]$table
cat("Classification Accuracy", c_mat_lda[[1]]$overall['Accuracy'], "\n")
lda_title <- "Confusion Matrix for LDA and Slingshot on PGC Data; Lineage 1"
plot_cnf_mat(cnf_mat_lda, clust_order, lda_title)

c_mat_lda <- confusion_mat(ps_lda, clusts_to_ints(all_clusts, clust_order))
cnf_mat_lda <- c_mat_lda[[2]]$table
cat("Classification Accuracy", c_mat_lda[[2]]$overall['Accuracy'], "\n")
lda_title <- "Confusion Matrix for MBCVE and Slingshot on PGC Data; Lineage 2"
plot_cnf_mat(cnf_mat_lda, clust_order, lda_title)





