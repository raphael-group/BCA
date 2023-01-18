setwd("~/Desktop/Research/structured-priors/experiments/cell-cycle/")
source("utils.R")
library(cowplot) # for putting multiple plots on grid
# if(!exists("foo", mode="function")) source("util.R")


data <- readRDS("data/real-data/HSC-kowalczyk/aging-hsc-old_kowalczyk.rds")
all_clusts <- data$prior_information$groups_id$group_id
#good_inds <- clusts == clusts
X <- data$expression
t <- data$progressions$percentage
clusts <- all_clusts
unique(all_clusts)

# clusts[all_clusts == "ST-HSC"] <- "MPP"
# clusts[all_clusts == "ST-HSC"] <- "LT-HSC"

# num_mid_c <- floor(unname(table(all_clusts)[3]) / 2)
# mid_c_inds1 <- sort(sample(which(all_clusts == "ST-HSC"), num_mid_c, replace=FALSE))
# mid_c_inds2 <- setdiff(which(all_clusts == "ST-HSC"), mid_c_inds1)
# clusts[mid_c_inds1] <- "MPP"
# clusts[mid_c_inds2] <- "LT-HSC"



# geo_dist_true <- dynwrap::calculate_geodesic_distances(data, data$waypoint_cells)
# #ggplot(geo_dist_true) + geom_tile()
# heatmap(geo_dist_true, Colv = NA, Rowv = NA, scale="none")
# 
# percents <- data[["milestone_percentages"]][["percentage"]]

clust_order <- rev(c("MPP", "ST-HSC", "LT-HSC"))
int_clusts <- clusts_to_ints(all_clusts, clust_order)

#fwrite(as.list(all_clusts), "data/real-data/HSC-kowalczyk/clusters-HSC.csv")

############################## PCA ###################################
X_svd <- svd(scale(X, center=TRUE, scale=FALSE))
Y_pca <- X_svd$u %*% diag(X_svd$d)

fwrite(as.matrix(Y_pca[,1:50]), "data/real-data/HSC-kowalczyk/PCA-50-HSC.csv")

sds_pca <- slingshot(Y_pca[,1:50], clusts, start.clus = "LT-HSC")
ps_pca <- sds_pca@assays@data@listData[["pseudotime"]]


order_pca_structdr <- fread("data/real-data/HSC-kowalczyk/PCA-StructDR-order.csv")$V1
#names(order_pca_structdr) <- "Lineage1"
taus <- simple_score(order_pca_structdr, int_clusts)

taus <- all_branches_score(ps_pca, int_clusts)

(p_pca <- plot_sling_curve(Y_pca, clusts, sds_pca, "PC", 
                 "PCA with Slingshot on Mouse HSCs (n=873)", clust_order,
                 "Cell Type", 3, FALSE))

c_mat_pca <- confusion_mat(ps_pca, clusts_to_ints(clusts, clust_order))
cnf_mat_pca <- c_mat_pca[[1]]$table
cat("Classification Accuracy", c_mat_pca[[1]]$overall['Accuracy'], "\n")
pca_title <- "Confusion Matrix for PCA and Slingshot on HSC Data"
plot_cnf_mat(cnf_mat_pca, clust_order, pca_title)

########################################################################
Y_our <- as.matrix(our_method(X, clusts)[,1])

sds_our <- slingshot(Y_our, clusts, start.clus = "LT-HSC")
ps_our <- sds_our@assays@data@listData[["pseudotime"]]
taus <- all_branches_score(ps_our, int_clusts)




plot_sling_curve(Y_our, clusts, sds_our, "BCC", 
                 " with Slingshot on Mouse Hematopoetic Stem Cells (n=873)","Cell Type", clust_order,
                 3, TRUE)

###########################################
Y_wour <- wour_method(X, clusts)[,1:2]
#fwrite(Y_wour, "data/real-data/HSC-kowalczyk/BCA-HSC.csv")


sds_wour <- slingshot(Y_wour, all_clusts, start.clus = "LT-HSC")
ps_wour <- sds_wour@assays@data@listData[["pseudotime"]]


order_wour_structdr <- fread("data/real-data/HSC-kowalczyk/BCA-StructDR-order.csv")
taus <- simple_score(order_wour_structdr, int_clusts)
#taus <- all_branches_score(order_wour_structdr, int_clusts)


taus <- all_branches_score(ps_wour, int_clusts)

gene_names <- colnames(X)
marker_genes <- c("Itga6", "Itgb1bp1", "Itgb2")
g <- X[,colnames(X) == high_cor_genes[5]]
y <- cor(ps_wour, X)
high_cor_genes <- gene_names[abs(y) > 0.6]

# [1] "Plac8"     "Cd34"      "Mpo"       "Nkg7"      "Sell"     
# [6] "Cd48"      "Serpinb1a" "Mcm5"      "Gpr97"     "Lig1"     
# [11] "Mcm2"      "Tk1"       "BC035044"  "Pbk"       "Top2a"    
# [16] "Uhrf1" 


(p_bca <- plot_sling_curve(Y_wour, all_clusts, sds_wour, "BCC", 
                 "BCA with Slingshot on Mouse HSCs (n=873)",
                 clust_order, "Cell Type", 3))

  

c_mat_mbcve <- confusion_mat(ps_wour, clusts_to_ints(all_clusts, clust_order))
cnf_mat_mbcve <- c_mat_mbcve[[1]]$table
cat("Classification Accuracy", c_mat_mbcve[[1]]$overall['Accuracy'], "\n")
mbcve_title <- "Confusion Matrix for MBCVE and Slingshot on HSC Data"
plot_cnf_mat(cnf_mat_mbcve, clust_order, mbcve_title)

#########################################################################
############ UMAP ######################
Y_umap <- umap(X)
fwrite(Y_umap$layout, "data/real-data/HSC-kowalczyk/UMAP-HSC.csv")
sds_umap <- slingshot(Y_umap$layout, clusterLabels = clusts, start.clus = "MPP")
ps_umap <- sds_umap@assays@data@listData[["pseudotime"]]


order_umap_structdr <- fread("data/real-data/HSC-kowalczyk/UMAP-StructDR-order.csv")$V1
taus <- simple_score(order_umap_structdr, int_clusts)


taus <- all_branches_score(ps_umap, int_clusts)

plot_sling_curve(Y_umap$layout, clusts, sds_umap, "UMAP", 
                 "MBCVE with Slingshot on Human Primordial Germ Cell data (n=272)",
                 "Cell Type", clust_order, 3)

################### LDA ##################
Y_lda <- do.lda(Y_pca, as.factor(clusts), ndim=2)
fwrite(Y_lda$Y, "data/real-data/HSC-kowalczyk/LDA-HSC.csv")

sds_lda <- slingshot(Y_lda$Y, clusterLabels = all_clusts, start.clus = "LT-HSC")
ps_lda <- sds_lda@assays@data@listData[["pseudotime"]]

order_lda_structdr <- fread("data/real-data/HSC-kowalczyk/LDA-StructDR-order.csv")
taus <- simple_score(order_lda_structdr, int_clusts)


taus <- all_branches_score(ps_lda, int_clusts)

p_lda <- plot_sling_curve(Y_lda$Y, all_clusts, sds_lda, "LD", 
                 "mLDA with Slingshot on Mouse HSCs (n=873)",
                 "Cell Type", clust_order, 3)

c_mat_lda <- confusion_mat(ps_lda, clusts_to_ints(all_clusts, clust_order))
cnf_mat_lda <- c_mat_lda[[1]]$table
cat("Classification Accuracy", c_mat_lda[[1]]$overall['Accuracy'], "\n")
lda_title <- "Confusion Matrix for mLDA and Slingshot on HSC Data"
plot_cnf_mat(cnf_mat_lda, clust_order, lda_title)


p <- plot_grid(p_pca, p_bca, p_lda, ncol=3)
save_plot(paste0("pca-bca-lda-slingshot-Mouse-HSCs.pdf"), p, ncol=3)











