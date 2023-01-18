setwd("~/Desktop/Research/structured-priors/experiments/cell-cycle/")
source("utils.R")

data <- fread("data/mouse-gonad-hvg.csv")
clusters <- fread("data/mouse-gonad-clusters.csv")$V2

all_clusts <- clusters
X <- data
all_clusts[all_clusts == "early_supporting"] <- "ESGC"
all_clusts[all_clusts == "CoelEpi_Lhx9"] <- "CoelEpi"
unique(all_clusts)
clusts <- all_clusts

# Comment out the line below and rerun the script to obtain results with no occluded cell types.
# More elegant ways exist to handle this but this is fine in an R script. 
# clusts[all_clusts == "ESGC"] <- "Sertoli"

# clusts[all_clusts == "early_supporting"] <- "Gi"
# clusts[all_clusts == "early_supporting"] <- "sPAX8"


clust_order <- c("CoelEpi", "Gi", "sPAX8", "ESGC", "Sertoli") # remap these cluster names
int_clusts <- clusts_to_ints(all_clusts, clust_order)
clust_pairs <- compute_clust_pairs(clust_order)

# Color Palette for cell types
col_pal <- c("#F8766D", "#00BF7D", "#A3A500", "#00B0F6", "#E76BF3")

#######################################################

X_svd <- svd(scale(X, center=TRUE, scale=FALSE))
Y_pca <- X_svd$u %*% diag(X_svd$d)
sds_pca <- slingshot(Y_pca[,1:50], clusts, start.clus = clust_order[1])
ps_pca <- slingPseudotime(sds_pca)
pca_good_inds_1 <- !is.na(ps_pca[,1])
pca_good_inds_2 <- !is.na(ps_pca[,2])

# n = 1,840
(p_pca <- plot_sling_curve(Y_pca, all_clusts, sds_pca, "PC", 
                           "", clust_order,
                           "Cell Type", 3, FALSE))
ggsave(paste0("plots/gonad-plots/pca-2d-mouse-gonad.pdf"), p_pca, width=11, height=8.5)

(pca_ovl_box_1 <- plot_boxplot_subset(ps_pca[,1], all_clusts, pca_good_inds_1))
ggsave(paste0("plots/gonad-plots/pca-mouse-gonad-lin1-box.pdf"), pca_ovl_box_1, width=11, height=8.5)
(pca_ovl_box_2 <- plot_boxplot_subset(ps_pca[,2], all_clusts, pca_good_inds_2))
ggsave(paste0("plots/gonad-plots/pca-mouse-gonad-lin2-box.pdf"), pca_ovl_box_2, width=11, height=8.5)

pca_ovl_plots <- plot_grid(pca_ovl_plot1, pca_ovl_plot2, ncol=2)


(fig_pca <- plot_ly(x = Y_pca[,1], y = Y_pca[,2], z = Y_pca[,3], type="scatter3d", color = all_clusts))

#######################################################################################
Y_bca <- BCA(X, clusts)[,1:3]
sds_bca <- slingshot(Y_bca, clusts, start.clus = clust_order[1])
ps_bca <- slingPseudotime(sds_bca)
bca_good_inds_1 <- !is.na(ps_bca[,1])
bca_good_inds_2 <- !is.na(ps_bca[,2])


(p_bca <- plot_sling_curve(Y_bca, all_clusts, sds_bca, "BCC", 
                           "", clust_order,
                           "Cell Type", 3, FALSE))
ggsave(paste0("plots/gonad-plots/bca-2d-mouse-gonad.pdf"), p_bca, width=11, height=8.5)

(bca_ovl_box_1 <- plot_boxplot_subset(ps_bca[,1], all_clusts, bca_good_inds_1))
ggsave(paste0("plots/gonad-plots/bca-mouse-gonad-lin1-box-ES-to-Sert.pdf"), bca_ovl_box_1, width=11, height=8.5)
(bca_ovl_box_2 <- plot_boxplot_subset(ps_bca[,2], all_clusts, bca_good_inds_2))
ggsave(paste0("plots/gonad-plots/bca-mouse-gonad-lin2-box-ES-to-Sert.pdf"), bca_ovl_box_2, width=11, height=8.5)



########################################################################################
Y_lda <- do.lda(Y_pca[,1:50], as.factor(clusts), ndim=3)
sds_lda <- slingshot(Y_lda$Y, clusterLabels = clusts, start.clus = clust_order[1])
ps_lda <- slingPseudotime(sds_lda)

(p_lda <- plot_sling_curve(Y_lda$Y, all_clusts, sds_lda, "LD", 
                           "", clust_order,
                           "Cell Type", 3, FALSE))
ggsave(paste0("plots/gonad-plots/lda-2d-mouse-gonad.pdf"), p_lda, width=11, height=8.5)

lda_ovls <- consecutive_OVL(ps_lda, all_clusts, clust_pairs)
(lda_ovl_box <- plot_boxplot(ps_lda, all_clusts, clust_order, lda_ovls, ""))
ggsave(paste0("plots/gonad-plots/lda-mouse-gonad-box-ES-to-Sert.pdf"), lda_ovl_box, width=11, height=8.5)

# p2 <- plot_boxplot(ps_ours, all_clusts, clust_order, ours_ovls, "BCA mouse gonad; ES relabeled as Sertoli")+ coord_flip()
# p <- plot_grid(p2, p3, ncol=2)
# p
# save_plot(paste0("plots/pseudotime-boxplots-mouse-gonad-ES-to-sPAX8.pdf"), p, ncol=2)


X$clusters <- clusters
centroids <- X %>% group_by(clusters) %>% summarise_all("median") %>% remove_rownames %>% column_to_rownames(var="clusters")
dist(centroids)











library(plot3D)

plot3d(Y_wour, col=get_colors(clusts, unique(colors_to_use)$fill), size=8)
legend3d("topright", legend = c('Coel', 'ESGC', 'Gi', 'Sertoli', 'sPAX8'), pch = 16, col = unique(colors_to_use)$fill, cex=1, inset=c(0.02))
plot3d.SlingshotDataSet(SlingshotDataSet(sds_wour), lwd = 5, add = TRUE)


sling_3d_plot <- function(Y, clusts, cols, sds, axis_labels, clust_order, leg_pos="topright") {
  plot3d(Y, col=get_colors(clusts,cols, clust_order), size=8, xlab="", ylab="", zlab="", cex=3)
  # legend3d(leg_pos, legend = clust_order, pch = 16, col = cols, inset=c(0.002))
  plot3d.SlingshotDataSet(SlingshotDataSet(sds), lwd = 5, add = TRUE)
  title3d(xlab = axis_labels[1], ylab = axis_labels[2], zlab = axis_labels[3])
}

# 1 Coel; 2 ESGC; 3 Gi; 4 Sertoli; 5 sPAX8

par3d(cex=1.0)
pc_axes <- paste("PC", c("1", "2", "3"))
bca_axes <- paste("BCC", c("1", "2", "3"))
lda_axes <- paste("LD", c("1", "2", "3"))
sling_3d_plot(Y_pca[,1:3], all_clusts, col_pal, sds_pca, pc_axes, clust_order) 
# have to manually rotate image and then save; 
# there are better 3d packages none of them support adding slingshot curves to them
um <- par3d()$userMatrix
view3d(userMatrix = um)
rgl.postscript("plots/gonad-plots/pca-sling-mouse-gonad-3d.pdf", fmt = "pdf")

sling_3d_plot(Y_bca, all_clusts, col_pal, sds_bca, bca_axes, clust_order)
um <- par3d()$userMatrix
view3d(userMatrix = um)
rgl.postscript("plots/gonad-plots/bca-sling-mouse-gonad-3d-ES-to-Sert.pdf", fmt = "pdf")

sling_3d_plot(Y_lda$Y, all_clusts, col_pal, sds_lda, lda_axes, clust_order)
um <- par3d()$userMatrix
view3d(userMatrix = um)
rgl.postscript("plots/gonad-plots/lda-sling-mouse-gonad-3d-ES-to-Sert.pdf", fmt = "pdf")


