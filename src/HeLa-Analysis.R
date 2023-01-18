setwd("~/Desktop/Research/structured-priors/experiments/cell-cycle/")
source("utils.R")

library(Revelio)
myData <- createRevelioObject(rawData = revelioTestData_rawDataMatrix,
                              cyclicGenes = revelioTestData_cyclicGenes)
myData <- getCellCyclePhaseAssignInformation(dataList = myData) # need to run this to obtain their inferred cell phase assignments
myData <- getPCAData(dataList = myData)
myData <- getOptimalRotation(dataList = myData) # Compute Revelio Embedding aka DCs

pcData <- as.data.frame(t(myData@transformedData$pca$data)[,1:50]) # PCs
colnames(pcData) <- NULL # Necessary for plotting slingshot curve
revProj <- as.data.frame(t(myData@transformedData$dc$data))[,1:2]  # DCs
colnames(revProj) <- NULL # Necessary for plotting slingshot curve
cell_phase_assignments <- myData@cellInfo$ccPhase
clusters <- as.vector(cell_phase_assignments)
cell_phase_order <- c("G1.S", "S", "G2", "G2.M", "M.G1") # Arbitrary ordering for plotting confusion mats


# PCA + SLingshot
sds_pca <- slingshot(pcData, clusters)
ps_pca <- sds_pca@assays@data@listData[["pseudotime"]]
best_t_pca <- compute_score_cyclic(ps_pca, clusters) 
(p_pca <- plot_sling_curve(pcData, clusters, sds_pca, "PC", 
                 "", cell_phase_order, "Cell Phase", 3, FALSE) + guides(fill="none", color="none"))
ggsave(paste0("plots/hela-plots/pca-sling-hela.pdf"), p_pca, width=11, height=8.5)
best_order_pca <- c("G2", "G2.M", "M.G1", "G1.S", "S") # Found from compute_score_cyclic (printed out)

c_mat_pca <- confusion_mat(ps_pca, clusts_to_ints(clusters, best_order_pca))
cnf_mat_pca <- c_mat_pca[[1]]$table
cat("Classification Accuracy", c_mat_pca[[1]]$overall['Accuracy'], "\n")
pc_title <- "Confusion Matrix for PCA and Slingshot on HeLa Cell Data"
plot_cnf_mat(cnf_mat_pca, cell_phase_order, pc_title)

# Revelio + SLingshot

sds_rev <- slingshot(revProj, clusterLabels = clusters)
ps_rev <- sds_rev@assays@data@listData[["pseudotime"]]
best_t_rev <- compute_score_cyclic(ps_rev, clusters) 
p_rev <- plot_sling_curve(revProj, clusters, sds_rev, "DC", 
                 "",  cell_phase_order, "Cell Phase", 3, FALSE) + guides(fill="none", color="none")
ggsave(paste0("plots/hela-plots/revelio-sling-hela.pdf"), p_rev, width=11, height=8.5)

best_order_rev <- c("G1.S", "S", "G2", "G2.M", "M.G1")
c_mat_rev <- confusion_mat(ps_rev, clusts_to_ints(clusters, best_order_rev))
cnf_mat_rev <- c_mat_rev[[1]]$table
cat("Classification Accuracy", c_mat_rev[[1]]$overall['Accuracy'], "\n")
rev_title <- "Confusion Matrix for Revelio and Slingshot on HeLa Cell Data"
plot_cnf_mat(cnf_mat_rev, cell_phase_order, rev_title)

#  BCA + SLingshot
X <- myData@DGEs$scaledData[myData@geneInfo$variableGenes,]
Y_bca <- BCA(t(X), clusters)

sds_bca <- slingshot(Y_bca, clusterLabels = clusters)
ps_bca <- sds_bca@assays@data@listData[["pseudotime"]]
best_t <- compute_score_cyclic(ps_bca, clusters) 
p_bca <- plot_sling_curve(Y_bca, clusters, sds_bca, "BCC", 
                 "", cell_phase_order, "Cell Phase", 3, FALSE) + guides(fill="none", color="none")
ggsave(paste0("plots/hela-plots/bca-sling-hela.pdf"), p_bca, width=11, height=8.5)
best_order_mbcve <- c("G2", "G2.M", "M.G1", "G1.S", "S")

c_mat_mbcve <- confusion_mat(ps_our, clusts_to_ints(clusters, best_order_mbcve))
cnf_mat_mbcve <- c_mat_mbcve[[1]]$table
cat("Classification Accuracy", c_mat_mbcve[[1]]$overall['Accuracy'], "\n")
mbcve_title <- "Confusion Matrix for MBCVE and Slingshot on HeLa Cell Data"
plot_cnf_mat(cnf_mat_mbcve, cell_phase_order, mbcve_title)


###################### LDA #######################
Y_lda <- do.lda(as.matrix(pcData), as.factor(clusters))
sds_lda <- slingshot(Y_lda$Y, clusterLabels = clusters)
ps_lda <- sds_lda@assays@data@listData[["pseudotime"]]
best_t <- compute_score_cyclic(ps_lda, clusters) 
best_order_lda <- c("G2.M", "M.G1", "G1.S", "S", "G2")
int_clusts <- clusts_to_ints(clusters, best_order_lda)
taus <- all_branches_score(ps_lda, int_clusts)

(p_lda <- plot_sling_curve(Y_lda$Y, clusters, sds_lda, "LD", 
                 "", cell_phase_order,
                 "Cell Phase", 3))
ggsave(paste0("plots/hela-plots/lda-sling-hela.pdf"), p_lda, width=11, height=8.5)


####################### UMAP ########################
Y_umap <- umap(X)
sds_umap <- slingshot(Y_umap$layout, clusterLabels = clusts)
ps_umap <- sds_umap@assays@data@listData[["pseudotime"]]
best_t <- compute_score_cyclic(ps_umap, clusters) 
taus <- all_branches_score(ps_umap, int_clusts)
best_order_umap <- c("G2", "G2.M", "M.G1", "G1.S", "S")

plot_sling_curve(Y_umap$layout, clusters, sds_umap, "UMAP", 
                 "", cell_phase_order,
                 "Cell Phase", 3)




