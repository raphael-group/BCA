# if(!exists("foo", mode="function")) source("util.R")
setwd("~/Desktop/Research/structured-priors/experiments/cell-cycle/")
source("utils.R")
library(dynwrap)
library(dyneval)
library(dynmethods)

lda_eval <- function(data) {
  X <- data$expression
  clusts <- data$prior_information$groups_id$group_id
  X_svd <- svd(scale(X, center=TRUE, scale=FALSE))
  Y_pca <- X_svd$u %*% diag(X_svd$d)
  Y_lda <- do.lda(Y_pca[,1:50], as.factor(clusts))$Y
  Y_ours <- wour_method(X, clusts)[,1:2]
  rownames(Y_lda) <- rownames(X)
  rownames(Y_pca) <- rownames(X)
  
  Y_lda_clusts <- as.data.frame(Y_lda)
  Y_lda_clusts$clusts <- clusts
  gdf <- Y_lda_clusts %>%
    group_by(clusts) %>%
    summarise_at(vars(V1, V2), funs(mean(., na.rm=TRUE))) %>% data.frame
  
  Y_lda_clusts2 <- as.data.frame(Y_ours)
  Y_lda_clusts2$clusts <- clusts
  gdf2 <- Y_lda_clusts2 %>%
    group_by(clusts) %>%
    summarise_at(vars(V1, V2), funs(mean(., na.rm=TRUE))) %>% data.frame
  
  rownames(gdf) <- gdf$clusts
  gdf <- gdf %>% select(-clusts) %>% as.matrix
  
  rownames(gdf2) <- gdf2$clusts
  gdf2 <- gdf2 %>% select(-clusts) %>% as.matrix
  
  trajectory <- add_dimred_projection(
    data,
    milestone_network = data$milestone_network,
    dimred = Y_ours[,1:2], 
    dimred_milestones = gdf2
  )
  
  cat("correlation is ", calculate_metrics(data, trajectory, metrics = "correlation",expression_source = data$expression)$correlation)
  
  
  
  if ("dynplot" %in% rownames(installed.packages())) {
    dynplot::plot_dimred(trajectory, color_cells = "grouping", expression_source = as.matrix(dataset$expression))
  } 
}

data <- readRDS("data/real-data/aging-hsc-old_kowalczyk.rds")
data2 <- readRDS("data/real-data/germline-human-both_guo.rds")
data3 <- readRDS("data/real-data/human-embryos_petropoulos.rds")

lda_eval(data)
lda_eval(data2)
lda_eval(data3)



X <- data$expression
clusts <- data$prior_information$groups_id$group_id

def <- definition(
  method = def_method(
    id = "naive"
  ),
  wrapper = def_wrapper(
    input_required = "expression",
    input_required = "groups_id",
    input_optional = "start_id"
  )
)

run_fun <- function(expression, priors, parameters, seed, verbose) {
  pseudotime <- as.vector(priors$groups_id[,2])
  
  
  
  dynwrap::wrap_data(cell_ids = rownames(expression)) %>%
    dynwrap::add_linear_trajectory(pseudotime = pseudotime)
}
ti_naive <- create_ti_method_r(def, run_fun, package_loaded = "dplyr")
trajectory <- infer_trajectory(data, ti_naive(), give_priors="groups_id")



trajectory <- add_cluster_graph(data2, data2$milestone_network)
calculate_metrics(data2, trajectory, metrics = "correlation",expression_source = data2$expression)$correlation


X_svd <- svd(scale(X, center=TRUE, scale=FALSE))
Y_pca <- X_svd$u %*% diag(X_svd$d)
Y_lda <- do.lda(Y_pca[,1:50], as.factor(clusts))$Y
rownames(Y_lda) <- rownames(data$expression)

Y_lda_clusts <- as.data.frame(Y_lda)
Y_lda_clusts$clusts <- clusts
gdf <- Y_lda_clusts %>%
  group_by(clusts) %>%
  summarise_at(vars(V1, V2), funs(mean(., na.rm=TRUE))) %>% data.frame

rownames(gdf) <- gdf$clusts
gdf <- gdf %>% select(-clusts) %>% as.matrix

trajectory <- add_dimred_projection(
  data,
  milestone_network = data$milestone_network,
  dimred = Y_lda, 
  dimred_milestones = gdf
)

calculate_metrics(data, trajectory, metrics = "correlation",expression_source = data$expression)$correlation

traj <- dynwrap::infer_trajectories(data, ti_slingshot(), give_priors=c("groups_id"))

dataset <- data <- dyntoy::generate_dataset(
   num_cells = 99,
   num_features = 101,
   model = "tree",
   normalise = FALSE
)
 model <- dynwrap::infer_trajectory(data, ti_slingshot())
 
 pca <- irlba::prcomp_irlba(dataset$expression, n = 5)
 
library(slingshot)
sds <- slingshot(Y_pca[,1:50], clusts)

if ("dynplot" %in% rownames(installed.packages())) {
  dynplot::plot_dimred(trajectory,  color_cells = "grouping", expression_source = as.matrix(dataset$expression))
}
 
 install.packages("/Users/alexanderstrzalkowski/Desktop/Research/structured-priors/experiments/cell-cycle/understand-dynverse/ti_slingshot/package", 
                  repos = NULL, 
                  type = "source")
 library(tislingshot)
 library(dyntoy)
