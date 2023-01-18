library(tidyverse)
library(dplyr)
library(ggplot2)
library(slingshot)
library(data.table)
library(scales)
library(viridis)
library(glue)
library(caret) # for confusion mat
library(umap)
library(Rdimtools)
library(plyr) # for mapvalues function (esentials functions as a dictionary mapping; no idea how efficient)
library(cowplot) # for putting multiple plots on grid
library(bayestestR) # for overlap (implementation of OVL using KDE)
library(ggpubr) # for stat_compare_means + ggboxplot 
library(stringr)
library(rgl)

set.seed(0)



BCA <- function(X, clusts) {
  X <- scale(X, center=TRUE, scale=FALSE)
  M <- as_tibble(as.matrix(X)) %>% add_column(clusters = clusts) %>% group_by(clusters) %>% summarise_all("mean") %>% select(-clusters)
  
  clust_sizes <- as_tibble(as.matrix(X)) %>% add_column(clusters = clusts) %>% group_by(clusters) %>% dplyr::summarize(n())
  size_Cks <- diag(clust_sizes$`n()`)
  M_tilde <- size_Cks %*% data.matrix(M)
  Y <- X %*% svd(M_tilde)$v
}

BCA_mat<- function(X, clusts) {
  X <- scale(X, center=TRUE, scale=FALSE)
  M <- as_tibble(as.matrix(X)) %>% add_column(clusters = clusts) %>% group_by(clusters) %>% summarise_all("mean") %>% select(-clusters)
  
  clust_sizes <- as_tibble(as.matrix(X)) %>% add_column(clusters = clusts) %>% group_by(clusters) %>% summarize(n())
  size_Cks <- diag(clust_sizes$`n()`)
  M_tilde <- size_Cks %*% data.matrix(M)
  svd(M_tilde)$v
}

#' Assign a color to each cell based on some value
#' 
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

add_ids <- function(df) {
 df <- left_join(
    df,
    data$milestone_percentages %>%
      group_by(.data$cell_id) %>%
      arrange(desc(.data$percentage)) %>%
      filter(dplyr::row_number() == 1) %>%
      select(.data$cell_id, .data$milestone_id, .data$percentage),
    "cell_id"
  )
  df
}

clusts_to_ints <- function(clusts, order_of_clusts=c("MPP", "ST-HSC", "LT-HSC")) {
  ints=1:length(order_of_clusts)
  names(ints)=order_of_clusts
  int_clusts <- ints[clusts]
  int_clusts
}

plot_traj <- function(df, c, title) {
  c <- sym(c)
  ggplot(df, aes(x=V1, y=V2)) + geom_point(aes(colour=!!c)) + ggtitle(title)
}

get_score <- function(labels, guess_order) {
  guess_cluster_order <- labels[order(guess_order)]
  cor(guess_cluster_order, sort(guess_cluster_order), method="kendall")
}

get_lin_inds <- function(lins, which_lin) {
  lin_i <- lins[,which_lin]
  !is.na(lin_i)
}

all_branches_score <- function(ps, int_clusters) {
  num_lineages <- dim(ps)[2]
  taus <- c()
  for (l in 1:num_lineages) {
    curr_lineage_inds <- get_lin_inds(ps, l)
    labels <- unname(int_clusters[curr_lineage_inds])
    inferred_pseudotime <- unname(ps[curr_lineage_inds,l])
    tau <- get_score(labels, inferred_pseudotime)
    taus <- c(taus, tau)
    cat(glue("Kendall's Tau for Lineage {l} is {tau}."), "\n")
  }
  taus
}

simple_score <- function(ps, int_clusters) {
    labels <- unname(int_clusters)
    inferred_pseudotime <- ps
    tau <- get_score(labels, inferred_pseudotime)
    cat(glue("Kendall's Tau for Lineage 1 is {tau}."), "\n")
  tau
}

compute_score_cyclic <- function(ps, clusters) {
  cc_phases <- c("G1.S", "S", "G2", "G2.M", "M.G1")
  num_phases <- length(cc_phases)
  cyc_perm_phases <- list()
  cyc_perm_phases[[1]] <- cc_phases
  for(i in 2:num_phases) {
    tmp_phases <- c(cc_phases[i:num_phases], cc_phases[1:i-1])
    cyc_perm_phases[[i]] <- tmp_phases
  }
  cyc_taus <- list()
  for (i in 1:num_phases){
    curr_int_clusts <- clusts_to_ints(clusters, cyc_perm_phases[[i]])
    t_i <- all_branches_score(ps, curr_int_clusts)
    cyc_taus[[i]] <- t_i
  }
  cat("Best Kendall's Tau over all Cyclic Permuations is:", max(unlist(cyc_taus)), "\n")
  cat("Best Ordering is", cyc_perm_phases[[order(unlist(cyc_taus))[num_phases]]], "\n")
  max(unlist(cyc_taus))
}

confusion_mat <- function(ps, int_clusters) {
  num_lineages <- dim(ps)[2]
  c_mats <- list()
  for (l in 1:num_lineages) {
    curr_lineage_inds <- get_lin_inds(ps, l)
    # labels <- unname(int_clusters[curr_lineage_inds])
    # inferred_pseudotime <- unname(ps[curr_lineage_inds,l])
    
    labels <- int_clusters[curr_lineage_inds]
    inferred_pseudotime <- unname(ps[curr_lineage_inds,l])
    
    guess_cluster_order <- names(labels[order(inferred_pseudotime)])
    known_cluster_order <- names(sort(labels)) # sort by integer mapping than map back to str name
    c_mats[[l]] <- confusionMatrix(as.factor(guess_cluster_order), as.factor(known_cluster_order))
  }
  c_mats
}


dimred_plot <- function(Y, clusts, axis_title, clust_order, legend_title="Cell Type", point_size=3){
  plt_df <- data.frame(Y, "cl" = clusts)
  # Note that the columns must be Dim.insert_col_num to work especially in the next two
  # lines otherwise it returns a bizzare error 
  colnames(plt_df) <- c(paste0("Dim.", 1:(dim(plt_df)[2]-1)), "Cell_Type")
  p <- ggplot(plt_df, aes(x = Dim.1, y = Dim.2)) + 
    geom_point(aes(fill = Cell_Type), shape=21, size=point_size)+
    labs(x= paste0(axis_title, " 1"), y = paste0(axis_title, " 2"), fill=legend_title) + 
    scale_fill_discrete(limits = clust_order) + 
    theme_classic(base_size=24)
  p
}

plot_sling_curve <- function(Y, clusters, sds, axis_title, my_title, clust_order, legend_title="Cell Type", point_size=3,  do_save=FALSE) {
  p <- dimred_plot(Y, clusters, axis_title, clust_order, legend_title, point_size)
  curves <- slingCurves(sds, as.df = TRUE)
  p <- p + geom_path(data = curves %>% arrange(Order),aes(group = Lineage), inherit.aes = TRUE) 
  p <- p + ggtitle(my_title) 
  if (do_save) {
    ggsave(file = paste0(my_title, ".png"))
  }
  else{
    p
  }
}

########## Utilities for OVL experiments #####################

# Motivation for OVL functions that end in *_subset. 
# There was a need to plot OVL values between adjacent pairs of cell type pseudotimes
# for the case of multi-lineage (>1) trajectories from slingshot. 
# Since subsets of cells can exclude cell types this screws with the ordering 
# thus some special precautions had to be taken to make everything work. 


# Trick taken from https://stackoverflow.com/questions/65979751/r-iterate-over-consecutive-pairs-in-a-list
compute_clust_pairs <- function(clust_order) {
  Map(c, clust_order[-length(clust_order)], clust_order[-1])
}

# Given two vectors of strings: 
# * `f` i.e. from values and 
# * `to_dict` i.e. to values.
# Map the values in a vector `xs` that contains values in `f`, to values in `to_dict`,
# sort and map back. 
# Example: You want to order the following strings: "banna", "alligator", "kangaroo". 
# Input: f = c("banna", "alligator", "kangaroo"), to_dict = 1:3, xs = c("kangaroo", "banna", "alligator");
# Output: c("banna", "alligator", "kangaroo")
sort_str_by_int_map <- function(xs, f, to_dict=1:5) {
  mapvalues(sort(mapvalues(xs, from = f, to = to_dict)), from = to_dict, to = f)
}

# Function copied from (for 3d ploting): 
# http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization
get_colors <- function(groups, group.col, clust_order){
  groups <- factor(groups, levels=clust_order)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

consecutive_OVL <- function(p, clusts, clust_pairs) {
  df <- data.frame(pseudotime = as.vector(p), clusters = clusts)
  consec_ovls <- rep(0, length(clust_pairs))
  i <- 1
  for(pair in clust_pairs) {
    c1 <- df %>% 
      filter(clusters == pair[1]) %>% select(pseudotime)
    
    c2 <- df %>% 
      filter(clusters == pair[2]) %>% select(pseudotime)
    
    consec_ovls[i] <- overlap(c1, c2)
    i <- i + 1
  }
  consec_ovls
}

consecutive_OVL_subset <- function(ps, all_clusts, inds, clust_order) {
  consecutive_OVL(ps[inds], all_clusts[inds], compute_clust_pairs(clust_order))
}

plot_boxplot <- function(p, clusts, clust_order, ovls, title, num_digits=2) {
  # scale pseudotimes from 0 to 100 to make easierr to visualize
  norm_scaling <- 100
  p_scaled <- norm_scaling * (p - min(p)) / max(p)
  df <- data.frame(Pseudotime = as.vector(p_scaled), clusters = clusts)
  df$clusters <- factor(df$clusters, levels = clust_order)
  # geom_bracket requires specification of where to place the left side (min) of all brackets as well as
  # the right side (max) this is specified using min_seq and max_seq. 
  min_seq <- clust_order[1:length(clust_order)-1] 
  max_seq <- clust_order[2:length(clust_order)]
  # consec_clusts must be in the form: list(c("1", "2"), c("2", "3"),...), 
  # where each c(.) is a consecutive pair
  # consec_clusts <- list(c("embryonic day 3", "embryonic day 4"), c("embryonic day 4", "embryonic day 5"), c("embryonic day 5", "embryonic day 6"), c("embryonic day 6", "embryonic day 7"))# apply(cbind(unq_clusts[-length(unq_clusts)], unq_clusts[-1]), 1, list)
  # 
  # note that the y.position in geom_bracket is hard coded, will try to make this smarter later. 
  ovls_rounded <- sapply(ovls, function(x) round(x, digits=num_digits))
  p <- ggboxplot(df, x="clusters", y="Pseudotime", add="jitter", color="clusters", show.legend=FALSE, orientation="horizontal") + 
    geom_bracket(label=ovls_rounded, y.position = 105, step.increase = 0.1, 
                 xmin=min_seq, xmax=max_seq, coord.flip = TRUE, show.legend=FALSE, label.size= 8) + 
    rotate_x_text(angle=15) + 
    ggtitle(title) + 
    theme_classic(base_size=32) + 
    guides(fill="none", color="none") + # last option removes legend 
    xlab("") +
    theme(axis.text.y=element_text(size=45))
}

plot_boxplot_subset <- function(ps, clusts, inds) {
  f <- c("CoelEpi", "Gi", "sPAX8", "ESGC", "Sertoli")
  cell_type_order <- sort_str_by_int_map(unique(clusts[inds]), f)
  ovls <- consecutive_OVL_subset(ps, clusts, inds, cell_type_order)
  p <- plot_boxplot(ps[inds], clusts[inds], cell_type_order, ovls, "")
  p
}





















# Adapted from: https://stackoverflow.com/questions/7421503/how-to-plot-a-confusion-matrix-using-heatmaps-in-r
plot_cnf_mat <- function(cnf_mat, known_order, title) {
  p_df <- as_tibble(cnf_mat) # Automatically turns it into right format for plotting
  p_df <- p_df %>%
    mutate(Reference = factor(Reference, levels = known_order), # alphabetical order by default
           Prediction = factor(Prediction, levels = rev(known_order)))
  
  p <- ggplot(p_df, aes(x=Reference, y=Prediction, fill=n)) +
    geom_tile(color="black") + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="Greens", direction=1) +
    guides(fill="none") + # removing legend for `fill`
    labs(title = title) + # using a title instead
    geom_text(aes(label=n), color="black") # printing values
  p
}
