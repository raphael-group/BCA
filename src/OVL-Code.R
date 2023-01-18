setwd("/Users/alexanderstrzalkowski/Desktop/Research/structured-priors/experiments/cell-cycle/gap-stats/")
source("utils.r")
library(dynwrap)
library(dyneval)
library(slingshot)
library(Rdimtools)
library(data.table)
library(cowplot) # for putting multiple plots on grid
library(stringr)
library(dplyr)
library(ggpubr) # for stat_compare_means + ggboxplot 
library("bayestestR") # for overlap (implementation of OVL using KDE)
library(stringr)

oldw <- getOption("warn")
options(warn = -1) # to silence warnings from lm

# "real/gold/human-embryos_petropoulos"

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

plot_boxplot <- function(p, clusts, clust_order, ovls, title, num_digits=3) {
  df <- data.frame(Pseudotime = as.vector(p), clusters = clusts)
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
    geom_bracket(label=ovls_rounded, y.position = 110, step.increase = 0.1, 
                 xmin=min_seq, xmax=max_seq, coord.flip = TRUE, show.legend=FALSE, label.size= 6) + 
    rotate_x_text(angle=15) + 
    ggtitle(title) + 
    theme_classic(base_size=24)
  #geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  #geom_jitter(position=position_jitter(width=.1, height=0)) + 
  #ggtitle(title)
}



real_gold_data_path <- "../reproduce-benchmark/data/dynverse_data_v2/real/gold/"
files <- list.files(path= real_gold_data_path, pattern = "\\.rds$")

n <- length(files) # 
f1s_pca <- rep(0, n)
f1s_ours <- rep(0, n)
f1s_lda <- rep(0, n)
 # seq_along(files) # c(1, 2, 8, 14, 15, 17, 19, 25)
for (i in c(8)) {
  #cat(paste0(i, "\n"))
  data <- readRDS(paste0(real_gold_data_path, files[8]))
  if (data$trajectory_type != "linear") next 
  clusts <- data$prior_information$groups_id$group_id
  # The ordering is embedded into the Dynverse datastructure in the milestone_network slot
  # consisting of a two column df: `from` and `to` indicating the of which clusters lead to others.
  # For linear trajectories the ordering is simply given by the whole `from` column and appending 
  # the last entry of `to`. (I'm prety sure this makes sense, and qualitatively at least holds).
  clust_order <- c(as.vector(data$milestone_network$from), data$milestone_network$to[length(data$milestone_network$to)])
  clust_pairs <- Map(c, clust_order[-length(clust_order)], clust_order[-1]) # https://stackoverflow.com/questions/65979751/r-iterate-over-consecutive-pairs-in-a-list
  
  
  norm_scaling <- 100
  
  use_start_clust <- FALSE
  Y_pca <- run_pca50(data)[,1:2]
  sds_pca <- run_slingshot(data, Y_pca, use_start_clust)
  p_pca <- slingPseudotime(sds_pca)
  p_pca <- norm_scaling * (p_pca - min(p_pca)) / max(p_pca)
  Y_lda <- run_lda(data)[,1:2]
  sds_lda <- run_slingshot(data, Y_lda, use_start_clust)
  p_lda <- slingPseudotime(sds_lda)
  p_lda <- norm_scaling * (p_lda- min(p_lda)) / max(p_lda)
  Y_ours <- run_ours(data)[,1:2]
  sds_ours <- run_slingshot(data, Y_ours, use_start_clust)
  p_ours <- slingPseudotime(sds_ours)
  p_ours <- norm_scaling * (p_ours - min(p_ours)) / max(p_ours)
  
  num_lineages <- c(dim(p_pca)[2], dim(p_lda)[2], dim(p_ours)[2])
  
  if (any(num_lineages > 1)){
    cat(paste0("Warning: PCA, LDA, or BCA caused slingshot to find more than one lineage in a linear Trajectory. Will not make plots."))
  }
  
  
  if(!any(num_lineages > 1)) {
    data_name <- str_split(data$id, '/')[[1]][3]
    cat(paste0(i, "\n"))
    
    pca_ovls <- consecutive_OVL(p_pca, clusts, clust_pairs)
    lda_ovls <- consecutive_OVL(p_lda, clusts, clust_pairs)
    ours_ovls <- consecutive_OVL(p_ours, clusts, clust_pairs)
    
    ovl_df <- data.frame(PCA_ovls=pca_ovls, BCA_ovls=ours_ovls, LDA_ovls=lda_ovls)
    row.names(ovl_df) <- unlist(lapply(clust_pairs, function(s) paste(s, collapse="-")))
    #fwrite(ovl_df, paste0("ovl-results/", "slinghot-ovl-", data_name, ".csv"))
    
    # tit_pca <- "Slingshot Pseudotimes by Cluster using PCA"
    # tit_bca <- "Slingshot Pseudotimes by Cluster using BCA"
    # tit_lda <- "Slingshot Pseudotimes by Cluster using mLDA"
    
    tit_pca <- ""
    tit_bca <- ""
    tit_lda <- ""
    
    (p1 <- plot_boxplot(p_pca, clusts, clust_order, pca_ovls, tit_pca) + rremove("legend"))
    (p2 <- plot_boxplot(p_ours, clusts, clust_order, ours_ovls, tit_bca) + rremove("legend"))
    (p3 <- plot_boxplot(p_lda, clusts, clust_order, lda_ovls, tit_lda) + rremove("legend"))
    p <- plot_grid(p1, p2, p3, ncol=3)
    p
    save_plot(paste0("plots/pseudotime-boxplots-", data_name, "-norm-100-coord-flip.pdf"), p, ncol=3)
    
    #f1s_pca[i] <- F1_macro(p_pca, clusts, clust_order)
    #f1s_ours[i] <- F1_macro(p_ours, clusts, clust_order)
    #f1s_lda[i] <- F1_macro(p_lda, clusts, clust_order)
  }

}
f1s_pca <- f1s_pca[f1s_pca != 0]
f1s_ours <- f1s_ours[f1s_ours != 0]
f1s_lda <- f1s_lda[f1s_lda != 0]
df <- data.frame(f1s_pca, f1s_ours, f1s_lda)

# data <- readRDS(paste0(real_gold_data_path, files[1]))
# clusts <- data$prior_information$groups_id$group_id




#data <- readRDS("human-embryos_petropoulos.rds")
clusts <- data$prior_information$groups_id$group_id
#clust_order <- unique(clusts) # this works for dataset 15 just because it is like day 3, 4, 5, ...


library(dplyr)
library(readr)
# combine all separate csvs of ovl results into df and melt to make easier for plotting boxplots
df <- list.files(path="ovl-results/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows  %>%
  melt(variable.name="Method", value.name="OVL")

my_comparisons <- list(c("PCA_ovls", "BCA_ovls"), c("PCA_ovls", "LDA_ovls"), c("BCA_ovls", "LDA_ovls"))
ggboxplot(df, x="Method", y="OVL", add="jitter") + 
  stat_compare_means(comparisons=my_comparisons, method="t.test", paired=TRUE, size=5) + 
  ggtitle("") + 
  scale_x_discrete(labels = c('PCA','BCA','mLDA')) + 
  theme_classic(base_size=24)


df %>% group_by(Method) %>% summarize(Mean = mean(OVL), Variance = var(OVL))




# ggboxplot(df, x = "dose", y = "len") +
#   geom_bracket(
#     xmin = c("0.5", "1"), xmax = c("1", "2"),
#     y.position = c(30, 35), label = c("***", "**"),
#     tip.length = 0.01
#   )



var_df <- data.frame(PCA=pca_var_sums, LDA=lda_var_sums, BCA=ours_var_sums)
plot_var_df <- melt(var_df, value.name="Total_Variance")
colnames(plot_var_df)[1] <- "Method"

plot_var_df %>%
  ggplot( aes(x=Method, y=Total_Variance)) +
  geom_boxplot() +scale_y_log10()+ 
  ylab("log(Total Variance)") + ggtitle("Total variance computed for each cluster for each real dataset in PCA/LDA/BCA space (no ptime)")


compute_clust_stats <- function(ps, clusts, just_var=TRUE) {
  df <- data.frame(pseudotime = sort(ps), clusters = clusts, cell_number=seq(1:length(ps)))
  names(df)[names(df) == "Lineage1"] <- "ptime"
  if (just_var) {
    out_df <- df %>% group_by(clusters) %>% 
      summarise(var_pseudotime = var(pseudotime))
  } else {
    out_df <- df %>% group_by(clusters) %>% 
      summarise(mean_pseudotime = mean(pseudotime), var_pseudotime = var(pseudotime))
  }
  out_df
}
# mean_pseudotime = mean(pseudotime), std_pseudotime = sd(pseudotime),


plot_pseudotimes <- function(ps, clusts, title){
  df <- data.frame(pseudotime = sort(ps), clusters = clusts, cell_number=seq(1:length(ps)))
  names(df)[names(df) == "Lineage1"] <- "ptime"
  
  ggplot(df, aes(x=cell_number, y=ps))+ 
    geom_point(aes(fill = clusters), shape=21, size=3) + ggtitle(title) + theme_classic() + 
    ylab("Pseudotime")
}

plot_pseudotime_1d <- function(ps, clusts, title) {
  df <- data.frame(pseudotime = sort(ps), clusters = clusts, cell_number=seq(1:length(ps)))
  names(df)[names(df) == "Lineage1"] <- "ptime"
  
  ggplot(df, aes(x=ps))+ 
    geom_dotplot(aes(fill = clusters)) + ggtitle(title) + theme_classic() + 
    xlab("Pseudotime")
}


plot_pseudotimes(p_pca, clusts, "PCA Pseudotimes")
plot_pseudotimes(p_lda, clusts, "LDA Pseudotimes")
plot_pseudotimes(p_ours, clusts, "BCA Pseudotimes")

plot_pseudotime_1d(p_pca, clusts, "PCA Pseudotimes")
plot_pseudotime_1d(p_lda, clusts, "LDA Pseudotimes")
plot_pseudotime_1d(p_ours, clusts, "BCA Pseudotimes")





p1_upper_quantile <- df %>% 
  filter(clusters == clust_order[1]) %>%
  summarise(third_quant = unname(quantile(ptime)[4])) # 75th quantile 


get_std_sum <- function(p, clusts) {
  df <- data.frame(pseudotime = p, clusters = clusts)
  sds <- df %>% group_by(clusters) %>% summarise(sd(Lineage1)) 
  sum(sds[,2])
}

get_var_sum <- function(p, clusts) {
  df <- data.frame(pseudotime = p, clusters = clusts)
  vars <- df %>% group_by(clusters) %>% summarise(var(Lineage1))
  sum(vars[,2])
}

get_std_sum_embed <- function(Y, clusts) {
  df <- data.frame(Y)
  df$clusters <- clusts
  sds <- df %>% group_by(clusters) %>% summarise(across(everything(), list(sd))) 
  sum(sds[,2:3]) # first column is clusts the next two are just the first two dimensions
}

get_var_sum_embed <- function(Y, clusts) {
  df <- data.frame(Y)
  df$clusters <- clusts
  vars <- df %>% group_by(clusters) %>% summarise(across(everything(), list(var))) 
  sum(vars[,2:3]) # first column is clusts the next two are just the first two dimensions
}

get_gap <- function(Y, p, clusts, clust_pairs) {
  df <- data.frame(pseudotime = p, clusters = clusts)
  names(df)[names(df) == "Lineage1"] <- "ptime"
  for(pair in clust_pairs) {
    c1 <- pair[1]
    p1_upper_quantile <- df %>% 
      filter(clusters == c1) %>%
      summarise(third_quant = unname(quantile(ptime, c(0.75)))) # 75th quantile 
    
    index1 <- which.min(abs(df$ptime-p1_upper_quantile$third_quant))
    point_1 <- Y[index1, ]
    
    c2 <- pair[2]
    p2_lower_quantile <- df %>% 
      filter(clusters == c2) %>%
      summarise(first_quant = unname(quantile(ptime, c(0.25)))) # 25th quantile 
    
    index2 <- which.min(abs(df$ptime-p2_lower_quantile$first_quant))
    point_2 <- Y[index2, ]
    
    cat(paste0("Distance between ", c1, " and ", c2, "\n"))
    cat(paste0(sqrt(sum((point_1 - point_2)^2)), "\n"))
  }
}




pc_title <- "PCA on Human DC cells (n=238)"
bc_title <- "BCA on Human DC cells (n=238)"
lda_title <- "mLDA on Human DC cells (n=238)"

cs <- c("blue", "green", "red")
cs <- c("#619CFF", "#00BA38", "#F8766D")
ls <- c("PreDC", "CDP", "MDP")

source("utils.r")

add_sling_curve <- function(p, sds, my_title) {
  curves <- slingCurves(sds, as.df = TRUE)
  p <- p + geom_path(data = curves %>% arrange(Order),aes(group = Lineage)) 
  p <- p + ggtitle(my_title) 
  p
}

(p_pca1 <- dimred_plot(Y_pca, clusts, "PC ", "", ls, cs, 5))
# p_pca1 <- add_sling_curve(p_pca1, sds_pca, pc_title)
(p_lda1 <- dimred_plot(Y_lda, clusts, "LD ", "", ls, cs, 5))
(p_ours1 <- dimred_plot(Y_ours, clusts, "BCC ", "", ls, cs, 5))

p2 <- plot_grid(p_pca1, p_ours1, p_lda1, ncol=3)
p2
save_plot(paste0("plots/dimplots-", data_name, "-.pdf"), p2, ncol=3)



overlap()






normalize <- function(xs, scaling=2) {
  scaling * ((xs - min(xs)) / (max(xs))) + 2
}

n <- 10000
x1 <- normalize(rnorm(n, 4)) 
x2 <- normalize(rnorm(n, 8))


overlapPlot(x1, x2, linetype = c(1,1), linecol = c("black", "black"), main="Zone 2")
Axis(side=1, labels=FALSE)
Axis(side=2, labels=FALSE)





plot_overlap <- function(p, c, which_pair = c("PreDC", "CDP")) {
  df <- data.frame(pseudotime = as.vector(p), clusters = c)
  c1 <- df %>% 
    filter(clusters == which_pair[1]) %>% select(pseudotime)
  
  c2 <- df %>% 
    filter(clusters == which_pair[2]) %>% select(pseudotime)
  
  plot(overlap(c1$pseudotime, c2$pseudotime))
}
   
plot_overlap(p_lda, clusts, c("ST-HSC", "LT-HSC"))
axis(side=1, labels = FALSE)
title(xlab="test", ylab="test")

use_start_clust <- FALSE
Y_pca <- run_pca50(data)[,1:2]
sds_pca <- run_slingshot(data, Y_pca, use_start_clust)
p_pca <- slingPseudotime(sds_pca)
p_pca <- norm_scaling * (p_pca - min(p_pca)) / max(p_pca)
Y_lda <- run_lda(data)[,1:2]
sds_lda <- run_slingshot(data, Y_lda, use_start_clust)
p_lda <- slingPseudotime(sds_lda)
p_lda <- norm_scaling * (p_lda- min(p_lda)) / max(p_lda)
Y_ours <- run_ours(data)[,1:2]
sds_ours <- run_slingshot(data, Y_ours, use_start_clust)
p_ours <- slingPseudotime(sds_ours)
p_ours <- norm_scaling * (p_ours - min(p_ours)) / max(p_ours)


data <- readRDS(paste0(real_gold_data_path, files[1]))

clusts <- data$prior_information$groups_id$group_id
