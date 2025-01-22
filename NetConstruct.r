# %%
library(phyloseq)
library(SPRING)
library(SpiecEasi)
library(caret)
library(pROC)
library(MLmetrics)
library(ggplot2)
library(tidyr)
library(ggraph)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
source("Packages\\stabENG.r")
source("Packages\\MyENG.r")
source("Packages\\stabENG.r")
rawdata <- readRDS("data\\DavarData1_substrate_phyloseq_1226_final_filtered.RDS")
otu_raw <- otu_table(rawdata)
otu_RA <- transform_sample_counts(otu_raw, function(x) x / sum(x) )
otu_RA <- filter_taxa(otu_RA, function(x) mean(x) > 1e-3, TRUE)
shared_otu <- rownames(otu_RA)
otu_Ab <- as.data.frame(t(otu_raw[shared_otu,]))
sam_info <- as.data.frame(sample_data(rawdata))
otu_tax <- as.data.frame(tax_table(rawdata))

# %%split otu_Ab by condition group
otu_Ab_Nplus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="plusN",]),]
otu_Ab_Nminus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="minusN",]),]
data_list <- list(Nplus = otu_Ab_Nplus, Nminus = otu_Ab_Nminus)
# by Group and Days
timestamps <- unique(sam_info$Days)
otu_Ab_Nplus_times <- list()
otu_Ab_Nminus_times <- list()
data_list_times <- list()
for (i in timestamps)
  {
    otu_Ab_Nplus_times[[i]] <- otu_Ab[rownames(otu_Ab) %in% 
      rownames(sam_info[sam_info$growthCondition=="plusN" & sam_info$Days == i,]),]
    otu_Ab_Nminus_times[[i]] <- otu_Ab[rownames(otu_Ab) %in%
      rownames(sam_info[sam_info$growthCondition=="minusN" & sam_info$Days == i,]),]
    data_list_times[[i]] <- list(Nplus = otu_Ab_Nplus_times[[i]], Nminus = otu_Ab_Nminus_times[[i]])
  }
# %%
i = "12"
data_list <- data_list_times[[i]]
cat('Calculating network on day ',i,'\n')
network_results <- stabENG(data_list, labels = shared_otu, var.thresh = 0.1, rep.num = 25,
  nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,
  lambda2.init=0.01,ebic.gamma=0.6)
network_Nplus <- network_results$opt.fit$Nplus # precision matrix estimates
network_Nminus <- network_results$opt.fit$Nminus # precision matrix estimates
# filter edge sparsity
network_Nplus[abs(network_Nplus) < 0.1] <- 0
network_Nminus[abs(network_Nminus) < 0.1] <- 0
diag(network_Nplus) = diag(network_Nminus) <- 0
# %% Plot network on Phylum level
Phylum_groups <- as.factor(otu_tax[rownames(network_Nplus),"Phylum"])
png(filename=paste0("Plots/BigDataFilter/Days_",i,"_Big_Network_Nplus_Phylum_Stab_Filtered_vsized.png"))
qgraph::qgraph(network_Nplus, 
  layout = "circle",
  edge.color = ifelse(network_Nplus > 0, "blue", "red"),
  title = "Stab Network Nplus by Phylum",
  vsize = 2.5,
  groups = Phylum_groups)
dev.off()
# %%
png(filename=paste0("Plots/BigDataFilter/Days_",i,"_Big_Network_Nminus_Phylum_Stab_Filtered_vsized.png"))
qgraph::qgraph(network_Nminus, 
  layout = "circle",
  edge.color = ifelse(network_Nminus > 0, "blue", "red"),
  title = "Stab Network Nminus by Phylum",
  vsize = 2.5,
  groups = Phylum_groups)
dev.off()

# Visualize network_Nplus (Circular Layout)
otu_tax_df <- tax_table(rawdata)[,1:5] %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
    dplyr::select(-OTU, everything(), OTU)

pairs <- list(
  c("Kingdom", "Phylum"),
  c("Phylum", "Class"),
  c("Class", "Order"),
  c("Order", "Family"),
  c("Family", "OTU"))
# Function to create edges for a single pair
create_edges <- function(pair, data)
  {
    from <- data[[pair[1]]]
    to <- data[[pair[2]]]
    data.frame(from = from, to = to)
  }

# Apply the function to all pairs and combine the results
edges <- unique(do.call(rbind, lapply(pairs, create_edges, data = otu_tax_df[otu_tax_df$OTU %in% shared_otu, ])))

# Extract lower triangular part
lower_tri <- lower.tri(network_Nplus, diag = FALSE)
# Get non-zero elements and their indices
non_zero <- which(lower_tri & network_Nplus != 0, arr.ind = TRUE)
# Create the new table
connect <- data.frame(
  from = rownames(network_Nplus)[non_zero[, 1]],
  to = colnames(network_Nplus)[non_zero[, 2]],
  score = network_Nplus[non_zero])

# create a vertices data.frame. One line per object of our hierarchy
vertices  <-  data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to))) , 
  value = runif(length(unique(c(as.character(edges$from), as.character(edges$to))))))
vertices$group  <-  edges$from[ match( vertices$name, edges$to ) ]

vertices$id <- NA
myleaves <- which(is.na(match(vertices$name, edges$from)))
vertices$value[myleaves] <- colSums(network_Nplus) / sum(network_Nplus)
nleaves <- length(myleaves)
vertices$id[myleaves] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)
mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)

ggraph::ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, width=0.1, aes(colour = after_stat(index))) +
  scale_edge_colour_gradient(low = "red", high = "blue") +
  geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=name, angle = angle, hjust=hjust, colour=group), size=2, alpha=1) +
  geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=value, alpha=0.2)) +
  scale_colour_manual(values= rep( brewer.pal(14,"Paired") , 30)) +
  scale_size_continuous(range = c(0.1,10) ) +
# set width by edge cor value?
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
ggsave(filename=paste0("Plots/BigDataFilter/Days_",i,"_Big_Nplus_Filtered_plot_circularized.pdf"), width = 12, height = 12, units = "in")

# %%Visualize Edge weights
cor_values_Nplus <- as.vector(network_Nplus)
cor_values_Nminus <- as.vector(network_Nminus)
cor_df <- data.frame(
  correlation = c(cor_values_Nplus, cor_values_Nminus),
  group = factor(rep(c("Nplus", "Nminus"), each = length(cor_values_Nplus)))
)
ggplot2::ggplot(cor_df, aes(x = correlation, fill = group)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  scale_fill_manual(values = c("Nplus" = "blue", "Nminus" = "red")) +
  theme_minimal() +
  labs(title = "Distribution of Correlation Values",
      x = "Correlation",
      y = "Frequency",
      fill = "Nitro Condition")

#%%
set.seed(10010)
synthesize_scaled_data <- function(dat, net)
  {
    graph <- (net != 0)*1     # SpiecEasi::make_graph('cluster', dat$n_OTUs, dat$n_edge)
    attr(graph, "class") <- "graph"
    Prec <- SpiecEasi::graph2prec(graph)
    Cor <- cov2cor(SpiecEasi::prec2cov(Prec))
    X <- SpiecEasi::synth_comm_from_counts(dat, mar = 2, distr = 'zinegbin', Sigma = Cor, n = nrow(dat))
    return(X)
  }
cat('Synthesize simulation data on day ',i,'\n')
Sim_list <- list()
for (j in 1:5)
  {
    Sim_list[[j]] <- list(
      Nplus = synthesize_scaled_data(otu_Ab_Nplus, network_Nplus),
      Nminus = synthesize_scaled_data(otu_Ab_Nminus, network_Nminus)
    )
  }
cat('Calculate simulation data network on day ',i,'\n')
Res_sim <- list()
for (j in 1:5)
  {
    Res_sim[[j]] <- stabENG(Sim_list[[j]], labels = shared_otu, var.thresh = 0.1, rep.num = 25,
      nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,
      lambda2.init=0.01,ebic.gamma=0.6)
    # filter edge sparsity
    Res_sim[[j]]$opt.fit$Nplus[abs(Res_sim[[j]]$opt.fit$Nplus) < 0.01] <- 0
    Res_sim[[j]]$opt.fit$Nminus[abs(Res_sim[[j]]$opt.fit$Nminus) < 0.01] <- 0
    diag(Res_sim[[j]]$opt.fit$Nplus) = diag(Res_sim[[j]]$opt.fit$Nminus) <- 0
  }
Sim_adj <- list()
for (j in 1:5)
  {
    Sim_adj[[j]] <- list(
      Nplus = (Res_sim[[j]]$opt.fit$Nplus !=0)*1,
      Nminus = (Res_sim[[j]]$opt.fit$Nminus !=0)*1
    )
  }
#%% Confusion matrices
calculate_mcc <- function(tp, tn, fp, fn)
  {
    numerator <- (tp * tn) - (fp * fn)
    denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

    # If the denominator is 0, return 0 to avoid NaN errors
    if (denominator == 0) {
      return(0)
    } else {
      return(numerator / denominator)
    }
  }
calculate_metrics <- function(true_adj, sim_adj)
  {
    # Flatten the matrices into vectors (upper triangular part, excluding diagonal)
    true_edges <- as.vector(true_adj[upper.tri(true_adj)])
    sim_edges <- as.vector(sim_adj[upper.tri(sim_adj)])
  
    # Confusion matrix
    cm <- confusionMatrix(as.factor(sim_edges), as.factor(true_edges), positive = "1")
    tn <- as.numeric(cm$table[1,1]) #true negatives
    fn <- as.numeric(cm$table[1,2]) #false negatives
    fp <- as.numeric(cm$table[2,1]) #false positives
    tp <- as.numeric(cm$table[2,2]) #true positives
  
    # Calculate TPR, FPR, Precision, Recall
    tpr <- tp / (tp + fn)  # Sensitivity / Recall
    fpr <- fp / (fp + tn)  # 1 - Specificity
    precision <- tp / (tp + fp)
    recall <- tpr
  
    # Calculate ROC and AUC
    roc_obj <- pROC::roc(true_edges, sim_edges)
    auc <- as.numeric(auc(roc_obj)) # plot here
  
    # F1 Score
    f1 <- F1_Score(sim_edges, true_edges)
  
    # MCC (Matthews correlation coefficient)
    mcc <- calculate_mcc(tp, tn, fp, fn)
  
    # Return metrics as a list
    return(list(TPR = tpr, FPR = fpr, Precision = precision, Recall = recall,
                F1 = f1, AUC = auc, MCC = mcc))
  }
cat('Calculate confusion matrices on day ',i,'\n')
confusion_results <- lapply(1:5, function(j)
  {
    true_adj_Nplus <- (network_Nplus !=0)*1
    true_adj_Nminus <- (network_Nminus !=0)*1
    Nplus_metrics <- calculate_metrics(true_adj_Nplus, Sim_adj[[j]]$Nplus)
    Nminus_metrics <- calculate_metrics(true_adj_Nminus, Sim_adj[[j]]$Nminus)
    return(list(Nplus = Nplus_metrics, Nminus = Nminus_metrics))
})


results_df <- do.call(rbind, lapply(confusion_results, as.data.frame))

results_df_long <- results_df %>%
  dplyr::select(starts_with("Nplus.") | starts_with("Nminus.")) %>%
  tidyr::pivot_longer(cols = everything(), 
              names_to = c("group", "metric"), 
              names_sep = "\\.",
              values_to = "value") %>%
  dplyr::mutate(matrix_id = rep(1:5, each = 14))

# %% barplot
for (metric_name in unique(results_df_long$metric)) {
  plot_data <- results_df_long %>%
    dplyr::filter(metric == metric_name)
  p <- ggplot(plot_data, aes(x = matrix_id, y = value, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste(metric_name, "for Simulated Networks"),
      x = "Matrix ID",
      y = metric_name) +
    theme_bw()
    ggsave(filename = paste0("Plots/BigDataFilter/Days_Big_Filtered_",i,"_", metric_name, "_barplot.png"), p)
}

# %% boxplot
for (metric_name in unique(results_df_long$metric)) {
  plot_data <- results_df_long %>% filter(metric == metric_name)

  p <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(title = paste(metric_name, "Distribution"),
        x = "Nitrogen Condition",
        y = metric_name) +
    theme_bw()
  ggsave(filename = paste0("Plots/BigDataFilter/Days_Big_Filtered_",i,"_", metric_name, "_boxplot.png"), p)
}
