# #%% requirements
# # !ebic.gamma: if too sparse change to 0.6
cat(R.version.string, "\nR Home:", R.home(), "\nR lib:", .libPaths(), "\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
library(devtools, quietly = TRUE)
library(BiocManager, quietly = TRUE)
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("SPRING", quietly = TRUE)) devtools::install_github("GraceYoon/SPRING")
if (!requireNamespace("SpiecEasi", quietly = TRUE)) devtools::install_github("zdk123/SpiecEasi")
if (!requireNamespace("stabJGL", quietly = TRUE)) devtools::install_github("camiling/stabJGL")
if (!requireNamespace("NetCoMi", quietly = TRUE)) devtools::install_github("stefpeschel/NetCoMi", ref = "develop",repos = c("https://cloud.r-project.org/",BiocManager::repositories()))
if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
if (!requireNamespace("MLmetrics", quietly = TRUE)) install.packages("MLmetrics")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr") # Likely installed with tidyverse
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse") # Installs a suite of packages
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel") # Usually base R, no need to install
if (!requireNamespace("mixedCCA", quietly = TRUE)) install.packages("mixedCCA")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("qgraph", quietly = TRUE)) install.packages("qgraph")
start_time <- Sys.time()
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
source("//home//14720078//ProjectCode//Packages//stabENG.r")
source("//home//14720078//ProjectCode//Packages//MyENG.r")
source("//home//14720078//ProjectCode//Packages//stabENG.r")
source("//home//14720078//ProjectCode//Packages//MyPlot.r")
cat("Packages loaded successfully\n")
start_time <- Sys.time()
rawdata <- readRDS("//home//14720078//ProjectCode//data//DavarData1_substrate_phyloseq_1226_final_filtered.rds")
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
timestamps <- as.character(sort(as.integer(levels(sam_info$Days))))
otu_Ab_Nplus_times <- list()
otu_Ab_Nminus_times <- list()
data_list_times <- list()
n_sim <- 25
for (i in timestamps)
  {
    otu_Ab_Nplus_times[[i]] <- otu_Ab[rownames(otu_Ab) %in% 
      rownames(sam_info[sam_info$growthCondition=="plusN" & sam_info$Days == i,]),]
    otu_Ab_Nminus_times[[i]] <- otu_Ab[rownames(otu_Ab) %in%
      rownames(sam_info[sam_info$growthCondition=="minusN" & sam_info$Days == i,]),]
    data_list_times[[i]] <- list(Nplus = otu_Ab_Nplus_times[[i]], Nminus = otu_Ab_Nminus_times[[i]])
}
# %%
network_list_raw <- list()
network_list <- list()
network_pcor_raw <- list()
network_pcor <- list()
confusion_results_df <- list()
for (i in timestamps)
{
  plot_path = paste0("//home//14720078//ProjectCode//Plots//BigDataDaysFilter//Day_",i)
  data_list <- data_list_times[[i]]
  cat('Calculating network on day ',i,'\n')
  network_results <- stabENG(data_list, labels = shared_otu, var.thresh = 0.1, rep.num = 25,
    nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,
    lambda2.init=0.01,ebic.gamma=0.2) #nCores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
  network_list_raw[[i]]$Nplus <- network_results$opt.fit$Nplus # precision matrix estimates
  network_list_raw[[i]]$Nminus <- network_results$opt.fit$Nminus # precision matrix estimates
  network_pcor_raw[[i]]$Nplus <- network_results$opt.fit.pcor$Nplus
  network_pcor_raw[[i]]$Nminus <- network_results$opt.fit.pcor$Nminus
  # filter edge sparsity for simulation
  network_list[[i]]$Nplus <- network_list_raw[[i]]$Nplus
  network_list[[i]]$Nplus[abs(network_list[[i]]$Nplus) < 0.01] <- 0
  network_list[[i]]$Nminus <- network_list_raw[[i]]$Nminus
  network_list[[i]]$Nminus[abs(network_list[[i]]$Nminus) < 0.1] <- 0
  diag(network_list[[i]]$Nplus) = diag(network_list[[i]]$Nminus) <- 0
  # filter network_pcor
  network_pcor[[i]]$Nplus <- network_pcor_raw[[i]]$Nplus
  network_pcor[[i]]$Nplus[abs(network_pcor[[i]]$Nplus) < 0.01] <- 0
  network_pcor[[i]]$Nminus <- network_pcor_raw[[i]]$Nminus
  network_pcor[[i]]$Nminus[abs(network_pcor[[i]]$Nminus) < 0.01] <- 0
  # network_pcor[[i]]$Nplus[abs(network_pcor[[i]]$Nplus) < 0.01] <- 0
  # network_pcor[[i]]$Nminus[abs(network_pcor[[i]]$Nminus) < 0.01] <- 0
  diag(network_pcor[[i]]$Nplus) = diag(network_pcor[[i]]$Nminus) <- 0
  # %% Plot network on Phylum level
  Phylum_groups <- as.factor(otu_tax[rownames(network_list[[i]]$Nplus),"Phylum"])
  png(filename=paste0(plot_path,"_network_Nplus_Phylum_Stab_Filtered_vsized.png"))
  qgraph::qgraph(network_list[[i]]$Nplus, 
    layout = "circle",
    edge.color = ifelse(network_list[[i]]$Nplus > 0, "blue", "red"),
    title = "Stab Network Nplus by Phylum",
    vsize = 2.5,
    groups = Phylum_groups)
  dev.off()
  
  png(filename=paste0(plot_path,"_network_Nminus_Phylum_Stab_Filtered_vsized.png"))
  qgraph::qgraph(network_list[[i]]$Nminus, 
    layout = "circle",
    edge.color = ifelse(network_list[[i]]$Nminus > 0, "blue", "red"),
    title = "Stab Network Nminus by Phylum",
    vsize = 2.5,
    groups = Phylum_groups)
  dev.off()
  # %%Visualize Edge weights (pCor or Prec?)
  cor_values_Nplus <- as.vector(network_pcor[[i]]$Nplus)
  cor_values_Nminus <- as.vector(network_pcor[[i]]$Nminus)
  cor_values_Nplus_filtered <- cor_values_Nplus[cor_values_Nplus != 0]
  cor_values_Nminus_filtered <- cor_values_Nminus[cor_values_Nminus != 0]

  ggplot(data.frame(Type = c(rep("N+", length(cor_values_Nplus_filtered)), rep("N-", length(cor_values_Nminus_filtered))), 
    Value = c(cor_values_Nplus_filtered, cor_values_Nminus_filtered)), 
    aes(x = Value, color = Type)) + 
    geom_density() + 
    labs(x = "Correlation Value", y = "Density")
  ggsave(filename=paste0(plot_path,"_correlation_distribution.png"))
  
  # Visualize network_list[[i]]$Nplus (Circular Layout)
  cir_plot(network_list[[i]]$Nplus,otu_tax, shared_otu, paste0(plot_path,"_network_Nplus_Circular"))
  cir_plot(network_list[[i]]$Nminus,otu_tax, shared_otu, paste0(plot_path,"_network_Nminus_Circular"))
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
  for (j in 1:n_sim)
    {
      Sim_list[[i]][[j]] <- list(
        Nplus = synthesize_scaled_data(otu_Ab_Nplus, network_list[[i]]$Nplus),
        Nminus = synthesize_scaled_data(otu_Ab_Nminus, network_list[[i]]$Nminus)
      )
    }
  cat('Calculate simulation data network on day ',i,' start.\n')
  Res_sim <- list()
  for (j in 1:n_sim)
    {
      Res_sim[[i]][[j]] <- stabENG(Sim_list[[i]][[j]], labels = shared_otu, var.thresh = 0.1, rep.num = 25,
        nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,
        lambda2.init=0.01,ebic.gamma=0.2) #nCores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
      cat('Simulation Progress on Day',i,':', 100*j/n_sim,'% done.\n')
      # filter edge sparsity
      Res_sim[[i]][[j]]$opt.fit$Nplus[abs(Res_sim[[i]][[j]]$opt.fit$Nplus) < 0.01] <- 0
      Res_sim[[i]][[j]]$opt.fit$Nminus[abs(Res_sim[[i]][[j]]$opt.fit$Nminus) < 0.01] <- 0
      diag(Res_sim[[i]][[j]]$opt.fit$Nplus) = diag(Res_sim[[i]][[j]]$opt.fit$Nminus) <- 0
      # plot sim network
      png(filename=paste0(plot_path,"_simulation_network_",j,"_Nplus_Phylum_Stab_Filtered_vsized.png"))
      qgraph::qgraph(Res_sim[[i]][[j]]$opt.fit$Nplus, 
        layout = "circle",
        edge.color = ifelse(Res_sim[[i]][[j]]$opt.fit$Nplus > 0, "blue", "red"),
        title = "Stab Network Nplus by Phylum",
        vsize = 2.5,
        groups = Phylum_groups)
      dev.off()
      
      png(filename=paste0(plot_path,"_simulation_network_",j,"_Nminus_Phylum_Stab_Filtered_vsized.png"))
      qgraph::qgraph(Res_sim[[i]][[j]]$opt.fit$Nminus, 
        layout = "circle",
        edge.color = ifelse(Res_sim[[i]][[j]]$opt.fit$Nminus > 0, "blue", "red"),
        title = "Stab Network Nminus by Phylum",
        vsize = 2.5,
        groups = Phylum_groups)
      dev.off()
  }
  Sim_adj <- list()
  for (j in 1:n_sim)
    {
      Sim_adj[[i]][[j]] <- list(
        Nplus = (Res_sim[[i]][[j]]$opt.fit$Nplus !=0)*1,
        Nminus = (Res_sim[[i]][[j]]$opt.fit$Nminus !=0)*1
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
      fp <- as.numeric(cm$table[1,2]) #false positives
      fn <- as.numeric(cm$table[2,1]) #false negatives
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
  confusion_results <- lapply(1:n_sim, function(j)
    {
      true_adj_Nplus <- (network_list[[i]]$Nplus !=0)*1
      true_adj_Nminus <- (network_list[[i]]$Nminus !=0)*1
      Nplus_metrics <- calculate_metrics(true_adj_Nplus, Sim_adj[[i]][[j]]$Nplus)
      Nminus_metrics <- calculate_metrics(true_adj_Nminus, Sim_adj[[i]][[j]]$Nminus)
      return(list(Nplus = Nplus_metrics, Nminus = Nminus_metrics))
  })
  
  
  results_df <- do.call(rbind, lapply(confusion_results, as.data.frame))
  
  results_df_long <- results_df %>%
    dplyr::select(starts_with("Nplus.") | starts_with("Nminus.")) %>%
    tidyr::pivot_longer(cols = everything(), 
                names_to = c("group", "metric"), 
                names_sep = "\\.",
                values_to = "value") %>%
    dplyr::mutate(matrix_id = rep(1:n_sim, each = 14)) %>%
    dplyr::mutate(times = rep(i, n_sim*14))
  # append results_df_long to confusion_results_df
  confusion_results_df <- rbind(confusion_results_df, results_df_long)
  # # %% boxplot
  # metrics_to_plot <- c("TPR", "FPR", "Precision", "Recall", "F1", "AUC", "MCC")
  # filtered_df <- results_df_long %>%
  #   filter(metric %in% metrics_to_plot)
  
  # # Create boxplots
  # p <- ggplot(filtered_df, aes(x = group, y = value, color = group)) +
  #   geom_boxplot() +
  #   facet_wrap(~ metric, nrow = 1, scales = "free_y") + # facet_wrap for one row
  #   labs(title = "Distribution of Confusion Metrics",
  #        x = "Group",
  #        y = "Value",
  #        color = "Group") +
  #   theme_minimal()
  # ggsave(filename = paste0(plot_path,"_Filtered_Confusion_boxplot.png"), p)
}

# Filter for TPR, FPR, and F1 metrics
selected_metrics <- c("TPR", "FPR", "F1")
filtered_data <- confusion_results_df %>% filter(metric %in% selected_metrics)

# Plot
png(filename = "//home//14720078//ProjectCode//Plots//BigDataDaysFilter//boxplot_confusion_matrix.png")
ggplot(filtered_data, aes(x = factor(as.integer(times)), y = value)) + # fill = group
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black", position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Nplus" = "#00BFC4", "Nminus" = "#F8766D")) +  # Custom colors
  facet_wrap(~metric, scales = "free_y") +  # Create separate plots for each metric
  labs(x = "Times", y = "Value", title = "Confusion Metrics") +
  ylim(0, 1) +  # Set y-axis limits
  theme_minimal(base_size = 14) +  # Clean theme with adjusted font size
  theme(
    legend.position = "top",  # Move legend to top
    panel.grid.major = element_line(color = "gray90"),  # Add grid lines
    panel.border = element_rect(color = "black", fill = NA),  # Add border
    strip.text = element_text(face = "bold", size = 12)  # Make facet labels bold
  )
dev.off()
# save network_list and network_pcor as csv
write.csv(network_pcor, "//home//14720078//ProjectCode//DataImage//network_pcor.csv")
save.image("//home//14720078//ProjectCode//DataImage//big1226_Days_network_results_Big_Days_Filtered.RData")
cat("All done: ", Sys.time() - start_time, "\n")