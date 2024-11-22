# %%
library(phyloseq)
library(SPRING)
library(SpiecEasi)
library(qgraph)
library(igraph)
library(caret)     # For confusionMatrix function
library(pROC)      # For ROC and AUC calculation
library(MLmetrics) # For F1 Score, MCC
library(ggplot2)
library(tidyr)
library(ggraph)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
source("Packages\\stabENG.r")
rawdata <- readRDS("data\\phyloseqDavarHoagIntersect1009.RDS")
otu_Ab <- as.data.frame(t(otu_table(rawdata)))
sam_info <- as.data.frame(sample_data(rawdata))
otu_tax <- as.data.frame(tax_table(rawdata))
shared_otu <- colnames(otu_Ab)
# split otu_Ab by meta
otu_Ab_Nplus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="plusN",]),]
otu_Ab_Nminus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="minusN",]),]
data_list <- list(Nplus = otu_Ab_Nplus, Nminus = otu_Ab_Nminus)

# %% 
source("Packages\\stabENG.r")
network_results <- stabENG(data_list, labels = shared_otu, var.thresh = 0.1, rep.num = 20,
  nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,
  lambda2.init=0.01,ebic.gamma=0.2)
network_Nplus <- network_results$opt.fit$Nplus # precision matrix estimates
network_Nminus <- network_results$opt.fit$Nminus # precision matrix estimates
diag(network_Nplus) = diag(network_Nminus) <- 0
network_Nplus_pcor <- network_results$opt.fit.pcor$Nplus
network_Nminus_pcor <- network_results$opt.fit.pcor$Nminus