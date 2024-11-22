# %%
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
library(devtools, quietly = TRUE)
library(BiocManager, quietly = TRUE)
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("SPRING", quietly = TRUE)) devtools::install_github("GraceYoon/SPRING")
if (!requireNamespace("SpiecEasi", quietly = TRUE)) devtools::install_github("zdk123/SpiecEasi")
if (!requireNamespace("stabJGL", quietly = TRUE)) devtools::install_github("camiling/stabJGL")
devtools::install_github("stefpeschel/NetCoMi", ref = "develop",repos = c("https://cloud.r-project.org/",BiocManager::repositories()))
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