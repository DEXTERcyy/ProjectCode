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
source("Packages\\stabENG.r")
rawdata <- readRDS("data\\DavarData1_substrate_phyloseq_1226_final_filtered.RDS")
otu_raw <- otu_table(rawdata)
otu_RA <- transform_sample_counts(otu_raw, function(x) x / sum(x) )
otu_RA <- filter_taxa(otu_RA, function(x) mean(x) > 1e-3, TRUE)
shared_otu <- rownames(otu_RA)
otu_Ab <- as.data.frame(t(otu_raw[shared_otu,]))
sam_info <- as.data.frame(sample_data(rawdata))
otu_tax <- as.data.frame(tax_table(rawdata))

# split otu_Ab by meta
#timestamps <- unique(sam_info$Days)
otu_Ab_Nplus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="plusN" & sam_info$Days == 4,]),]
otu_Ab_Nminus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="minusN" & sam_info$Days == 4,]),]
data_list <- list(Nplus = otu_Ab_Nplus, Nminus = otu_Ab_Nminus)

# %%
network_results <- stabENG(data_list, labels = shared_otu, var.thresh = 0.1, rep.num = 20,
  nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,
  lambda2.init=0.01,ebic.gamma=0.2, parallelize=T)
network_Nplus <- network_results$opt.fit$Nplus # precision matrix estimates
network_Nminus <- network_results$opt.fit$Nminus # precision matrix estimates
diag(network_Nplus) = diag(network_Nminus) <- 0
network_list <- list(network = list(Nplus = network_Nplus, Nminus = network_Nminus))
# %%
GroupNetworkBoot(data_list = data_list, groupNetwork = network_list,
  nboots = 10, bootSeed = 1234, ncores = 16)

# %%
# network_Nplus[abs(network_Nplus) < 0.01] <- 0
# network_Nminus[abs(network_Nminus) < 0.01] <- 0
network_Nplus_pcor <- network_results$opt.fit.pcor$Nplus
network_Nminus_pcor <- network_results$opt.fit.pcor$Nminus

save.image("DataImage\\big1226_network_results_stabENG.RData")

# %% Plot network on Phylum level
load(file = "DataImage\\big1226_network_results_stabENG.RData")
Phylum_groups <- as.factor(otu_tax[rownames(network_Nplus),"Phylum"])
png(filename="Plots/BigData/big1226_network_Nplus_Phylum_Stab.png")
qgraph::qgraph(network_Nplus, 
  layout = "circle",
  edge.color = ifelse(network_Nplus > 0, "blue", "red"),
  title = "Stab Network Nplus by Phylum",
  groups = Phylum_groups)
dev.off()

png(filename="Plots/BigData/big1226_network_Nminus_Phylum_Stab.png")
qgraph::qgraph(network_Nminus, 
  layout = "circle",
  edge.color = ifelse(network_Nminus > 0, "blue", "red"),
  title = "Stab Network Nminus by Phylum",
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

  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
ggsave(filename="Plots/BigData/big1226_Nplus_plot_circularized.pdf", width = 12, height = 12, units = "in")

# Visualize Edge weights
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

