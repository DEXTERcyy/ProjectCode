# %% Install necessary packages
load(file = "DataImage\\big1226_Days_network_results_Big_Days_Filtered.RData")
if(!requireNamespace("discordant", quietly = TRUE)) install.packages("discordant")
if(!requireNamespace("Biobase", quietly = TRUE)) install.packages("Biobase")
if(!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
# %% Modified function to analyze differential correlations with better handling for small feature sets
analyzeFisherCorrelations <- function(assoMat1, assoMat2, countMat1, countMat2,
  output_dir = output_dir,
  p_adjust_method = "BH", # Benjamini-Hochberg FDR control
  alpha = 0.1,           # Significance level for adjusted p-values
  diff_threshold = 0.001) {# Optional minimum absolute difference in correlation

    # Check names consistency
    feature_names <- colnames(assoMat1)
    # Get sample sizes
    n1 <- nrow(countMat1)
    n2 <- nrow(countMat2)

    if (n1 <= 3 || n2 <= 3) {
    stop("Sample sizes (rows in countMat1/countMat2) must be greater than 3 for Fisher's z-test.")
    }

    cat("Comparing networks with", length(feature_names), "features.\n")
    cat("Sample sizes: n1 =", n1, ", n2 =", n2, "\n")

    # --- Perform Pairwise Comparison using Fisher's z-test ---

    # Get indices for the upper triangle (excluding the diagonal)
    upper_tri_indices <- which(upper.tri(assoMat1, diag = FALSE), arr.ind = TRUE)
    num_tests <- nrow(upper_tri_indices)

    # Initialize results storage
    results_list <- vector("list", num_tests)

    cat("Performing", num_tests, "pairwise correlation comparisons...\n")

    # Pre-calculate standard error denominator part (constant for all tests)
    se_denom_part <- sqrt(1/(n1 - 3) + 1/(n2 - 3))

    # Clamp function to keep correlations slightly away from +/- 1 for z-transform
    clamp_r <- function(r) {
    # Check if r is NA or NaN before clamping
    r_clamped <- ifelse(is.na(r) | is.nan(r), NA, pmin(pmax(r, -0.99999), 0.99999))
    return(r_clamped)
    }

    # Fisher's z-transform function
    fisher_z <- function(r) {
    0.5 * log((1 + r) / (1 - r))
    }

    for (k in 1:num_tests) {
    i <- upper_tri_indices[k, 1]
    j <- upper_tri_indices[k, 2]

    feat1_name <- feature_names[i]
    feat2_name <- feature_names[j]

    r1 <- assoMat1[i, j]
    r2 <- assoMat2[i, j]

    # Handle potential NA values in correlations
    if (is.na(r1) || is.na(r2)) {
    z1 <- NA
    z2 <- NA
    z_stat <- NA
    p_value <- NA
    } else {
    # Clamp correlations away from exactly +/- 1
    r1_clamped <- clamp_r(r1)
    r2_clamped <- clamp_r(r2)

    # Apply Fisher's z-transformation
    z1 <- fisher_z(r1_clamped)
    z2 <- fisher_z(r2_clamped)

    # Calculate the Z statistic for the difference
    z_stat <- (z1 - z2) / se_denom_part

    # Calculate the two-sided p-value from the standard normal distribution
    p_value <- 2 * pnorm(-abs(z_stat))
    }

    results_list[[k]] <- data.frame(
    Feature1 = feat1_name,
    Feature2 = feat2_name,
    Corr1 = r1,
    Corr2 = r2,
    Diff = r1 - r2,
    Z_stat = z_stat,
    p_value = p_value,
    stringsAsFactors = FALSE
    )
    }

    # Combine results into a single data frame
    results_df <- do.call(rbind, results_list)

    # Remove rows where p-value calculation failed (due to NAs in correlations)
    results_df <- na.omit(results_df) 

    if(nrow(results_df) == 0) {
    warning("No valid correlation pairs found or all resulted in NA p-values. Returning empty data frame.")
    return(data.frame())
    }

    cat("Performed comparisons. Adjusting p-values using method:", p_adjust_method, "...\n")

    # --- Adjust P-values for Multiple Testing ---
    results_df$p_adjusted <- p.adjust(results_df$p_value, method = p_adjust_method)

    # --- Filter Significant Results ---
    significant_results <- subset(results_df, 
    p_adjusted < alpha & abs(Diff) >= diff_threshold)

    # Sort by adjusted p-value
    significant_results <- significant_results[order(significant_results$p_adjusted), ]

    cat("Found", nrow(significant_results), "significant differential connections (edges) at alpha =", alpha, 
    "and difference threshold =", diff_threshold, "\n")

    # Save Results to File
    if(nrow(significant_results) > 0) {
      output_file <- file.path(output_dir, "significant_differential_connections.csv")
      write.csv(significant_results, file = output_file, row.names = FALSE)
      cat("Significant differential connections saved to:", output_file, "\n")
    } else {
      warning("No significant differential connections found based on the specified criteria.\n")
    }
    # --- Return Results ---
    return(significant_results)
}
# %%
analyzeDiscordantCorrelations <- function(assoMat1, assoMat2, countMat1, countMat2,
  discordThresh = 0.6, plot = TRUE, 
  output_dir = "DiffNetResults//Discordant//Timeseries",
  ncomp = NULL, # Number of components (auto-calculated if NULL)
  seed = 123) { # Added seed for reproducibility

# --- Input Validation and Setup ---

# Ensure necessary packages are loaded or checked
if (!requireNamespace("Biobase", quietly = TRUE)) {
stop("Package 'Biobase' is required. Please install it.", call. = FALSE)
}
# Load fully as we use ExpressionSet etc.
library(Biobase) 

# Check for discordant 
if (!requireNamespace("discordant", quietly = TRUE)) {
stop("Package 'discordant' is required but not installed.", call. = FALSE)
}

# Check for plotting package if needed
RColorBrewer_available <- requireNamespace("RColorBrewer", quietly = TRUE)
if(plot && !RColorBrewer_available) {
warning("Package 'RColorBrewer' not found. Plotting will use basic colors.", call. = FALSE)
}

# Check input matrix types and basic structure
if (!is.matrix(assoMat1) || !is.matrix(assoMat2) || !is.numeric(assoMat1) || !is.numeric(assoMat2)) {
stop("assoMat1 and assoMat2 must be numeric matrices.")
}
if (nrow(assoMat1) != ncol(assoMat1) || nrow(assoMat2) != ncol(assoMat2)) {
stop("Association matrices (assoMat1, assoMat2) must be square.")
}
if (!is.matrix(countMat1) || !is.matrix(countMat2)) {
stop("countMat1 and countMat2 must be matrices (samples x features).")
}

# Check dimensions consistency
n_features <- ncol(countMat1)
if (n_features != ncol(countMat2)) {
stop("Count matrices must have the same number of features (columns).")
}
if (nrow(assoMat1) != n_features || ncol(assoMat1) != n_features || 
nrow(assoMat2) != n_features || ncol(assoMat2) != n_features) {
stop("Association matrix dimensions must match each other and the number of features in count matrices.")
}

# Check names consistency (CRITICAL)
# Assume colnames of countMat hold the canonical feature names
feature_names <- colnames(countMat1) 
if (is.null(feature_names)) {
stop("countMat1 must have column names (feature names).")
}
if (!identical(feature_names, colnames(countMat2))) {
stop("Column names of countMat1 and countMat2 must be identical.")
}
# Check assoMat names against feature_names
if (is.null(rownames(assoMat1)) || is.null(colnames(assoMat1)) || 
!identical(feature_names, rownames(assoMat1)) || !identical(feature_names, colnames(assoMat1)) ||
is.null(rownames(assoMat2)) || is.null(colnames(assoMat2)) ||
!identical(feature_names, rownames(assoMat2)) || !identical(feature_names, colnames(assoMat2))) {
stop("Row and column names of assoMat1 and assoMat2 must exist and be identical ",
"to each other and to the column names of the count matrices.")
}

# Create output directory
if(plot && !dir.exists(output_dir)) {
dir.create(output_dir, recursive = TRUE)
cat("Created output directory:", output_dir, "\n")
}

cat("Processing", n_features, "features...\n")

# --- 1. Prepare data for discordant ---

# Combine count data
combMat <- rbind(countMat1, countMat2)

# Assign unique sample names (important for ExpressionSet)
n_samples1 <- nrow(countMat1)
n_samples2 <- nrow(countMat2)
sample_names <- make.unique(c(rownames(countMat1), rownames(countMat2)))
if (length(sample_names) != (n_samples1 + n_samples2)) {
warning("Generating generic unique sample names as rownames were not unique/present.")
sample_names <- paste0("Sample_", 1:(n_samples1 + n_samples2))
}
rownames(combMat) <- sample_names
colnames(combMat) <- feature_names # Ensure colnames are correct

# Create ExpressionSet object (features as rows in assayData)
x_expr <- Biobase::ExpressionSet(assayData = t(combMat)) 
# Assign feature and sample names explicitly
Biobase::featureNames(x_expr) <- feature_names
Biobase::sampleNames(x_expr) <- sample_names

# Define groups and add to phenoData (optional but good practice)
groups <- factor(c(rep("Group1", n_samples1), rep("Group2", n_samples2)))
Biobase::pData(x_expr) <- data.frame(Group = groups, row.names = sample_names)

# --- 2. Extract and name correlation vectors ---

# Use upper triangle to avoid redundant pairs
upper_tri_idx <- which(upper.tri(assoMat1, diag = FALSE))
corrVector1 <- assoMat1[upper_tri_idx]
corrVector2 <- assoMat2[upper_tri_idx]

# Optimized function to get vector names using matrix indexing
.getVecNames <- function(feature_names) {
n <- length(feature_names)
# Get row and column indices for the upper triangle
idx <- which(upper.tri(matrix(NA, n, n), diag = FALSE), arr.ind = TRUE)
# Paste names based on indices (col index corresponds to outer loop in original code)
paste(feature_names[idx[, 1]], feature_names[idx[, 2]], sep = "_")
}

# Add names to correlation vectors
vector_names <- .getVecNames(feature_names)
if (length(corrVector1) != length(vector_names)) {
stop("Internal error: Length mismatch between correlation vectors and generated names.")
}
names(corrVector1) <- vector_names
names(corrVector2) <- vector_names

# --- 3. Run discordant analysis ---

# Initialize output matrices with NAs initially, named correctly
classMat <- matrix(NA_integer_, nrow=n_features, ncol=n_features, dimnames = list(feature_names, feature_names))
diffProbs <- matrix(NA_real_, nrow=n_features, ncol=n_features, dimnames = list(feature_names, feature_names))
diag(classMat) <- 1 # Self-comparison is always class 1 (same)
diag(diffProbs) <- 0 # Probability of self being different is 0

# Auto-calculate number of components if not provided
ncomp_calculated <- FALSE
if(is.null(ncomp)) {
ncomp_calculated <- TRUE
# Refined Heuristic: Consider samples and features
n_samples_total = n_samples1 + n_samples2
# Ensure ncomp is less than the number of features and effective samples/groups
max_comp_feat = max(2, floor(n_features / 5)) # Need enough features per component? Maybe not strict.
max_comp_samp = max(2, floor(n_samples_total / 3)) # Need enough samples per component
# PCA/SVD constraint: ncomp should be <= min(n_features, n_samples_total) - 1
# Often safer to be less than smaller group size - 1
max_comp_group = min(n_samples1, n_samples2) - 1

# Proposed ncomp: Min of sample/group constraints, capped reasonably
ncomp <- min(max_comp_samp, max_comp_group, 5) 
ncomp <- max(ncomp, 2) # Ensure at least 2 components

cat("Automatically calculated ncomp =", ncomp, "\n")

# Final sanity check for very small datasets
if(ncomp >= min(n_features, n_samples_total)) {
original_ncomp = ncomp
ncomp = max(2, min(n_features, n_samples_total) - 1) # Adjust down
cat("Adjusted ncomp to", ncomp, "due to small data dimensions (Features:", n_features, ", Samples:", n_samples_total,"). Original calculation was", original_ncomp, ".\n")
}
if (ncomp < 2) {
stop("Cannot determine a valid ncomp >= 2 based on data dimensions. Need more samples or features, or specify ncomp manually.")
}
}
ncomp_requested <- ncomp # Store the ncomp value we intend to use

# Run discordant analysis
cat("Running discordant analysis with ncomp =", ncomp_requested, "...\n")
set.seed(seed) # Set seed *before* the potentially stochastic analysis
discord <- discordant::discordantRun(corrVector1, corrVector2, x_expr, components=ncomp_requested)

# Post-run validation
if(is.null(discord) || !is.list(discord) || 
is.null(discord$classMatrix) || is.null(discord$discordPPMatrix) ||
!is.matrix(discord$classMatrix) || !is.matrix(discord$discordPPMatrix) ||
nrow(discord$classMatrix) != n_features || ncol(discord$discordPPMatrix) != n_features) {
stop("discordantRun completed but returned an invalid/incomplete object structure.")
}

# Get results and fill matrices (use upper triangle assignment for efficiency)
classMat[upper.tri(classMat)] <- discord$classMatrix[upper.tri(discord$classMatrix)]
diffProbs[upper.tri(diffProbs)] <- discord$discordPPMatrix[upper.tri(discord$discordPPMatrix)]

# Ensure symmetry
classMat[lower.tri(classMat)] <- t(classMat)[lower.tri(t(classMat))]
diffProbs[lower.tri(diffProbs)] <- t(diffProbs)[lower.tri(t(diffProbs))]

cat("Discordant analysis successful.\n")

# --- 5. Plotting (Optional) ---
if(plot) {
cat("Generating plot...\n")
plot_file <- file.path(output_dir, "correlation_comparison_plot.png")
png(plot_file, width=800, height=800, res=100)

# Initialize plot variables
plot_colors <- "black" # Default if something goes wrong
main_title <- "Correlation Comparison"
legend_elements <- list(labels = NULL, colors = NULL, pch = NULL)

tryCatch({ # Wrap plotting in tryCatch to avoid stopping the function if plotting fails
# Get pair classes for coloring
if(!is.null(discord$classMatrix)) {
# Extract from matrix to build pair classes vector
pairClasses <- classMat[upper.tri(classMat)]

# Define colors
if (RColorBrewer_available) {
base_cols <- RColorBrewer::brewer.pal(max(3, max(pairClasses, na.rm=TRUE)), "Set1") 
color_map <- c("1" = "grey80", "2" = base_cols[2], "3" = base_cols[1]) # 1=grey, 2=Blue, 3=Red
} else {
color_map <- c("1" = "grey80", "2" = "blue", "3" = "red") # Basic fallback colors
}

plot_colors <- color_map[as.character(pairClasses)]
plot_colors[is.na(plot_colors)] <- "black" # Color NAs black

main_title <- paste("Correlation Comparison (Discordant Method, ncomp=", ncomp_requested, ")", sep="")
potential_labels <- c("Same (Class 1)", "Diff. Magnitude (Class 2)", "Diff. Sign (Class 3)")
potential_colors <- color_map[c("1", "2", "3")]

# Only include legend entries for classes present
present_classes <- unique(na.omit(pairClasses))
legend_mask <- c(1 %in% present_classes, 2 %in% present_classes, 3 %in% present_classes)

if(any(legend_mask)) {
legend_elements$labels <- potential_labels[legend_mask]
legend_elements$colors <- potential_colors[legend_mask]
legend_elements$pch <- rep(19, sum(legend_mask))
}
} else {
warning("Discordant object structure unexpected for plotting, using default colors.", call.=FALSE)
plot_colors <- "black"
main_title <- "Correlation Comparison (Discordant Analysis)"
}

# Create the scatter plot
plot(corrVector1, corrVector2, 
xlab = "Correlation (N-plus)", 
ylab = "Correlation (N-minus)", 
main = main_title,
col = plot_colors, 
pch = 19, 
cex = 0.6,
xlim = c(-1, 1), ylim = c(-1, 1))
abline(h = 0, v = 0, lty = 2, col = "grey50") # Zero lines
abline(a = 0, b = 1, lty = 3, col = "blue")   # y=x line (identity)

# Add legend if elements were defined
if(!is.null(legend_elements$labels)) {
legend("topleft", 
legend = legend_elements$labels, 
col = legend_elements$colors, 
pch = legend_elements$pch, 
bg = "white",
cex=0.8)
}

}, error = function(e_plot) {
cat("Error during plotting:", e_plot$message, "\n")
# Try to ensure device is closed if error occurred during plotting
if(dev.cur() != 1) dev.off() 
})

# Ensure graphics device is off, even if tryCatch didn't catch error but png() was called
if(exists("plot_file") && dev.cur() != 1 && names(dev.cur()) == "png") {
dev.off()
}
cat("Plot generation finished. File:", plot_file, "\n")
} # End if(plot)

# --- 6. Return Results ---
cat("Analysis complete.\n")

return(list(
classMatrix = classMat, 
discordanceProbabilityMatrix = diffProbs,
discordObject = discord, # Full object from discordant::discordantRun
ncomp_used = ncomp_requested,
settings = list(
discordThresh = discordThresh,
seed = seed, 
ncomp_requested = ncomp_requested, # The ncomp value used
ncomp_calculated = ncomp_calculated # Flag if ncomp was auto-calculated
) 
))
}
# %% 
diffnetResults <- list()
for (i in timestamps){
  output_dir_Discordant <- paste0("DiffNetResults//Discordant//Day_",i)
  output_dir_Fisher <- paste0("DiffNetResults//Fisher//Day_",i)
  
  if(!dir.exists(output_dir_Discordant)) dir.create(output_dir_Discordant, recursive = TRUE)
  if(!dir.exists(output_dir_Fisher)) dir.create(output_dir_Fisher, recursive = TRUE)
  

  countMat1 <- otu_Ab_Nplus_times[[i]]
  countMat2 <- otu_Ab_Nminus_times[[i]]
  assoMat1 <- network_pcor_raw[[i]]$Nplus
  assoMat2 <- network_pcor_raw[[i]]$Nminus
  
  diffnetResults[[i]] <- list()
  
  # Discordant
  diffnetResults[[i]]$Discordant <- analyzeDiscordantCorrelations(
    assoMat1, assoMat2, countMat1, countMat2,
    output_dir = output_dir_Discordant
  )
  
  # Fisher
  diffnetResults[[i]]$Fisher <- analyzeFisherCorrelations(
    assoMat1, assoMat2, countMat1, countMat2, 
    output_dir = output_dir_Fisher
  )
}
