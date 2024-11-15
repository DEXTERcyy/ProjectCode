myNetAnalysis <- function(data, data2,sampleSize = c(118,127), 
  alpha = c(0.05,0.05), lfdrThresh = c(0.2,0.2), sparsMethod = "t-test")
  {
    set.seed(10010)
    assoMat1 <- data
    assoMat2 <- data2
    counts1 <- counts2 <- NULL
    countsJointOrig <- NULL
    countsOrig1 <- countsOrig2 <- NULL
    groups <- NULL
    dissFunc = "signed"
    dissFuncPar = NULL
    simFunc = simFuncPar = NULL
    weighted = TRUE
    adjust = "none" #"lfdr"
  
  
    sparsReslt <- .sparsify(assoMat = assoMat1,
      countMat = NULL,
      sampleSize = sampleSize[1],
      measure = "spieceasi",
      measurePar = NULL,
      assoType = "partialCorr",
      sparsMethod = sparsMethod,
      thresh = NULL,
      alpha = alpha[1],
      adjust = adjust,
      lfdrThresh = lfdrThresh[1],
      trueNullMethod = NULL,
      nboot = NULL,
      assoBoot = FALSE,
      cores = 18,
      softThreshType = "unsigned",
      softThreshPower = NULL,
      softThreshCut = NULL,
      logFile = NULL,
      kNeighbor = 3L,
      knnMutual = FALSE,
      verbose = FALSE,
      seed = 10010)

    assoEst1 <- assoMat1
    assoMat1 <- sparsReslt$assoNew
    power1 <- sparsReslt$power
    dissEst1 <- dissScale1 <- NULL

    dissMat1 <- .transToDiss(x = assoMat1,  dissFunc= dissFunc,
            dissFuncPar = dissFuncPar)

    if (sparsMethod == "softThreshold") {
      simMat1 <- sparsReslt$simMat
      adjaMat1 <- assoMat1
      assoBoot1 <- NULL

    } else {
      simMat1 <- .transToSim(x = dissMat1, simFunc = simFunc,
              simFuncPar = simFuncPar)

      adjaMat1 <- .transToAdja(x = simMat1, weighted = weighted)

      if (sparsMethod == "bootstrap" && !is.null(assoBoot) &&
        !is.list(assoBoot) && assoBoot == TRUE) {
        assoBoot1 <- sparsReslt$assoBoot
      } else {
      assoBoot1 <- NULL
      }
    }

    sparsReslt <- .sparsify(assoMat = assoMat2,
            countMat = NULL,
            sampleSize = sampleSize[2],
            measure = "spieceasi",
            measurePar = NULL,
            assoType = "partialCorr",
            sparsMethod = sparsMethod,
            thresh = NULL,
            alpha = alpha[2],
            adjust = adjust,
            lfdrThresh = lfdrThresh[2],
            trueNullMethod = NULL,
            nboot = NULL,
            assoBoot = FALSE,
            cores = 18,
            softThreshType = "unsigned",
            softThreshPower = NULL,
            softThreshCut = NULL,
            logFile = NULL,
            kNeighbor = 3L,
            knnMutual = FALSE,
            verbose = FALSE,
            seed = 10010)

    assoEst2 <- assoMat2
    assoMat2 <- sparsReslt$assoNew
    power2 <- sparsReslt$power
    dissEst2 <- dissScale2 <- NULL

    dissMat2 <- .transToDiss(x = assoMat2, dissFunc = dissFunc,
              dissFuncPar = dissFuncPar)

    if (sparsMethod == "softThreshold") {
      simMat2 <- sparsReslt$simMat
      adjaMat2 <- assoMat2
      assoBoot2 <- NULL
    } else {
      simMat2 <- .transToSim(x = dissMat2, simFunc = simFunc,
              simFuncPar = simFuncPar)

      adjaMat2 <- .transToAdja(x = simMat2, weighted = weighted)

      if (sparsMethod == "bootstrap" && !is.null(assoBoot) &&
        !is.list(assoBoot) && assoBoot == TRUE) {
        assoBoot2 <- sparsReslt$assoBoot}
      else {
        assoBoot2 <- NULL
      }
    }
    # Create edge list
    g <- igraph::graph_from_adjacency_matrix(adjaMat1, weighted = weighted,
      mode = "undirected", diag = FALSE)

    if (is.null(igraph::E(g)$weight)) {
      isempty1 <- TRUE
      edgelist1 <- NULL
    } else {
      isempty1 <- FALSE

      edgelist1 <- data.frame(igraph::get.edgelist(g))
      colnames(edgelist1) <- c("v1", "v2")
      if (!is.null(assoMat1)) {
        edgelist1$asso <- sapply(1:nrow(edgelist1), function(i) {
          assoMat1[edgelist1[i, 1], edgelist1[i, 2]]
        })
      }

      edgelist1$diss <- sapply(1:nrow(edgelist1), function(i) {
        dissMat1[edgelist1[i, 1], edgelist1[i, 2]]
      })

      if (all(adjaMat1 %in% c(0,1))) {
        edgelist1$sim <-sapply(1:nrow(edgelist1), function(i) {
          simMat1[edgelist1[i, 1], edgelist1[i, 2]]
        })
      }

      edgelist1$adja <- sapply(1:nrow(edgelist1), function(i) {
        adjaMat1[edgelist1[i, 1], edgelist1[i, 2]]
      })
    }

    # Create edge list
    g <- igraph::graph_from_adjacency_matrix(adjaMat2, weighted = weighted, 
          mode = "undirected", diag = FALSE)
  
    if (is.null(igraph::E(g)$weight)) {
      isempty2 <- TRUE
      edgelist2 <- NULL
    } else {
      isempty2 <- FALSE

      edgelist2 <- data.frame(igraph::get.edgelist(g))
      colnames(edgelist2) <- c("v1", "v2")

      if (!is.null(assoMat2)) {
        edgelist2$asso <- sapply(1:nrow(edgelist2), function(i) {
          assoMat2[edgelist2[i, 1], edgelist2[i, 2]]
        })
      }

      edgelist2$diss <- sapply(1:nrow(edgelist2), function(i) {
        dissMat2[edgelist2[i, 1], edgelist2[i, 2]]
      })

      if (all(adjaMat2 %in% c(0, 1))) {
        edgelist2$sim <- sapply(1:nrow(edgelist2), function(i) {
          simMat2[edgelist2[i, 1], edgelist2[i, 2]]
        })
      }

      edgelist2$adja <- sapply(1:nrow(edgelist2), function(i) {
        adjaMat2[edgelist2[i, 1], edgelist2[i, 2]]
      })
    }

    if (isempty1 && verbose > 0) {
      message("\nNetwork 1 has no edges.")
    }
    if (isempty2 && verbose > 0) {
      message("Network 2 has no edges.")
    }

    #=============================================================================
    output <- list()

    output$edgelist1 <- edgelist1
    output$edgelist2 <- edgelist2
    output$assoMat1 <- assoMat1
    output$assoMat2 <- assoMat2
    output$dissMat1 <- dissMat1
    output$dissMat2 <- dissMat2
    output$simMat1 <- simMat1
    output$simMat2 <- simMat2
    output$adjaMat1 <- adjaMat1
    output$adjaMat2 <- adjaMat2

    output$assoEst1 <- assoEst1
    output$assoEst2 <- assoEst2
    output$dissEst1 <- dissEst1
    output$dissEst2 <- dissEst2
    output$dissScale1 <- dissScale1
    output$dissScale2 <- dissScale2

    output$assoBoot1 <- assoBoot1
    output$assoBoot2 <- assoBoot2

    output$countMat1 <- countsOrig1
    output$countMat2 <- countsOrig2
    if (!is.null(countsJointOrig)) output$countsJoint <- countsJointOrig
    output$normCounts1 <- counts1
    output$normCounts2 <- counts2
    output$groups <- groups
    output$sampleSize <- sampleSize
    output$softThreshPower <- list(power1 = power1, 
      power2 = power2) # calculated power

    output$call = match.call()
    output$twoNets <- TRUE
    class(output) <- "microNet"
    return(output)
}
# %%
# load(file = "DataImage\\network_results_stabENG.RData")
source(file="C:\\Users\\dexter\\Documents\\VU_UVA\\ResearchProject\\ProjectCode\\Packages\\NetComi\\R\\dot-sparsify.R")
source(file="C:\\Userss\\dexter\\Downloads\\VU_UVA\\ResearchProject\\ProjectCode\\Packages\\NetComi\\R\\transform.R")
# %%
NetRes <- myNetAnalysis(network_Nplus_pcor,network_Nminus_pcor)
netAnalyze(NetRes,
  clustMethod = "cluster_fast_greedy")