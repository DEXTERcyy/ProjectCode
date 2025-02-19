GroupNetworkBoot <- function(
  data_list, 
  groupNetwork, 
  nboots = 100, 
  bootSeed,
  ...)     {
  
  
  #----------------Random seed for replicability --------
  if(!missing(bootSeed)) set.seed(bootSeed)
  
  # ---------------Check if data.frames have names------
 
 # If no names are present, these are added and a warning is given
 if(is.null(names(data_list))) {
   names(data_list) = paste0("G", 1:length(data_list))
   warning("data.frames in list have no names associated with them, so are coerced to G1, G2, ..., Gn.")
 }
 
  # ---------------Save important variables-------------
  G = length(data_list) # number of samples/groups
  Gnames = names(data_list) # Get sample data.frame names
  Gn = data.frame(G = 1:G) # Data.frame including sample sizes
  Gn$n = sapply(data_list, nrow) 
  
  nvar = ncol(groupNetwork$network[[1]])
  edges = colnames(groupNetwork$network[[1]])
  
  tracker = round(seq(from = 0, to = 100, length.out = nboots), 1)
  
  
  start_time = Sys.time()
  
  
  # --------------- Check data -------------------------
  if(any(sapply(data_list, ncol) != nvar))
  {
     stop("Error: all datasets (data_list) must include the same variables of the original network (groupNetwork)")
  }

  # ---------------Define empty arguments---------------

  labels = colnames(groupNetwork$network[[1]])

  # ---------------Create bootstrapping list------------
  output <- list(data = data_list,
                 sample = groupNetwork,
                 boot = list())
  
  
  # ---------------Run bootstrapping--------------------
  for(i in 1:nboots)      {
     
     # Bootstrap samples
     boot_dat = list()
     
     for(j in 1:G)     {
        
        ## Save group
        boot_dat[[j]] <- data_list[[j]][sample(1:Gn$n[j], size = Gn$n[j], replace = TRUE),]
        
     }
     
     ## Name boot_dat data.frames
     names(boot_dat) <- Gnames
     
     # ---------------Run EstimateGroupNetwork-------
     boot_network <- 
        myENG(boot_dat, n = Gn$n, weights = "equal", penalty = "fused",
      lambda1 = lambda1, lambda2 = lambda2, maxiter=500, tol=1e-5, rho=1, truncate=1e-5)
  
     
     ## Save bootstrap network edges to output list
     output$boot[[length(output$boot) + 1]] <- boot_network$network
     names(output$boot) <- paste("b", 1:i, sep = "")
     
     
     ## Print tracker to show progress
     
     # Calculate remaining time
     cur_time <- Sys.time()
     time_completed <- as.numeric(cur_time - start_time)
     time_remaining <- round((time_completed / i) * (nboots - i))
     hours_remaining <- round(time_remaining / 60)
     minutes_remaining <- time_remaining - hours_remaining * 60
     
     # Print tracker
     cat(paste("\r", tracker[i], "% (~", hours_remaining, " Hours ", 
               minutes_remaining, " Minutes remaining)        ", sep = ""))
     
     
     # Remove temporary variables
     rm(boot_dat)
     rm(boot_network)
     rm(time_remaining)
     rm(time_completed)
     rm(hours_remaining)
     rm(minutes_remaining)
     rm(cur_time)
     
  }
  
  # Return Output
  return(output)
  
}