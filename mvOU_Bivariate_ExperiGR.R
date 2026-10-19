library(tidyverse)
library(mvSLOUCH)
library(doParallel)
library(geiger)
library(phytools)
library(PCMBaseCpp)

#----------------------------------------------------------#
#                       0 Preparation                      #
#----------------------------------------------------------#

### 1. Load Data ---
# Data
load("data_bivariate_mvOU_experiGR.RData")
# Load tree
load("tree_mvSLOUCH.RData")
# Load regime list
load("Regime_list.RData")


### 2. Check if the tree and data are good ---
# Get the data for two example traits
data_mvSLOUCH <- data_bivariate_mvOU_experiGR[,1:2]
data_mvSLOUCH <- as.matrix(data_mvSLOUCH)
# Check the tree
mvSLOUCH::phyltree_paths(tree_mvSLOUCH)$tree_height # = 1
# Check if rownames map to tip labels
isTRUE(all.equal(tree_mvSLOUCH$tip.label,rownames(data_mvSLOUCH)))
table(row.names(data_mvSLOUCH) == tree_mvSLOUCH$tip.label)

### 3. Define Atype  ---
Atype_custom <- list(
  Diagonal = matrix(
    c("+", "0",
      "0", "+"),
    nrow = 2,
    ncol = 2,
    byrow = T
  ),
  UpperTri = matrix(
    c("+", NA,
      "0", "+"),
    nrow = 2,
    ncol = 2,
    byrow = T
  ),
  LowerTri = matrix(
    c("+", "0",
      NA, "+"),
    nrow = 2,
    ncol = 2,
    byrow = T
  ),
  Any = matrix(
    c("+", NA,
      NA, "+"),
    nrow = 2,
    ncol = 2,
    byrow = T
  )
)

### 4. Define regime type ---
Rtype <- c("regime_trait1", "regime_trait2", "regime_trait1_trait2")

### 5. Test: Conduct an example model fitting for 1 time
OUOUmodel_start_regime <- ouchModel(tree_mvSLOUCH, data_mvSLOUCH, 
                                    regimes = regime_df$regime_trait1_trait2, 
                                    Atype="DecomposablePositive", # Assuming A type as DecomposablePositive
                                    Syytype="Diagonal", diagA=NULL, 
                                    maxiter=c(25,50))


#----------------------------------------------------------------------------------------#
#                1 Fit mvOU for all traits ; Syy = Diagonal ; without M.error            #
#----------------------------------------------------------------------------------------#

### 1. Define predictor names ---
predictors <- c("Metabolic_capacity", "Metabolic_versatility", 
                "Metabolic_completeness", "GC_percentage", "Genome_size", "rRNA_operon_copy_number", 
                "Niche_breadth", "Oxygen", "Temperature", "Salinity", "pH")

### 2. Initiate the result list ---
mvOU_BAYOU_SyyDiagonal_niter100 <- list()

### 3. Setup parallel work only once ---
ncores <- 100
cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  .libPaths("/home/haochh/Rlib")
  library(mvSLOUCH)
  library(PCMBaseCpp)
  library(tidyverse)
})

# Export objects that do not change across predictors
clusterExport(
  cl,
  varlist = c("tree_mvSLOUCH", "Rtype", "Atype_custom"),
  envir = environment()
)

### 4. Start mvOU modelling, make sure cluster stops even if error occurs ---
on.exit(stopCluster(cl), add = TRUE)

for (pre in predictors) {
  print(pre)
  var <- c("Growth_rate_experiment", pre)
  
  data_mvSLOUCH <- data_bivariate_mvOU_experiGR %>%
    select(all_of(var)) %>%
    as.matrix()
  
  regime_df <- Regime_list[[pre]]
  
  # Export only objects that change for each predictor
  clusterExport(
    cl,
    varlist = c("data_mvSLOUCH", "regime_df"),
    envir = environment()
  )
  
  niter <- 100
  
  mvOU_experiGR_Regime_BAYOU_SyyDiagonal_niter100 <- foreach(
    n = 1:niter,
    .packages = c("mvSLOUCH", "PCMBaseCpp"),
    .errorhandling = "pass"
  ) %dopar% {
    
    iter_res <- list()
    iter_res$iteration <- n
    
    ## ---------------------------
    ## 1) initial fits
    ## ---------------------------
    seed1_null <- 2025000 + n
    set.seed(seed1_null)
    
    start_null <- ouchModel(
      phyltree = tree_mvSLOUCH,
      mData = data_mvSLOUCH,
      regimes = NULL,
      Atype = "DecomposablePositive",
      Syytype = "Diagonal",
      diagA = NULL,
      maxiter = c(25, 50)
    )
    
    start_regime <- vector("list", length(Rtype))
    names(start_regime) <- Rtype
    
    for (i in seq_along(Rtype)) {
      regime_name <- Rtype[i]
      regime_vec <- regime_df[[regime_name]]
      
      seed1_regime <- 2025000 + n * 100 + i
      set.seed(seed1_regime)
      
      start_regime[[regime_name]] <- ouchModel(
        phyltree = tree_mvSLOUCH,
        mData = data_mvSLOUCH,
        regimes = regime_vec,
        Atype = "DecomposablePositive",
        Syytype = "Diagonal",
        diagA = NULL,
        maxiter = c(25, 50)
      )
    }
    
    iter_res$start <- list(
      null = start_null,
      regime = start_regime
    )
    
    ## ---------------------------
    ## 2) final fits using starts
    ## ---------------------------
    iter_res$final <- list()
    
    for (j in seq_along(Atype_custom)) {
      model_A <- Atype_custom[[j]]
      nameA <- names(Atype_custom)[j]
      
      iter_res$final[[nameA]] <- list()
      
      ## null model
      seed2_null <- 3035000 + n * 100 + j
      set.seed(seed2_null)
      
      fit_null <- ouchModel(
        phyltree = tree_mvSLOUCH,
        mData = data_mvSLOUCH,
        regimes = NULL,
        Atype = "Any",
        Syytype = "Diagonal",
        diagA = NULL,
        maxiter = c(100, 1000),
        start_point_for_optim = list(
          A = start_null$FinalFound$ParamsInModel$A,
          Syy = start_null$FinalFound$ParamsInModel$Syy
        ),
        parameter_signs = list(signsA = model_A)
      )
      
      iter_res$final[[nameA]][["null"]] <- list(
        model = fit_null,
        seed1 = seed1_null,
        seed2 = seed2_null
      )
      
      ## regime models
      for (i in seq_along(Rtype)) {
        regime_name <- Rtype[i]
        regime_vec <- regime_df[[regime_name]]
        
        seed2_regime <- 3035000 + n * 10000 + j * 100 + i
        set.seed(seed2_regime)
        
        fit_regime <- ouchModel(
          phyltree = tree_mvSLOUCH,
          mData = data_mvSLOUCH,
          regimes = regime_vec,
          Atype = "Any",
          Syytype = "Diagonal",
          diagA = NULL,
          maxiter = c(100, 1000),
          start_point_for_optim = list(
            A = start_regime[[regime_name]]$FinalFound$ParamsInModel$A,
            Syy = start_regime[[regime_name]]$FinalFound$ParamsInModel$Syy
          ),
          parameter_signs = list(signsA = model_A)
        )
        
        iter_res$final[[nameA]][[regime_name]] <- list(
          model = fit_regime,
          seed1 = 2025000 + n * 100 + i,
          seed2 = seed2_regime
        )
      }
    }
    
    iter_res
  }
  
  print(pre)
  
  mvOU_BAYOU_SyyDiagonal_niter100[[pre]] <-
    mvOU_experiGR_Regime_BAYOU_SyyDiagonal_niter100
}

### 5. Stop the parallel work ---
stopCluster(cl)

### 6. Get the results ---
mvOU_BAYOU_SyyDiagonal_niter100


#-------------------------------------------------------------------------------------#
#                2 Fit mvOU for all traits ; Syy = Diagonal ; with M.error            #
#-------------------------------------------------------------------------------------#

### 1. Define predictor names ---
predictors <- c("Metabolic_capacity", "Metabolic_versatility", 
                "Metabolic_completeness", "GC_percentage", "Genome_size", 
                "rRNA_operon_copy_number", "Niche_breadth", "Oxygen", 
                "Temperature", "Salinity", "pH")

### 2. Define measurement-error scenarios ---
### same error for both variables in each bivariate model
ME_list <- list(
  low  = 0.01,
  mid  = 0.05,
  high = 0.10
)

### 3. store results by predictor, then by M.error scenario ---
mvOU_BAYOU_SyyDiagonal_niter100_ME <- list()

for (pre in predictors) {
  print(pre)
  
  mvOU_BAYOU_SyyDiagonal_niter100_ME[[pre]] <- list()
  
  
  var <- c("Growth_rate_experiment", pre)
  
  
  data_mvSLOUCH <- data_bivariate_mvOU_experiGR %>%
    select(all_of(var))
  data_mvSLOUCH <- as.matrix(data_mvSLOUCH)
  

  regime_df <- Regime_list[[pre]]
  
  for (me_name in names(ME_list)) {
    
    M.error_current <- ME_list[[me_name]]
    
    # Setup for parallel work
    ncores <- 100
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    clusterEvalQ(cl, {
      .libPaths("/home/haochh/Rlib")
      library(mvSLOUCH)
      library(PCMBaseCpp)
      library(tidyverse)
    })
    
    # export current measurement error to workers
    clusterExport(
      cl,
      varlist = c("tree_mvSLOUCH", "data_mvSLOUCH", "regime_df", 
                  "Rtype", "Atype_custom", "M.error_current","me_name"),
      envir = environment()
    )
    
  
    niter <- 100
    
    # use a temporary object specific to this M.error scenario
    mvOU_experiGR_Regime_BAYOU_SyyDiagonal_niter100_ME <- foreach(
      n = 1:niter,
      .packages = c("mvSLOUCH", "PCMBaseCpp","tidyverse"),
      .errorhandling = "pass"
    ) %dopar% {
      
    
      library(dplyr)
      library(mvSLOUCH)
      library(PCMBaseCpp)
      
      iter_res <- list()
      iter_res$iteration <- n
      
      
      iter_res$M.error_name <- me_name
      iter_res$M.error_value <- M.error_current
      
      ## ---------------------------
      ## 1) initial fits
      ## ---------------------------
      seed1_null <- 2025000 + n
      set.seed(seed1_null)
      
      start_null <- ouchModel(
        phyltree = tree_mvSLOUCH,
        mData = data_mvSLOUCH,
        regimes = NULL,
        M.error = M.error_current,   ### ERROR
        Atype = "DecomposablePositive",
        Syytype = "Diagonal",
        diagA = NULL,
        maxiter = c(25, 50)
      )
      
      start_regime <- vector("list", length(Rtype))
      names(start_regime) <- Rtype
      
      for (i in seq_along(Rtype)) {
        regime_name <- Rtype[i]
        regime_vec <- regime_df[[regime_name]]
        
        seed1_regime <- 2025000 + n * 100 + i
        set.seed(seed1_regime)
        
        start_regime[[regime_name]] <- ouchModel(
          phyltree = tree_mvSLOUCH,
          mData = data_mvSLOUCH,
          regimes = regime_vec,
          M.error = M.error_current,   ### ERROR
          Atype = "DecomposablePositive",
          Syytype = "Diagonal",
          diagA = NULL,
          maxiter = c(25, 50)
        )
      }
      
      iter_res$start <- list(
        null = start_null,
        regime = start_regime
      )
      
      ## ---------------------------
      ## 2) final fits using starts
      ## ---------------------------
      iter_res$final <- list()
      
      for (j in 1:length(Atype_custom)) {
        model_A <- Atype_custom[[j]]
        nameA <- names(Atype_custom)[j]
        
        iter_res$final[[nameA]] <- list()
        
        ## null model: fit once per A
        seed2_null <- 3035000 + n * 100 + j
        set.seed(seed2_null)
        
        fit_null <- ouchModel(
          phyltree = tree_mvSLOUCH,
          mData = data_mvSLOUCH,
          regimes = NULL,
          M.error = M.error_current,   ### ERROR
          Atype = "Any",
          Syytype = "Diagonal",
          diagA = NULL,
          maxiter = c(100, 1000),
          start_point_for_optim = list(
            A = start_null$FinalFound$ParamsInModel$A,
            Syy = start_null$FinalFound$ParamsInModel$Syy
          ),
          parameter_signs=list(signsA=model_A)
        )
        
        iter_res$final[[nameA]][["null"]] <- list(
          model = fit_null,
          seed1 = seed1_null,
          seed2 = seed2_null
        )
        
        ## regime models
        for (i in seq_along(Rtype)) {
          regime_name <- Rtype[i]
          regime_vec <- regime_df[[regime_name]]
          
          seed2_regime <- 3035000 + n * 10000 + j * 100 + i
          set.seed(seed2_regime)
          
          fit_regime <- ouchModel(
            phyltree = tree_mvSLOUCH,
            mData = data_mvSLOUCH,
            regimes = regime_vec,
            M.error = M.error_current,   ### ERROR
            Atype = "Any",
            Syytype = "Diagonal",
            diagA = NULL,
            maxiter = c(100, 1000),
            start_point_for_optim = list(
              A = start_regime[[regime_name]]$FinalFound$ParamsInModel$A,
              Syy = start_regime[[regime_name]]$FinalFound$ParamsInModel$Syy
            ),
            parameter_signs=list(signsA=model_A)
          )
          
          iter_res$final[[nameA]][[regime_name]] <- list(
            model = fit_regime,
            seed1 = 2025000 + n * 100 + i,
            seed2 = seed2_regime
          )
        }
      }
      
      iter_res
    }
    
    stopCluster(cl)
    
    # save results by predictor and M.error scenario
    mvOU_BAYOU_SyyDiagonal_niter100_ME[[pre]][[me_name]] <- 
      mvOU_experiGR_Regime_BAYOU_SyyDiagonal_niter100_ME
    
    # progress message includes M.error scenario
    print(paste("Finished predictor:", pre, "| M.error:", me_name, "| value:", M.error_current))
  }
}

### 4. Get the results ---
mvOU_BAYOU_SyyDiagonal_niter100_ME


#----------------------------------------------------------------------------------------#
#                3 Fit mvOU for all traits ; Syy = UpperTri ; without M.error            #
#----------------------------------------------------------------------------------------#

### 1. Define predictor names ---
predictors <- c("Metabolic_capacity", "Metabolic_versatility", 
                "Metabolic_completeness", "GC_percentage", "Genome_size", "rRNA_operon_copy_number", 
                "Niche_breadth", "Oxygen", "Temperature", "Salinity", "pH")

### 2. Initiate the result list ---
mvOU_BAYOU_SyyUpperTri_niter100 <- list()

### 3. Setup parallel work only once ---
ncores <- 100
cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  .libPaths("/home/haochh/Rlib")
  library(mvSLOUCH)
  library(PCMBaseCpp)
  library(tidyverse)
})

# Export objects that do not change across predictors
clusterExport(
  cl,
  varlist = c("tree_mvSLOUCH", "Rtype", "Atype_custom"),
  envir = environment()
)

### 4. Start mvOU modelling, make sure cluster stops even if error occurs ---
on.exit(stopCluster(cl), add = TRUE)

for (pre in predictors) {
  print(pre)
  var <- c("Growth_rate_experiment", pre)
  
  data_mvSLOUCH <- data_bivariate_mvOU_experiGR %>%
    select(all_of(var)) %>%
    as.matrix()
  
  regime_df <- Regime_list[[pre]]
  
  # Export only objects that change for each predictor
  clusterExport(
    cl,
    varlist = c("data_mvSLOUCH", "regime_df"),
    envir = environment()
  )
  
  niter <- 100
  
  mvOU_experiGR_Regime_BAYOU_SyyUpperTri_niter100 <- foreach(
    n = 1:niter,
    .packages = c("mvSLOUCH", "PCMBaseCpp"),
    .errorhandling = "pass"
  ) %dopar% {
    
    iter_res <- list()
    iter_res$iteration <- n
    
    ## ---------------------------
    ## 1) initial fits
    ## ---------------------------
    seed1_null <- 2025000 + n
    set.seed(seed1_null)
    
    start_null <- ouchModel(
      phyltree = tree_mvSLOUCH,
      mData = data_mvSLOUCH,
      regimes = NULL,
      Atype = "DecomposablePositive",
      Syytype = "UpperTri",
      diagA = NULL,
      maxiter = c(25, 50)
    )
    
    start_regime <- vector("list", length(Rtype))
    names(start_regime) <- Rtype
    
    for (i in seq_along(Rtype)) {
      regime_name <- Rtype[i]
      regime_vec <- regime_df[[regime_name]]
      
      seed1_regime <- 2025000 + n * 100 + i
      set.seed(seed1_regime)
      
      start_regime[[regime_name]] <- ouchModel(
        phyltree = tree_mvSLOUCH,
        mData = data_mvSLOUCH,
        regimes = regime_vec,
        Atype = "DecomposablePositive",
        Syytype = "UpperTri",
        diagA = NULL,
        maxiter = c(25, 50)
      )
    }
    
    iter_res$start <- list(
      null = start_null,
      regime = start_regime
    )
    
    ## ---------------------------
    ## 2) final fits using starts
    ## ---------------------------
    iter_res$final <- list()
    
    for (j in seq_along(Atype_custom)) {
      model_A <- Atype_custom[[j]]
      nameA <- names(Atype_custom)[j]
      
      iter_res$final[[nameA]] <- list()
      
      ## null model
      seed2_null <- 3035000 + n * 100 + j
      set.seed(seed2_null)
      
      fit_null <- ouchModel(
        phyltree = tree_mvSLOUCH,
        mData = data_mvSLOUCH,
        regimes = NULL,
        Atype = "Any",
        Syytype = "UpperTri",
        diagA = NULL,
        maxiter = c(100, 1000),
        start_point_for_optim = list(
          A = start_null$FinalFound$ParamsInModel$A,
          Syy = start_null$FinalFound$ParamsInModel$Syy
        ),
        parameter_signs = list(signsA = model_A)
      )
      
      iter_res$final[[nameA]][["null"]] <- list(
        model = fit_null,
        seed1 = seed1_null,
        seed2 = seed2_null
      )
      
      ## regime models
      for (i in seq_along(Rtype)) {
        regime_name <- Rtype[i]
        regime_vec <- regime_df[[regime_name]]
        
        seed2_regime <- 3035000 + n * 10000 + j * 100 + i
        set.seed(seed2_regime)
        
        fit_regime <- ouchModel(
          phyltree = tree_mvSLOUCH,
          mData = data_mvSLOUCH,
          regimes = regime_vec,
          Atype = "Any",
          Syytype = "UpperTri",
          diagA = NULL,
          maxiter = c(100, 1000),
          start_point_for_optim = list(
            A = start_regime[[regime_name]]$FinalFound$ParamsInModel$A,
            Syy = start_regime[[regime_name]]$FinalFound$ParamsInModel$Syy
          ),
          parameter_signs = list(signsA = model_A)
        )
        
        iter_res$final[[nameA]][[regime_name]] <- list(
          model = fit_regime,
          seed1 = 2025000 + n * 100 + i,
          seed2 = seed2_regime
        )
      }
    }
    
    iter_res
  }
  
  print(pre)
  
  mvOU_BAYOU_SyyUpperTri_niter100[[pre]] <-
    mvOU_experiGR_Regime_BAYOU_SyyUpperTri_niter100
}

### 5. Stop the parallel work ---
stopCluster(cl)

### 6. Get the results ---
mvOU_BAYOU_SyyUpperTri_niter100



