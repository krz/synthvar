# This function selects variables and calculates the synthetic control
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#library(mRMRe)

synthvar <- function(df = NULL,
                     time = NULL, # pre intervention time
                     treated = NULL,  # treated unit (integer)
                     control= NULL,  # control unit vector
                     method = "mrmre")  # default method = mrmre
  {
  # sanity checks
  if(is.null(df) == TRUE) {stop("\n No input data frame provided")}
  if(is.null(time) == TRUE) {stop("\n No time argument provided")}
  if(is.null(treated) == TRUE) {stop("\n No treated argument provided")}
  if(is.null(control) == TRUE) {stop("\n No control argument provided")}
  if(!("time" %in% colnames(df))) {stop("\n Please provide a time column")}
  if(!("unit" %in% colnames(df))) {stop("\n Please provide a unit column")}

  # geometry checks
  if(length(control) < 2) {stop("\n Please provide at least 2 control units")}
  if(length(treated) > 2) {stop("\n Please provide only one treated unit")}

  # split the df into a list of individual data frames, controls only, pre-intervention time points only
  controls <- list()
  for(i in 1:length(time)) {
    tmp <- df[df$time == time[i] & df$unit %in% control, ]
    tmp$time <- NULL
    tmp$unit <- NULL
    tmp <- tmp[,colSums(is.na(tmp))<nrow(tmp)]  # remove NA columns
    controls[[i]] <- tmp
  }
  # perform variables selection for each df in controls list
  # store the top variables in a list
  vars <- list()
  if(method == "mrmre") {
    for(i in 1:length(controls)){
      dd <- mRMR.data(data = controls[[i]])
      feats <- mRMR.ensemble(data = dd, target_indices = 1, solution_count = 5, feature_count = 5)
      vi <- data.frame('importance'=feats@mi_matrix[nrow(feats@mi_matrix),])
      vi$feature <- rownames(vi)
      row.names(vi) <- NULL
      vi <- na.omit(vi)
      vi <- vi[order(vi$importance, decreasing=TRUE),]
      vars[[i]] <- vi$feature[1:4]  # 4 top features
    }
  } else {
    stop(paste("\n", "unknown method:", method))
  }
  # final list of variables + outcome variable:
  impvars <- c(colnames(df)[1], unique(unlist(vars)))
  # make data frames with outcome and important variables
  # make treated vector
  X1 <- stack(df[(df$unit == treated) & (df$time %in% time), impvars])$values
  # make control matrix
  tmp <- df[(df$unit %in% control) & (df$time %in% time), c(impvars)]
  X0 <- sapply(1:length(tmp), function(i) unlist(tmp[i:(i+2), ]))

  #tmp <- df[(df$unit %in% control) & (df$time %in% time), c(impvars)]
  #l <- list()
  #for(i in 1:length(control)) {
  #  l[[i]] <- stack(tmp[i:(i+2), ])$values
  #}
  #X0 <- do.call(cbind, l)

  # scale X
  big.dataframe <- cbind(X0, X1)

  # big.dataframe <- big.dataframe[complete.cases(big.dataframe), ]
  # replace NA with 0
  big.dataframe[is.na(big.dataframe)] <- 0

  divisor <- sqrt(apply(big.dataframe, 1, var))
  scaled.matrix <-
    t(t(big.dataframe) %*% ( 1/(divisor) *
                               diag(rep(dim(big.dataframe)[1], 1)) ))

  X0.scaled <- scaled.matrix[, c(1:(dim(X0)[2]))]
  if(is.vector(X0.scaled) == TRUE) {X0.scaled <- t(as.matrix(X0.scaled))}
  X1.scaled <- scaled.matrix[, dim(scaled.matrix)[2]]

  # use quadprog to find synthetic controls
  Rinv <- solve(chol(t(X0.scaled) %*% X0.scaled))
  C <- cbind(rep(1, ncol(X0.scaled)), diag(ncol(X0.scaled)))
  b <- c(1,rep(0, ncol(X0.scaled)))
  d <- t(X1.scaled) %*% X0.scaled
  sol <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
  #round(sol$solution, 4)


  # print results
  cat("\n****************",
      "\n****************",
      "\n****************",
      "\n\nMSPE (LOSS V):", res$value,
      #      "\n\nLOSS (W):", loss.w,
      "\n\constrained solution:\n", round(as.numeric(sol$solution), 5),
      "\n\unconstrained solution:\n", round(as.numeric(sol$unconstrained.solution), 5),
      "\n\n"
  )

  optimize.out <- list(
    solution = round(as.numeric(sol$solution), 5),
    solution.unconstrained = round(as.numeric(sol$unconstrained.solution), 5),
    loss = res.value,
    lagrangian = res$Lagrangian
  )

  return(invisible(optimize.out))

}

