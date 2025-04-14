FindOutliers <- function(
  object,
  vars,
  var.limit = list(),
  batch = NULL,
  type = c("auto", "both", "lower", "higher"),
  cut.1 = "median",
  cut.2 = "mad",
  n = 3,
  log = FALSE,
  ...
){
  for ( varx in vars ){
    if ( is.null(var.limit[[varx]]) ){
      var.limit[[varx]] <- c(min= NA, max = NA)
    }
  }
  param_value = Seurat::FetchData(object, vars = vars)
  batch <- as.vector(FetchData(object, vars = batch)[[1]])
  names(batch) <- rownames(param_value)
  outliers_list = list()
  for ( paramx in vars ){
    param_vector <- param_value[[paramx]]
    if ( min(param_vector) == max(param_vector) ) next
    names(param_vector) <- rownames(param_value)
    outliers_list[[paramx]] <- isOutlier(param_vector, batch = batch,
                                         limit = var.limit[[paramx]],
                                         n = n, type = type, cut.1 = cut.1,
                                         cut.2 = cut.2, log = log, ... )
  }
  return( outliers_list )
}

isOutlier <- function(
  metric,
  n = 3,
  type = c("both", "lower", "higher"),
  cut.1 = "median",
  cut.2 = "mad",
  log = FALSE,
  subset = NULL,
  batch = NULL,
  limit = NULL, 
  min.diff = NA
){
  if (log) {
    metric <- log2(metric)
    limit <- ifelse( !is.null(limit), log2(limit), limit)
  }
  type <- match.arg(type)
  N <- length(metric)
  if (nobatch <- is.null(batch)) {
    batch <- rep("1", N)
  } else {
    if (length(batch) != N) {
      stop("length of 'batch' must equal length of 'metric'")
    }
    batch <- as.character(batch)
  }
  if (!is.null(subset)) { t
    M <- metric[subset]
    B <- batch[subset]
  } else {
    M <- metric
    B <- batch
  }
  if (any(na.drop <- is.na(M))) {  
    M <- M[!na.drop]
    B <- B[!na.drop]
    warning("missing values ignored during outlier detection")
  }
  by.batch <- split(M, B)
  all_batches <- sort(unique(batch))
  empty <- rep(NA_real_, length(all_batches))
  names(empty) <- all_batches
  cur.1 <- empty
  cur.1[names(by.batch)] <- switch (cut.1,
    "median" = { unlist(lapply(by.batch, median)) },
    "mean" = { unlist(lapply(by.batch, mean)) }
  )
  cur.2 <- empty
  cur.2[names(by.batch)] <- switch (cut.2,
    "sd" = {unlist(mapply(sd, x=by.batch))  },
    "median" = { unlist(lapply(by.batch, median)) },
    "mad" = { unlist(mapply(mad, x=by.batch, center=cur.1[names(by.batch)], SIMPLIFY=FALSE)) }
  )
  diff.val <- pmax(min.diff, n * cur.2, na.rm = TRUE)
  batch.min <- sapply(by.batch, min) 
  lower.limit <- cur.1 - diff.val 
  for ( batchx in names(batch.min) ){
    lower.limit[batchx] <- ifelse( lower.limit[batchx] < batch.min[batchx],
                                   batch.min[batchx], lower.limit[batchx] )
    if ( !is.na(limit[1]) ) lower.limit[batchx] <- limit[1]
  }
  upper.limit <- cur.1 + diff.val
  batch.max <- sapply(by.batch, max)
  for ( batchx in names(batch.max) ){
    upper.limit[batchx] <- ifelse( upper.limit[batchx] > batch.max[batchx],
                                   batch.max[batchx], upper.limit[batchx] )
    if ( !is.na(limit[2]) ) upper.limit[batchx] <- limit[2]
  }
  if (type == "lower") {
    upper.limit[] <- Inf
  } else if (type == "higher") {
    lower.limit[] <- -Inf
  }
  collected <- (metric < lower.limit[batch] | upper.limit[batch] < metric)
  names(collected) <- names(metric)
  all.threshold <- rbind(lower=lower.limit, higher=upper.limit)
  output <- list(outliers=collected, thresholds=all.threshold)
  thresholds <- output$thresholds
  if (nobatch) {
    thresholds <- drop(thresholds)
  }
  if (log) {
    thresholds <- 2^thresholds
  }

  outliers <- output$outliers
  attr(outliers, "thresholds") <- thresholds
  outliers
}

RemoveDoublets <-function(
    object,
    doublet.rate,
    pN=0.25,
    pc.num=1:30
  ){
    sweep.res.list <- paramSweep_v3(object, PCs = pc.num, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
    bcmvn <- find.pK(sweep.stats) 
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    homotypic.prop <- modelHomotypic(object$seurat_clusters) 
    nExp_poi <- round(doublet.rate*ncol(object)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    seu.scored <- doubletFinder_v3(object, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                                   nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    cname <-colnames(seu.scored[[]])
    DF<-cname[grep('^DF',cname)]
    seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")
    seu.removed <- subset(seu.scored, subset = doublet != 1)
    p1 <- DimPlot(seu.scored, group.by = DF)
    res.list <- list("plot"=p1, "obj"=seu.removed)
    return(res.list) 
}