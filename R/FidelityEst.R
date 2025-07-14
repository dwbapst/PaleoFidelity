#' Compositional measures of live-dead fidelity
#'
#' FidelityEst estimates compositional fidelity by comparing two matching
#' (live and dead) matrices with community abundance data. The function returns
#' fidelity measures for individual sites, mean measures across sites, and means
#' for groups of sites. The function also returns sample-standardized, bias-corrected
#' and 'perfect fidelity' estimates.
#'
#' @details FidelityEst assesses compositional fidelity using
#' measures of correlation/associations/similarity.
#'
#' (1) x - a measure of correlation/association: Spearman, Kendall, or Pearson
#' (2) y - an abundance-based index of similarity such as Bray or Jaccard-Chao
#'
#' Because fidelity measures are sensitive to under-sampling and unbalanced sampling,
#' FidelityEst function attempts to correct for sampling biases by (1) estimating data-specific
#' biases or (2) standardizing sampling coverage. In the first approach, the bias is estimated
#' using a resampling protocol under the perfect fidelity (PF) model, in which
#' pooled data (live + dead) are randomly partitioned into replicate pairs of samples
#' (using sample sizes of original samples), thus creating sample pairs derived from
#' single underlying rank abundance species distributions (i.e., the perfect fidelity).
#' For an unbiased estimator, the resampled fidelity measures should indicate perfect
#' fidelity (e.g., Spearman rho = 1). The offset between the expected observed
#' PF value (1 - PF) provides a data-specific estimate of sampling bias. The adjusted
#' fidelity measure is then given by Adjusted = Observed + (1 - PF). Replicate samples
#' produce a distribution of PF values and resulting adjusted fidelity measures,
#' from which confidence intervals and significance tests can be derived.
#' In the second approach, fidelity measures are computed for sample standardized data,
#' where all samples are subsampled to a sample size given by the smallest sample.
#' Replicate resampling produces a distribution of sample-standardized fidelity estimates
#' used to generate confidence intervals and means for standardized fidelity estimates.
#'
#' @param live A matrix with counts of live-collected specimens (rows=sites, columns=species
#' or other variable). Dimensions and 'rownames' and 'colnames' of 'live' and 'dead' matrices
#' must match exactly.
#'
#' @param dead A matrix with counts of dead-collected specimens (rows=sites, columns=species
#' or other variables). Dimensions and 'rownames' and 'colnames' of 'live' and 'dead' matrices
#' must match exactly.
#'
#' @param gp An optional univariate factor defining groups of sites. The length of 'gp' must
#' equal number of rows of 'live' and 'dead' matrices.
#'
#' @param cor.measure A character string (default='spearman') defining correlation measure
#' (passed on to \code{\link[stats]{cor}} function) used to estimate live-dead correlations.
#'
#' @param sim.measure A character string (default='chao') defining similarity measure (passed
#' on to \code{\link[vegan]{vegdist}}) used to estimate live-dead similarity.
#'
#' @param n.filters An integer used to filter out small samples (default n.filters=0,
#' all samples kept).
#'
#' @param t.filters An integer used to filter out rare taxa (default t.filters=1,
#' all taxa with at least one occurrence kept).
#'
#' @param report Logical (default report = FALSE) (suppresses notes, warnings and data summary).
#'
#' @param iter An integer defining number of replicate samples
#' for perfect fidelity resampling model (default iter = 10).
#'
#' @param iter2 An integer defining number of resampling iteration for subsampling
#' standardization (default iter2 = 10).
#'
#' @param min.sam An integer defining number of specimens for
#' sample standardization (default = 30).
#'
#' @param CI A numerical value (default = 0.95) defining confidence limits for
#' adjusted and sample-standardized estimates of fidelity based on percentiles of
#' resampled estimates of sample-standardized, bias-corrected, or 'perfect fidelity'
#' measures of correlation/similarity.
#'
#' @param rm.zero Logical (default rm.zero = FALSE) removes double 0's when
#' computing correlation measure.
#'
#' @param tfsd A character string (default='none') specifying data standardization
#' or transformations (applicable only for some similarity measure).The following options are
#' available: 'none' (or any unused character string) - no standardization/transformation,
#' 'total' - relative abundance, 'wisconsin' - double relativization,
#' 'r4' - 4th root transformation, 'log' - ecological log-transformation,
#' 'total4' - 4th root transformation of relative abundances.
#'
#' @return A list containing the following components:
#'   \item{x}{Live-dead correlation coefficients for each site}
#'   \item{y}{Live-dead similarity coefficients for each live-dead comparison}
#'   \item{xc}{Adjusted live-dead correlation coefficients, confidence intervals,
#'    and estimated live-dead correlation coefficients for perfect fidelity model}
#'   \item{yc}{Adjusted live-dead similarity coefficients, confidence intervals,
#'    and estimated live-dead similarity coefficients for perfect fidelity model}
#'   \item{xs}{Sample-standardized live-dead correlation coefficients,
#'   confidence intervals, and standardized sample size}
#'   \item{ys}{Sample-standardized live-dead similarity coefficients,
#'   confidence intervals, and standardized sample size}
#'   \item{x.stats}{Statistical summary for raw correlation coefficients for
#'   all data and for each group when 'gp' factor provided}
#'   \item{y.stats}{Statistical summary for raw similarity coefficients for
#'   all data and for each group when 'gp' factor provided}
#'   \item{xc.stats}{Statistical summary for adjusted correlation coefficients for
#'   all data and for each group when 'gp' factor provided}
#'   \item{yc.stats}{Statistical summary for adjusted similarity coefficients for
#'   all data and for each group when 'gp' factor provided}
#'   \item{xs.stats}{Statistical summary for sample-standardized correlation
#'   coefficients for all data and for each group when 'gp' factor provided}
#'   \item{ys.stats}{Statistical summary for sample-standardized similarity
#'   coefficients for all data and for each group when 'gp' factor provided}
#'   \item{x.pf.dist}{Distributions of randomized correlation values under perfect
#'     fidelity model for each of the live-dead sample comparison}
#'   \item{y.pf.dist}{Distributions of randomized correlation values under perfect
#'     fidelity model for each of the live-dead sample comparison}
#'   \item{xc.dist}{Distributions of model adjusted correlation values for
#'    each of the samples}
#'    \item{yc.dist}{Distributions of model adjusted similarity values for
#'    each of the samples}
#'   \item{live}{The post-processed version of 'live' data matrix used in analyses}
#'   \item{dead}{The post-processed version of 'dead' data matrix used in analyses}
#'   \item{gp}{The post-processed version of 'gp', when 'gp' factor provided}
#'   \item{values}{A list with values of parameters used in the analysis}
#'
#' @examples
#'
#' data(FidData)
#' out1 <- FidelityEst(live = FidData$live[6:9,], dead = FidData$dead[6:9,],
#'                     gp = FidData$habitat[6:9], cor.measure='spearman',
#'                     sim.measure='bray', n.filters=30, iter=99, rm.zero=FALSE, tfsd='none')
#' SJPlot(out1, gpcol=c('forestgreen', 'coral3'))
#'
#' @export
#'
#' @importFrom vegan vegdist
#'

FidelityEst <- function (
    live, dead, 
    gp = NULL, 
    cor.measure = "spearman",
    sim.measure = "chao", 
    n.filters = 0, 
    t.filters = 1,
    report = FALSE, 
    iter = 10, 
    iter2 = 10, 
    min.sam = 30, 
    CI = 0.95,
    rm.zero = TRUE, 
    tfsd = "total",
    perf.rescale = "scale",
    ignore.na = TRUE
    ){
    
    out <- FidelitySummary(
        live, dead, 
        gp, 
        report = report, 
        output = TRUE,
        n.filters = n.filters, 
        t.filters = t.filters)
        
    if (length(out) == 2) {
        live <- out$live
        dead <- out$dead
        }
    
    if (length(out) == 3) {
        live <- out$live
        dead <- out$dead
        gp <- out$gp
        }
  
    my.cor.F <- function(x, y) {
        if (rm.zero) {
            if (sum(colSums(rbind(x, y)) == 0) > 0) {
                good.taxa <- which(colSums(rbind(x, y)) > 0)
                x <- x[good.taxa, drop = FALSE]
                y <- y[good.taxa, drop = FALSE]
                }
            }
        stats::cor(x, y, method = cor.measure)
        }

    cor.e <- mapply(
        function(x, y){
            my.cor.F(as.vector(x), as.vector(y))
            }, 
        as.data.frame(t(live)), 
        as.data.frame(t(dead))
        )
    
    tfsd_standardize <- function(input, tfsd){
        if(all(tfsd != 
               c("total", "total4", "wisconsin", "log", "r4", "none"))){
            stop("tfsd is not an allowed option")
            }
        #
        if (tfsd == "total") {
            out <- vegan::decostand(input, "total")
            }
        if (tfsd == "total4") {
            out <- vegan::decostand(input, "total")^(1/4)
            }
        if (tfsd == "wisconsin") {
            out <- vegan::wisconsin(input)
            }
        if (tfsd == "log") {
            out <- vegan::decostand(input, "log")
            }
        if (tfsd == "r4") {
            out <- input^0.25
            }
        if(tfsd == "none"){
            out <- input
            }
        if(is.null(out)){
            stop("bad transform")
            }
        #if(any(rowSums(out) < 1)){
        #    stop("why are there empty rows in out")}
        return(out)
        }
        
    x3 <- tfsd_standardize(live, tfsd)
    x4 <- tfsd_standardize(dead, tfsd)    

    sim.e <- mapply(
        function(x, y){ 
            1 - vegan::vegdist(rbind(x, y), method = sim.measure)
            }, 
        as.data.frame(t(x3)), 
        as.data.frame(t(x4))
        )
    
    fid.sam <- cbind(cor.e, sim.e)
    colnames(fid.sam) <- c(cor.measure, sim.measure)
    rownames(fid.sam) <- rownames(live)
    raw.stats <- cbind(mean(cor.e), mean(sim.e))
    colnames(raw.stats) <- c(cor.measure, sim.measure)
    rownames(raw.stats) <- "all samples"
  
    if (length(gp) == 0){
        mean.gp <- NA
        }
    
    if (length(gp) > 0 & nrow(live) > 1) {
        mean.gp <- cbind(
            tapply(cor.e, gp, mean), 
            tapply(sim.e, gp, mean)
            )
        colnames(mean.gp) <- c(cor.measure, sim.measure)
        raw.stats <- rbind(raw.stats, mean.gp)
        }
  
    raw.x <- cbind(raw.stats[, 1])
    raw.y <- cbind(raw.stats[, 2])
  
    if (iter2 >= 10) {
    
        fid.subsam <- function(x, y, min) {
            a <- suppressWarnings(vegan::rrarefy(x, sample = min))
            b <- suppressWarnings(vegan::rrarefy(y, sample = min))
            cor.ab <- mapply(
                function(x, y){
                    my.cor.F(as.vector(x), as.vector(y))
                    }, 
                as.data.frame(t(a)), 
                as.data.frame(t(b))
            )
            #
            #if(any(is.na(cor.ab))){
            #    # NAs mainly occur when too few species are sampled, or
            #        # all species are sampled with identical abundance
            #    warning(
            #        "Some rarefied correlation coefs returned as NA, converted to zero instead"
            #        )
            #    cor.ab[is.na(cor.ab)] <- -999
            #    #stop(paste0("can't calculate correlation for samples: ", 
            #    #            paste0(which(is.na(cor.ab),",")))
            #    }
            #
            
            #if(any(rowSums(a))<1){
            #    stop("why are there empty rows in a")}
            #if(any(rowSums(b))<1){
            #    stop("why are there empty rows in b")}               
            
            c <- tfsd_standardize(x, tfsd)
            d <- tfsd_standardize(y, tfsd)
            
            sim.ab <- mapply(
              function(x, y){
                1 - vegan::vegdist(
                    rbind(x,y), 
                    method = sim.measure)
              }, 
              as.data.frame(t(c)), 
              as.data.frame(t(d))
            )
            
            #if(length(cor.ab)<ncol(a)){stop("woopsie")}
            #if(length(sim.ab)<ncol(b)){stop("oopsie")}
            
            out <- cbind(cor.ab, sim.ab)
            return(out)
            }
  
    outSS <- array(NA, dim = c(nrow(live), 2, iter2))
    
    for (i in 1:iter2) {
        outSS[, , i] <- fid.subsam(
            live, dead,
            min = min.sam
            )
        }
        
    #if CI is not NA
    #if(!is.na(CI)){

    LCL <- (1 - CI)/2
    UCL <- 1 - LCL

    SSCOR <- cbind(
        est = rowMeans(rbind(outSS[, 1, ])), 
        t(apply(
            rbind(outSS[,1, ]), 1, 
            stats::quantile, 
            prob = c((1 - CI)/2, 1-(1-CI)/2),
            na.rm = ignore.na
            )), 
        n.std = min.sam)
    
    SSSIM <- cbind(
        est = rowMeans(rbind(outSS[, 2, ])), 
        t(apply(
            rbind(outSS[,2, ]), 1, 
            stats::quantile, 
            prob = c((1 - CI)/2, 1 - (1 - CI)/2),
            na.rm = ignore.na
            )), 
        n.std = min.sam)
    
    colnames(SSCOR)[1:3] <- c(
        paste(cor.measure, ".SS", sep = ""), 
        paste(LCL), paste(UCL)
        )
    
    colnames(SSSIM)[1:3] <- c(
        paste(sim.measure, ".SS", sep = ""), 
        paste(LCL), paste(UCL)
        )

    rownames(SSCOR) <- rownames(live)
    rownames(SSSIM) <- rownames(live)

    ifelse(
        nrow(live) > 1, 
        xss <- colMeans(outSS[, 1, ]),
        xss <- mean(outSS[, 1, ])
        )

    ifelse(
        nrow(live) > 1, 
        yss <- colMeans(outSS[, 2, ]),
        yss <- mean(outSS[, 2, ])
        )

    xss2 <- stats::quantile(xss, 
        prob = c(LCL, 0.5, UCL),
        na.rm = ignore.na)
    yss2 <- stats::quantile(yss, 
        prob = c(LCL, 0.5, UCL),
        na.rm = ignore.na)

    xs.stat <- cbind(mean = mean(xss), LCL = xss2[1],
                          median = xss2[2], UCL = xss2[3])
    rownames(xs.stat) <- "all samples"

    ys.stat <- cbind(mean = mean(yss), LCL = yss2[1],
                     median = yss2[2], UCL = yss2[3])
    rownames(ys.stat) <- "all samples"
        
    if (length(gp) > 0 & nrow(live) > 1) {
        xsg <- apply(
          outSS[, 1, ], 2, 
          function(x) tapply(x, gp, mean))
        xsg2 <- t(rbind(
          apply(xsg, 1, mean), 
          apply(xsg, 1, 
                stats::quantile, 
                prob = c(LCL, 0.5, UCL),
                na.rm = ignore.na)
          ))
        colnames(xsg2) <- c("mean", "LCL", "median", "UCL")
        
        xs.stat <- rbind(`all samples` = xs.stat, xsg2)
        
        ysg <- apply(outSS[, 2, ], 2, 
                   function(x) tapply(x, gp, mean)
                   )
        
        ysg2 <- t(rbind(
          apply(ysg, 1, mean), 
          apply(ysg, 1, stats::quantile, 
                prob = c(LCL, 0.5, UCL),
                na.rm = ignore.na)
          ))
        
        colnames(ysg2) <- c("mean", "LCL", "median", "UCL")
        
        ys.stat <- rbind(`all samples` = ys.stat, ysg2)
        }

    colnames(xs.stat)[c(2, 4)] <- c(paste(LCL), paste(UCL))
    colnames(ys.stat)[c(2, 4)] <- c(paste(LCL), paste(UCL))

    
    ss.cor.rep <- t(outSS[, 1, ])
    ss.sim.rep <- t(outSS[, 2, ])

    colnames(ss.cor.rep) <- rownames(live)
    colnames(ss.sim.rep) <- rownames(live)
    }

    if (iter >= 10) {
        FidPerfModel <- function(live, dead) {
            pooled <- live + dead
            rlive <- as.vector(
                vegan::rrarefy(x = pooled, sum(live))
                )
            rdead <- as.vector(
                vegan::rrarefy(x = pooled, sum(dead))
                )
        r.cor.e = my.cor.F(rlive, rdead)
        rx3 <- tfsd_standardize(rlive, tfsd)
        rx4 <- tfsd_standardize(rdead, tfsd)

        r.sim.e = 1 - vegan::vegdist(t(cbind(rx3, rx4)),
            method = sim.measure)
        return(c(r.cor.e, r.sim.e))
        }
    pf.output <- array(0, dim = c(iter, nrow(live), 2), 
        dimnames = list(1:iter, rownames(live), c(cor.measure, sim.measure)))
    for (i in 1:iter) {
        for (j in 1:nrow(live)) {
            pf.output[i, j, ] <- FidPerfModel(rbind(live[j,
            ]), rbind(dead[j, ]))
            }
        }
    
    perfidest <- apply(pf.output, c(2, 3), stats::median)
    
    cor.obs.rep <- matrix(fid.sam[, 1], iter, 
        length(fid.sam[, 1]), byrow = T)
    sim.obs.rep <- matrix(fid.sam[, 2], iter, 
        length(fid.sam[, 2]), byrow = T)
    
    # perfect fidelity rescaling
    if(perf.rescale == "shift"){
        cor.m.adj <- (1 - pf.output[, , 1]) + cor.obs.rep
        sim.m.adj <- (1 - pf.output[, , 2]) + sim.obs.rep
    
        cor.m.adj[cor.m.adj > 1] <- 1
        sim.m.adj[sim.m.adj > 1] <- 1
        }
    
    if(perf.rescale == "scale"){
        cor.m.adj <- (cor.obs.rep + 1)/2
        cor.m.adj <- cor.m.adj/((pf.output[, , 1]+1)/2)
        cor.m.adj <- (cor.m.adj*2)-1
        
        sim.m.adj <- sim.obs.rep/pf.output[, , 2]
        }
    
    cor.adj.sam <- apply(cor.m.adj, 2, mean)
    sim.adj.sam <- apply(sim.m.adj, 2, mean)
    
    cor.adj.sam.CI <- apply(cor.m.adj, 2, 
        stats::quantile, 
        prob = c(LCL, UCL),
        na.rm = ignore.na
        )
    sim.adj.sam.CI <- apply(sim.m.adj, 2, 
        stats::quantile, 
        prob = c(LCL, UCL),
        na.rm = ignore.na
        )
    
    corrected <- cbind(cor.adj.sam, sim.adj.sam)
    corrected.mean <- c(mean(cor.m.adj), mean(sim.m.adj))
    
    names(corrected.mean) <- c(cor.measure, sim.measure)
    xc.stat <- c(mean(corrected[, 1]), 
        stats::quantile(corrected[, 1], 
            prob = c(LCL, 0.5, UCL),
            na.rm = ignore.na
            ))
    yc.stat <- c(mean(corrected[, 2]), 
        stats::quantile(corrected[, 2], 
            prob = c(LCL, 0.5, UCL),
            na.rm = ignore.na
            ))
    
    if (length(gp) == 0) {
        xc.stat <- rbind(xc.stat)
        yc.stat <- rbind(yc.stat)
        colnames(xc.stat) <- c("mean", paste(LCL),
                             "median", paste(UCL))
        colnames(yc.stat) <- c("mean", paste(LCL),
                             "median", paste(UCL))
        rownames(xc.stat) <- "all samples"
        rownames(yc.stat) <- "all samples"
        }
    if (length(gp) > 0 & nrow(live) > 1) {
      statsgp.x1 <- tapply(corrected[, 1], gp, mean)
      statsgp.x2 <- tapply(corrected[, 1], gp, stats::quantile,
        prob = c(LCL, 0.5, UCL),
        na.rm = ignore.na)
      statsgp.x3 <- cbind(statsgp.x1, matrix(unlist(statsgp.x2),
        length(levels(gp)), 3, byrow = T))
      statsgp.y1 <- tapply(corrected[, 2], gp, mean)
      statsgp.y2 <- tapply(corrected[, 2], gp, stats::quantile,
        prob = c(LCL, 0.5, UCL),
        na.rm = ignore.na)
      statsgp.y3 <- cbind(statsgp.y1, 
        matrix(unlist(statsgp.y2),
        length(levels(gp)), 3, byrow = T))
      xc.stat <- rbind(all.samples = xc.stat, statsgp.x3)
      yc.stat <- rbind(all.samples = yc.stat, statsgp.y3)
      colnames(xc.stat) <- c("mean", paste(LCL),
        "median", paste(UCL))
      colnames(yc.stat) <- c("mean", paste(LCL),
        "median", paste(UCL))
        }
    xc.stat <- as.matrix(xc.stat)
    yc.stat <- as.matrix(yc.stat)
    xc.sum <- cbind(corrected[, 1], t(cor.adj.sam.CI),
                         perfidest[, 1])
    colnames(xc.sum) <- c(paste(cor.measure, ".ADJ", sep = ""),
                       paste(LCL), paste(UCL), 
                       paste(cor.measure, ".PF", sep = ""))
    yc.sum <- cbind(corrected[, 2], t(sim.adj.sam.CI),
                         perfidest[, 2])
    colnames(yc.sum) <- c(paste(sim.measure, ".ADJ", sep = ""),
                       paste(LCL), paste(UCL), 
                       paste(sim.measure, ".PF", sep = ""))
    
    if (nrow(live) == 1) {
        raw.stats <- c(NA, NA)
        xc.stat <- NA
        yc.stat <- NA
        xs.stat <- NA
        ys.stat <- NA
        }
    }

    if (iter < 10) {
        xc.sum <- yc.sum <- xc.stat <- yc.stat <- cor.m.adj <- sim.m.adj <- "NA"
        pf.output <- array(dim = c(1, 1, 2))
        }

    if (iter2 < 10) {
        SSCOR <- SSSIM <- xs.stat <- ys.stat <- ss.cor.rep <- ss.sim.rep <- "NA"
        }
    
  out1 <- list(x = fid.sam[,1], y = fid.sam[,2], 
               xc = xc.sum, yc = yc.sum, 
               xs = SSCOR, ys = SSSIM, 
               x.stats = raw.x, y.stats = raw.y, 
               xc.stats = xc.stat, yc.stats = yc.stat,
               xs.stats = xs.stat, ys.stats = ys.stat, 
               x.pf.dist = pf.output[, , 1], y.pf.dist = pf.output[, , 2], 
               xc.dist = cor.m.adj, yc.dist = sim.m.adj, 
               xs.dist = ss.cor.rep, ys.dist = ss.sim.rep,
               live = live, dead = dead, 
               gp = gp, 
               values = list(measures = c(cor.measure, sim.measure),
                    data.transf = tfsd, remove.double.zeros = rm.zero,
                    PFiter = iter, SSIter = iter2, min.sam = min.sam)
               )

    class(out1) <- append(class(out1), "FidelityEst")
    return(out1)
    }
