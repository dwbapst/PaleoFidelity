#' Live-Dead differences in alpha diversity and evenness
#'
#' FidelityDiv provides estimates of differences in alpha diversity and alpha evenness
#' between pairs of live and dead samples from single sites (or data pooled across
#' multiple sites) using community abundance data. For datasets including multiple site,
#' the function also returns means of live-dead differences across all sites or groups of sites.
#'
#' @details FidelityDiv assesses live-dead offsets in evenness/diversity between pairs of
#' sympatric live and dead samples using a bivariate approach described in Olszewski and
#' Kidwell (2007). The estimates of offsets in alpha diversity and evenness are returned
#' as two separate objects.
#'
#' (1) x - Live-dead offsets in alpha diversity for individual sites. The difference is
#' measured as the difference between natural logarithms of sample-standardized
#' species richness of sympatric dead and live samples:
#'
#' DELTA S = ln(S[dead]) - ln(S[live])
#'
#' A negative value indicates that the sample-standardized alpha diversity of a live sample
#' exceeds alpha diversity of the sympatric dead sample and a positive value indicates that a dead
#' sample is more diverse than live sample. Confidence intervals and p-values for the two-tailed null hypothesis:
#' Delta S = 0 are reported. A vegan function "rrarefy" is used to perform sample standardization
#'
#' (2) y - Live-dead offsets in evenness for individual sites. The difference is measured as the
#' difference between estimates of Hurlbert's PIE for live and dead samples:
#'
#' DELTA[PIE] = PIE[DEAD] - PIE[LIVE]
#'
#' A negative value indicates that evenness of live samples exceeds evenness of the sympatric
#' dead sample and a positive value indicates that a dead sample is more even than live sample.
#' Confidence intervals and p values for null hypothesis: Delta PIE = 0 are reported.
#'
#' @param live A matrix with counts of live-collected specimens (rows=sites, columns=taxa).
#'  Dimensions of 'live' and 'dead' matrices must match exactly.
#'
#' @param dead A matrix with counts of dead-collected specimens (rows=sites, columns=taxa).
#'  Dimensions of 'live' and 'dead' matrices must match exactly.
#'
#' @param gp An optional univariate factor defining groups of sites. The length of gp must
#'  equal number of rows of 'live' and 'dead' matrices.
#'
#' @param report Logical (default=FALSE) to print report generated by the function FidelitySummary
#'
#' @param n.filters An integer used to filter out small samples (default = 2, all samples with n > 1 kept).
#'  Note that for samples with n < 2, Hurlbert's PIE cannot be computed.
#'
#' @param t.filters An integer used to filter out rare taxa (default = 1, taxa >= 1 occurrence kept)
#' Note that removing rare taxa (t.filters > 1) is not advisable when measuring evenness or standardized richness.
#'
#' @param iter An integer defining number of replicate samples (default iter=100)
#'
#' @param CI A numerical value (default = 0.5) defining confidence bars for individual sites.
#'  Note: 0.5 - plots bars representing inter-quartile ranges, 0.95 - plots 95% confidence bars, etc.
#'  Confidence bars are estimated as percentiles of subsampled estimates of Delta S and Delta PIE.
#'
#' @param CImean A numerical value (default = 0.99) defining confidence bars for means of all sites
#'  or groups of sites (if 'gp' factor is provided). Note: 0.5 - bars representing inter-quartile
#'  range, 0.99 - plots 99% confidence bars, etc. Confidence bars are estimated as percentiles of
#'  subsampled estimates of Delta S and Delta PIE based on n (=iter) replicate subsamples.
#'
#' @param outdata Logical (default = FALSE) to determine if data files should be included in the output
#'
#' @param messages Logical (default = FALSE) to enable printing notes generated internally
#' by FidelitySummary function.
#'
#'
#' @return A list containing the following components:
#'   \item{live}{The post-processed version of 'live' data matrix used in all analyses}
#'   \item{dead}{The post-processed version of 'dead' data matrix used in all analyses}
#'   \item{gp}{The post-processed version of 'gp' factor, when provided}
#'   \item{x}{DELTA S values for each live-dead comparisons (site-level differences in sample standardized
#'            species richness)}
#'   \item{y}{DELTA PIE values for each live-dead comparisons (site-level differences in evenness
#'            estimated as Hurlbert's PIE)}
#'   \item{xmean}{Grand mean of Delta S}
#'   \item{ymean}{Grand mean of Delta PIE}
#'   \item{xgp}{Group means of Delta S (when 'gp' factor provided)}
#'   \item{ygp}{Group means of Delta PIE (when 'gp' factor provided)}
#'   \item{p.values}{p.values for Null H: Delta.S.p=0 and Delta.PIE.p=0}
#'   \item{p.gps}{per-group p.values for Null H: Delta.S.p=0 and Delta.PIE.p=0}
#'
#' @examples
#'
#' my.fid <- FidelityDiv(FidData$live[6:9,], FidData$dead[6:9,],
#' FidData$habitat[6:9], n.filters=20, iter=100, CI=0.95)
#' my.fid$x # site-level estimates of Delta S with 95% CIs and p values
#' my.fid$p.gps # p values for means of groups
#' AlphaPlot(my.fid, col.gp=c('forestgreen', 'coral1'), bgpt='beige')
#'
#' @export
#'
#' @importFrom vegan rrarefy
#'
#' @references Olszewski, T.D., and Kidwell, S.M., 2007, The preservational fidelity of evenness
#'             in molluscan death assemblages. Paleobiology 33: 1:23.

FidelityDiv <- function(live, dead, gp=NULL, report=FALSE, n.filters=2, t.filters=1,
                        iter=100, CI=0.5, CImean=0.99, outdata=FALSE, messages=F)
{

  # 1. Data assessment and filtering

  if (messages) out <- FidelitySummary(live, dead, gp, report=report, output=TRUE,
                         n.filters=n.filters, t.filters=t.filters) # check/filter data
  if (!messages) out <- suppressMessages(FidelitySummary(live, dead, gp, report=report,
                                                        output=TRUE, n.filters=n.filters,
                                                        t.filters=t.filters)) # check/filter data
  if (length(out) == 2) {live <- out$live;  dead <- out$dead}
  if (length(out) == 3) {live <- out$live;  dead <- out$dead; gp <- out$gp}
  if (min(apply(live, 1, sum)) < 2)
    stop("Hulbert's PIE cannot be computed for 'live' samples with n<2,
         filter out small samples using n.filters argument")
  if (min(apply(dead, 1, sum)) < 2)
    stop("Hulbert's PIE cannot be computed for 'dead' samples with n<2,
         filter out small samples using n.filters argument")

  # 2. Alpha Diversity/Evenness
  pie.f <- function(x) (sum(x) / (sum(x) - 1)) * (1 - sum((x / sum(x)) ^ 2))
  delta.alpha <- function(x, y, min) {
    a <- suppressWarnings(vegan::rrarefy(x, sample=min))
    b <- suppressWarnings(vegan::rrarefy(y, sample=min))
    cbind(
        log(apply(b, 1, function(z) sum(z > 0))) - log(apply(a, 1, function(z) sum(z > 0))),
        apply(b, 1, pie.f) - apply(a, 1, pie.f)
        )
  }
  min.sam <- apply(cbind(rowSums(live), rowSums(dead)), 1, min)
  out1 <- array(NA, dim=c(nrow(live), 2, iter))
  for (i in 1:iter) out1[,,i] <- delta.alpha(live, dead, min.sam)
   p1DS <- apply(rbind(out1[,1,]), 1, function(x) sum(x > 0))
   p2DS <- apply(rbind(out1[,1,]), 1, function(x) sum(x < 0))
   p1DP <- apply(rbind(out1[,2,]), 1, function(x) sum(x > 0))
   p2DP <- apply(rbind(out1[,2,]), 1, function(x) sum(x < 0))
   p.DS <- 2 * apply(cbind(p1DS, p2DS), 1, min) / iter
   p.DP <- 2 * apply(cbind(p1DP, p2DP), 1, min) / iter
  if (sum(p.DS == 0) > 0) p.DS[p.DS == 0] <- 1 / iter
  if (sum(p.DP == 0) > 0) p.DP[p.DP == 0] <- 1 / iter
   DS <- cbind(n.std=min.sam, est=rowMeans(rbind(out1[,1,])),
              t(apply(rbind(out1[,1,]), 1, stats::quantile, prob=c((1 - CI) / 2, 1 - (1 - CI) / 2))),
              p=p.DS)
   DP <- cbind(n.std=min.sam, est=rowMeans(rbind(out1[,2,])),
              t(apply(rbind(out1[,2,]), 1, stats::quantile, prob=c((1 - CI) / 2, 1 - (1 - CI) / 2))),
              p=p.DP)

# 3.Means for all data and by groups (if 'gp' factor is provided)
  # All data
   outDS3 <- NULL
   outDP3 <- NULL
   p.GP <- NULL
if (nrow(live) > 1) {
  allS <- colMeans(rbind(out1[,1,]))
  allP <- colMeans(rbind(out1[,2,]))
  meanDS <- c(mean(allS), stats::quantile(allS, prob=c((1-CImean)/2, 1 - (1-CImean)/2)))
  meanDP <- c(mean(allP), stats::quantile(allP, prob=c((1-CImean)/2, 1 - (1-CImean)/2)))
  allpS <- 2 * min(sum(allS > 0), sum(allS < 0))
  ifelse(allpS == 0, Delta.S.p <- 1 / iter, Delta.S.p <- allpS / iter)
  allpP <- 2 * min(sum(allP > 0), sum(allP < 0))
  ifelse(allpP == 0, Delta.PIE.p <- 1 / iter, Delta.PIE.p <- allpP / iter)
  # By groups
  if (length(gp) > 0) {
   outDS <- sapply(as.data.frame(out1[,1,]), function(x) tapply(x, gp, mean))
   outDS2 <- cbind(n.sites=tapply(gp, gp, length), est=apply(outDS, 1, mean),
            t(apply(outDS, 1, stats::quantile, prob=c((1-CImean)/2, 1 - (1-CImean)/2))))
   outDP <- sapply(as.data.frame(out1[,2,]), function(x) tapply(x, gp, mean))
   outDP2 <- cbind(n.sites=tapply(gp, gp, length), est=apply(outDP, 1, mean),
                  t(apply(outDP, 1, stats::quantile, prob=c((1-CImean)/2, 1 - (1-CImean)/2))))
   gpDS1 <- apply(outDS, 1, function(x) sum(x<0))
   gpDS2 <- apply(outDS, 1, function(x) sum(x>0))
   gpDP1 <- apply(outDP, 1, function(x) sum(x<0))
   gpDP2 <- apply(outDP, 1, function(x) sum(x>0))
   p.gp.DS <- 2 * apply(cbind(gpDS1, gpDS2), 1, min) / iter
   p.gp.DP <- 2 * apply(cbind(gpDP1, gpDP2), 1, min) / iter
   if (sum(p.gp.DS == 0) > 0) p.gp.DS[p.gp.DS == 0] <- 1 / iter
   if (sum(p.gp.DP == 0) > 0) p.gp.DP[p.gp.DP == 0] <- 1 / iter
   if (length(gp) > 0) outDS3 <- outDS2
   if (length(gp) > 0) outDP3 <- outDP2
   if (length(gp) > 0)  p.GP <- cbind(p.Delta.S=p.gp.DS, p.Delta.PIE=p.gp.DP)
   outDS3 <- data.frame(outDS3, group=levels(gp), measure='Delta S')
   colnames(outDS3)[3:4] <- c(paste(1 - CImean, '%', sep=''), paste(CImean, '%', sep=''))
   outDP3 <- data.frame(outDP3, group=levels(gp), measure='Delta PIE')
   colnames(outDP3)[3:4] <- c(paste(1 - CImean, '%', sep=''), paste(CImean, '%', sep=''))
   rownames(outDP3) <- NULL
   rownames(outDS3) <- NULL
 }
}

if (nrow(live) == 1)  meanDS <- meanDP <- outDS3 <- outDP3 <- Delta.S.p <- Delta.PIE.p <- p.GP <- NA

if (outdata) out1 <- list(live=live, dead=dead, gp=gp, out=out1, x=DS, y=DP, xmean=meanDS, ymean=meanDP,
                          xgp=outDS3, ygp=outDP3, p.values=cbind(Delta.S.p, Delta.PIE.p), p.gps=p.GP)
if (!outdata) out1 <- list(gp=gp, x=DS, y=DP, xmean=meanDS, ymean=meanDP, xgp=outDS3, ygp=outDP3,
                           p.values=cbind(Delta.S.p, Delta.PIE.p), p.gps=p.GP)
  class(out1) <- append(class(out1),"FidelityDiv")
  return(out1)
}
