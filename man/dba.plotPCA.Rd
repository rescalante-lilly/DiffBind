\name{dba.plotPCA}
\alias{dba.plotPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PCA plot
}
\description{
Principal Component Analysis plot
}
\usage{
dba.plotPCA(DBA, attributes, minval, maxval,
           contrast, method=DBA$config$AnalysisMethod, 
           th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
           report, score, mask, sites, cor=FALSE,
           b3D=FALSE, vColors, dotSize, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}



\item{DBA}{
DBA object.
}

%\subsection{Optional/defaulted arguments}{}

\item{attributes}{
attribute or vector of attributes to use to color plotted points. Each unique combination of attribute values will be assigned a color. Chosen from:
\itemize{ 
\item \link{DBA_GROUP}
\item \link{DBA_ID}
\item \link{DBA_TISSUE}
\item \link{DBA_FACTOR}
\item \link{DBA_CONDITION}
\item \link{DBA_TREATMENT}
\item \link{DBA_REPLICATE}
\item \link{DBA_CONSENSUS}
\item \link{DBA_CALLER}
}
Note that DBA_GROUP is a special attribute which will result in samples from each group in a contrast being colored separately.
}

\item{minval}{
Set all scores less than this to minval
}
\item{maxval}{
Set all scores greater than this to maxval
}

\item{contrast}{
number of contrast to use for PCA; if present, plots a PCA based on a differential binding affinity analysis (see  \code{\link{dba.analyze}}). If mask is unspecified, only the samples in the contrast will be included.  See  \code{\link{dba.show}(DBA, bContrast=T)} to get contrast numbers. If missing, uses scores in the main binding matrix.
}
\item{method}{
method used for analysis (used in conjunction with contrast):
\itemize{ 
\item \link{DBA_EDGER}
\item \link{DBA_DESEQ}
\item \link{DBA_DESEQ2}
\item \link{DBA_EDGER_BLOCK}
\item \link{DBA_DESEQ_BLOCK}
\item \link{DBA_DESEQ2_BLOCK}
}
}

\item{th}{
significance threshold; all sites with FDR (or p-values, see bUsePval) less than or equal to this value will be included in the PCA, subject to maxVal. Used in conjunction with contrast.
}
\item{bUsePval}{
if TRUE, uses p-value instead of FDR for thresholding. Used in conjunction with contrast.
}
\item{report}{
report (obtained from  \code{\link{dba.report}}) specifying the data to be used . If this is present, the method, th, and bUsePval parameters are ignored.
}

\item{score}{
Score to use for count data. Only used when plotting the global binding matrix (no contrast specified). One of:
\itemize{
\item \link{DBA_SCORE_READS}
\item \link{DBA_SCORE_READS_MINUS}
\item \link{DBA_SCORE_READS_FOLD}
\item \link{DBA_SCORE_RPKM}
\item \link{DBA_SCORE_RPKM_FOLD}
\item \link{DBA_SCORE_TMM_READS_FULL}
\item \link{DBA_SCORE_TMM_READS_EFFECTIVE}
\item \link{DBA_SCORE_TMM_MINUS_FULL}
\item \link{DBA_SCORE_TMM_MINUS_EFFECTIVE}
}
}
\item{mask}{
mask indicating a subset of peaksets to use when using global binding matrix scores. If a contrast is specified, these peaksets will be included, but only the significantly differentially bound sites (using th, bUsePval, and/or report) will be included. See  \code{\link{dba.mask}}.
}
 \item{sites}{
logical vector indicating which sites to include in PCA. Only relevant when using global binding matrix (contrast is missing).
}
 \item{cor}{
a logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. Passed into princomp.
}
\item{b3D}{
logical indicating that three principal components should be plotted (requires package\{rgl\}). If FALSE, the first two principal components are plotted.
}

\item{vColors}{
vector of custom colors; is absent, default colors will be used.
}
\item{dotSize}{
size of dots to plot; is absent, a default  will be calculated.
}

\item{\dots}{
arguments passed to plot or plot3d (rgl).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
MODE: PCA plot using significantly differentially bound sites:

dba.plotPCA(DBA, attributes, minval, maxval,
                       contrast, method, th, bUsePval,
                       b3D=F, vColors, dotSize, ...)

MODE: PCA plot using global binding matrix:

dba.plotPCA(DBA, attributes, minval, maxval,
                       mask, sites,
                       b3D=F, vColors, dotSize, ...)
}
\value{
matrix with color legend
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
\note{
uses rgl package for 3D plots (if available)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{dba.analyze}},  \code{\link{dba.plotHeatmap}}
}

\examples{
data(tamoxifen_peaks)

# peakcaller scores PCA
dba.plotPCA(tamoxifen)

# raw count correlation PCA
data(tamoxifen_analysis)
dba.plotPCA(tamoxifen)

#PCA based on normalized data for all sites
dba.plotPCA(tamoxifen,contrast=1,th=1)

#PCA based on DB sites only
par(mfrow=c(1,2))
dba.plotPCA(tamoxifen,contrast=1)
dba.plotPCA(tamoxifen,contrast=1,attributes=DBA_TISSUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
