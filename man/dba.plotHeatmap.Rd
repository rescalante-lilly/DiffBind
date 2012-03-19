\name{dba.plotHeatmap}
\alias{dba.plotHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw a binding site heatmap
}
\description{
Draws a binding site heatmap
}
\usage{
dba.plotHeatmap(DBA, attributes=DBA$attributes, maxSites=1000, minval, maxval,
                contrast, method=DBA$config$AnalysisMethod, 
                th=.1, bUsePval=FALSE, report, score,
                mask, sites, sortFun,
                correlations=TRUE, olPlot=DBA_COR,
                margin=10, colScheme="Greens", distMethod="pearson",
                ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{

%\subsection{Required arguments}{}

\item{DBA}{
DBA object. 
}

%\subsection{Optional/defaulted arguments}{}

\item{attributes}{
attribute or vector of attributes to use for column labels:
\itemize{ 
\item DBA_ID
\item DBA_TISSUE
\item DBA_FACTOR
\item DBA_CONDITION
\item DBA_TREATMENT
\item DBA_REPLICATE
\item DBA_CONSENSUS
\item DBA_CALLER
\item DBA_CONTROL
}
}

\item{maxSites}{
maximum number of binding sites to use in heatmap. Only used when not drawing a correlation heatmap (correlations=FALSE)
}

\item{minval}{
Set all scores less than this to minval
}
\item{maxval}{
Set all scores greater than this to maxval
}

\item{contrast}{
number of contrast to report on; if present, draws a heatmap based on a differential binding affinity analysis (see dba.analyze). See dba.show(DBA, bContrast=T) to get contrast numbers.
}
\item{method}{
analysis method (used in conjunction with contrast):
\itemize{ 
\item DBA_EDGER
\item DBA_DESEQ
\item DBA_EDGER_BLOCK
\item DBA_DESEQ_BLOCK
}
}

\item{th}{
significance threshold; all sites with FDR (or p-values, see bUsePval) less than or equal to this value will be included in the report (subject to maxSites). Used in conjunction with contrast.
}
\item{bUsePval}{
logical indicating whether to use FDR (FALSE) or p-value (TRUE) for thresholding. Used in conjunction with contrast.
}
\item{report}{
report (obtained from dba.report) specifying the data to be used . If this is present, the method, th, and bUsePval parameters are ignored. Used in conjunction with contrast.
}
\item{score}{
Score to use for count data. Only used when plotting the global binding matrix (no contrast specified). One of:
\itemize{
\item DBA_SCORE_READS
\item DBA_SCORE_READS_MINUS
\item DBA_SCORE_READS_FOLD
\item DBA_SCORE_RPKM
\item DBA_SCORE_RPKM_FOLD
\item DBA_SCORE_TMM_READS_FULL
\item DBA_SCORE_TMM_READS_EFFECTIVE
\item DBA_SCORE_TMM_MINUS_FULL
\item DBA_SCORE_TMM_MINUS_EFFECTIVE
}
}

\item{mask}{
mask indicating a subset of peaksets to use when using global binding matrix (contrast is missing). See dba.mask.
}
  \item{sites}{
logical vector indicating which sites to include; first maxSites of these. Only relevant when using global binding matrix (contrast is missing).
}

 \item{sortFun}{
function taking a vector of scores and returning a single value. Only relevant when using global binding matrix (contrast is missing). If present, the global binding matrix will be sorted (descending) on the results, and the first maxSites used in the heatmap. Recommended sort function options include sd, mean, median, min.
}

\item{correlations}{
logical indicating that a correlation heatmap should be plotted (TRUE). If FALSE, a binding heatmap of scores/reads is plotted. This parameter can also be set to a correlation record; see dba.overlap(mode=DBA_OLAP_ALL), in which case a correlation heatmap is plotted based on the specified correlation record, using the statistic specified in olPlot.
}

 \item{olPlot}{
if correlations is specified as a dataframe returned by dba.overlap, indicates which statistic to plot. One of:
\itemize{
\item{DBA_COR}  {Correlation}
\item{DBA_OLAP}  {Percentage overlap}
\item{DBA_INALL}  {number of peaks common to both samples}
}
}
  \item{margin}{
margin size of plot
}
  \item{colScheme}{
Color scheme; see colorRampPalette {RColorBrewer} 
}
  \item{distMethod}{
distance method for clustering; see Dist {amap}.
}

\item{\dots}{
passed on to heatmap.2 (gplots), e.g. scale etc.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
MODE: Correlation Heatmap plot using statistics for global binding matrix:

dba.plotHeatmap(DBA, attributes=DBA$attributes, minval, maxval,
                      correlations, olPlot,
                      colScheme="Greens", distMethod="pearson", ...)
                      
MODE: Correlation Heatmap plot using statistics for significantly differentially bound sites:

dba.plotHeatmap(DBA, attributes=DBA$attributes, minval, maxval,
                      contrast, method=DBA_EDGER, th=.1, bUsePval=F,
                      overlaps, olPlot=DBA_COR,
                      colScheme="Greens", distMethod="pearson", ...)
                      
MODE: Binding heatmap plot using significantly differentially bound sites:

dba.plotHeatmap(DBA, attributes, maxSites, minval, maxval,
                      contrast, method, th, bUsePval, correlations=FALSE,
                      colScheme, distMethod, ...)


MODE: Binding heatmap plot using the global binding matrix:

dba.plotHeatmap(DBA, attributes, maxSites, minval, maxval,
                      mask, sites, correlations=FALSE, sortFun,
                      colScheme, distMethod, ...)
}
\value{
%%  ~Describe the value returned
if  correlations is not FALSE, the overlap/correlation matrix is returned.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(tamoxifen_peaks)
# peak overlap correlation heatmap
dba.plotHeatmap(tamoxifen)

data(tamoxifen_counts)
# counts correlation heatmap
dba.plotHeatmap(tamoxifen)

data(tamoxifen_analysis)
#correlation heatmap based on all normalized data
dba.plotHeatmap(tamoxifen,contrast=1,th=1)

#correlation heatmap based on DB sites only
dba.plotHeatmap(tamoxifen,contrast=1)

#binding heatmap based on DB sites
dba.plotHeatmap(tamoxifen,contrast=1,correlations=FALSE)

#binding heatmap based on 1,000 sites with highest variance
dba.plotHeatmap(tamoxifen,contrast=1,th=1,correlations=FALSE,sortFun=var)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
