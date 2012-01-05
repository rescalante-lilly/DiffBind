\name{dba.count}
\alias{dba.count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count reads in binding site intervals
}
\description{
Counts reads in binding site intervals
}
\usage{
dba.count(DBA, peaks, minOverlap=2, score=DBA_SCORE_READS_MINUS, bLog=FALSE,
          insertLength, minMaxval, bCalledMasks=TRUE, 
          bCorPlot=TRUE, bParallel=DBA$config$RunParallel) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object
}
%\subsection{Optional/defaulted arguments}{}
 \item{peaks}{
RangedData or matrix containing intervals to use. If missing, generates a consensus peakset using minOverlap parameter.
}
\item{minOverlap}{
only include peaks in at least this many peaksets when generating consensus peakset (i.e. when peaks parameter is missing).
}
\item{score}{
which score to use in the binding affinity matrix. Note that all raw read counts are maintained for use by dba.analyze, regardless of how this is set. One of:
\tabular{ll}{
DBA_SCORE_RPKM \tab RPKM for interval using only reads from ChIP\cr
DBA_SCORE_RPKM_FOLD \tab RPKM for interval from ChIP divided by RPKM for interval from control\cr
DBA_SCORE_READS \tab raw read count for interval using only reads from ChIP\cr
DBA_SCORE_READS_FOLD \tab raw read count for interval from ChIP divided by read count for interval from control\cr
DBA_SCORE_READS_MINUS \tab raw read count for interval from ChIP minus read count for interval from control\cr
}
}
\item{bLog}{
logical indicating whether log2 of score should be used (only applies to DBA_SCORE_RPKM_FOLD and DBA_SCORE_READS_FOLD).
}
\item{insertLength}{
if present, this value will be used as the length of the reads. 
Each read will be extended from its endpoint along the appropriate strand by this many bases. 
If missing, the read size indicated in the BAM/SAM/BED file will be used.
}
\item{minMaxval}{
value to use for filtering intervals with low read counts. 
Only intervals where at least one sample has at least minMaxval reads will be included. 
If missing, includes all intervals.
If peaks is NULL, will remove sites from existing DBA object without recounting.
}
\item{bCalledMasks}{
logical indicating whether to compute site masks for each peakset indicating which sites were originally identified as peaks(used by dba.report). 
}
\item{bCorPlot}{
logical indicating whether to plot a correlation heatmap for the counted data
}

\item{bParallel}{
if TRUE, use multicore to get counts for each read file in parallel
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
DBA object with binding affinity matrix based on read count scores.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark and Gordon Brown
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# These won't run unless you have the reads available in a BAM, SAM, or BED file
data(tamoxifen_peaks)
\dontrun{tamoxifen = dba.count(tamoxifen)}


# Count using a peakset made up of only peaks in all MCF7 replicates
data(tamoxifen_peaks)
mcf7Common = dba.overlap(tamoxifen,tamoxifen$masks$MCF7)
\dontrun{tamoxifen = dba.count(tamoxifen,peaks=mcf7Common$inAll)}
tamoxifen

# Change binding affinity scores
data(tamoxifen_counts)
tamoxifen = dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_READS)
head(tamoxifen$vectors)
tamoxifen = dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_RPKM_FOLD)
head(tamoxifen$vectors)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
