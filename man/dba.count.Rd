\name{dba.count}
\alias{dba.count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count reads in binding site intervals
}
\description{
Counts reads in binding site intervals.  Files must be one of bam,
bed and gzip-compressed bed.  File suffixes must be ".bam", ".bed", or ".bed.gz"
respectively.
}
\usage{
dba.count(DBA, peaks, minOverlap=2, score=DBA_SCORE_TMM_MINUS_EFFECTIVE, bLog=FALSE,
          insertLength, maxFilter, bRemoveDuplicates=FALSE, bScaleControl=TRUE, 
          bCalledMasks=TRUE, bCorPlot=TRUE, bParallel=DBA$config$RunParallel) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


\item{DBA}{
DBA object
}
%\subsection{Optional/defaulted arguments}{}
 \item{peaks}{
GRanges, RangedData, dataframe, or matrix containing intervals to use. If missing, generates a consensus peakset using minOverlap parameter. If NULL, changes the score used in the global binding matrix to the score type specified in the score parameter.
}
\item{minOverlap}{
only include peaks in at least this many peaksets when generating consensus peakset (i.e. when peaks parameter is missing).  If minOverlap is between zero and one, peak will be included from at least this proportion of peaksets.
}
\item{score}{
which score to use in the binding affinity matrix. Note that all raw read counts are maintained for use by dba.analyze, regardless of how this is set. One of:
\tabular{ll}{
DBA_SCORE_READS \tab raw read count for interval using only reads from ChIP\cr
DBA_SCORE_READS_FOLD \tab raw read count for interval from ChIP divided by read count for interval from control\cr
DBA_SCORE_READS_MINUS \tab raw read count for interval from ChIP minus read count for interval from control\cr
DBA_SCORE_RPKM \tab RPKM for interval using only reads from ChIP\cr
DBA_SCORE_RPKM_FOLD \tab RPKM for interval from ChIP divided by RPKM for interval from control\cr
DBA_SCORE_TMM_READS_FULL \tab TMM normalized (using edgeR), using ChIP read counts and Full Library size\cr
DBA_SCORE_TMM_READS_EFFECTIVE \tab TMM normalized (using edgeR), using ChIP read counts and Effective Library size\cr
DBA_SCORE_TMM_MINUS_FULL \tab TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size\cr
DBA_SCORE_TMM_MINUS_EFFECTIVE \tab TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Effective Library size\cr
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
\item{maxFilter}{
value to use for filtering intervals with low read counts. 
Only intervals where at least one sample has at least maxFilter reads will be included. 
If missing, includes all intervals.
If peaks is NULL, will remove sites from existing DBA object without recounting.
}
\item{bRemoveDuplicates}{
logical indicating if duplicate reads (ones that map to exactly the same genomic position) should be removed. If TRUE, any location where multiple reads map will  be counted as a single read.
}
\item{bScaleControl}{
logical indicating if the Control reads should be scaled based on relative library sizes. If TRUE, and there are more reads in the Control library than in the ChIP library, the number of Control reads for each peak will be multiplied by a scaling factor determined by dividing the total number of reads in the ChIP library by the total number of reads in the Control library. If this value is not an integer, the number of Control reads for each peak  will be the next highest integer.
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


# Count using a peakset made up of only peaks in all responsive MCF7 replicates
data(tamoxifen_peaks)
mcf7Common = dba.overlap(tamoxifen,tamoxifen$masks$MCF7&tamoxifen$masks$Responsive)
\dontrun{tamoxifen = dba.count(tamoxifen,peaks=mcf7Common$inAll)}
tamoxifen

# Change binding affinity scores
data(tamoxifen_counts)
tamoxifen = dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_READS)
head(tamoxifen$vectors)
tamoxifen = dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_RPKM_FOLD)
head(tamoxifen$vectors)
tamoxifen = dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_TMM_MINUS_FULL)
head(tamoxifen$vectors)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
