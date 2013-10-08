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
dba.count(DBA, peaks, minOverlap=2, score=DBA_SCORE_TMM_MINUS_FULL, bLog=FALSE,
          insertLength, filter=0, bRemoveDuplicates=FALSE, bScaleControl=TRUE, 
          bCalledMasks=TRUE, filterFun=max, bCorPlot=DBA$config$bCorPlot, 
          bLowMem, bUseSummarizeOverlaps=FALSE, readFormat = DBA_READS_DEFAULT, 
          bParallel=DBA$config$RunParallel) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


 \item{DBA}{
DBA object
}
 \item{peaks}{
If GRanges, RangedData, dataframe, or matrix, this parameter contains the intervals to use for counting. If character string, it specifies a file containing the intervals to use (with the first three columns specifying chromosome, startpos, endpos).If missing or a mask, generates a consensus peakset using minOverlap parameter (after applying the mask if present). If NULL, changes the score used in the global binding matrix to the score type specified in the score parameter without re-counting.
}
\item{minOverlap}{
only include peaks in at least this many peaksets when generating consensus peakset (i.e. when peaks parameter is missing).  If minOverlap is between zero and one, peak will be included from at least this proportion of peaksets.
}
\item{score}{
which score to use in the binding affinity matrix. Note that all raw read counts are maintained for use by \code{\link{dba.analyze}}, regardless of how this is set. One of:
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
If missing, the read size indicated in the BAM/BED file will be used.
}
\item{filter}{
value to use for filtering intervals with low read counts. 
Only intervals where at least one sample has a \code{score} at least this high will be included. 
If \code{peaks} is NULL, will remove sites from existing DBA object without recounting. 
If filter is a vector of values, \code{dba.count} will return a vector of the same length, indicating how many intervals will be retained for each specified \code{filter} level.
}
\item{bRemoveDuplicates}{
logical indicating if duplicate reads (ones that map to exactly the same genomic position) should be removed. If TRUE, any location where multiple reads map will  be counted as a single read. Note that if \code{bLowMem} is set, duplicates needs to have been already marked in all of the BAM files.
The built-in counting code may not correctly handle certain cases when the \code{bRemoveDuplicates} paramter is set to \code{TRUE}. These cases include paried-end data and datasets where the read length may differ within a single BAM file. In these cases, see the \code{bUseSummarizeOverlaps} parameter.
}
\item{bScaleControl}{
logical indicating if the Control reads should be scaled based on relative library sizes. If TRUE, and there are more reads in the Control library than in the ChIP library, the number of Control reads for each peak will be multiplied by a scaling factor determined by dividing the total number of reads in the ChIP library by the total number of reads in the Control library. If this value is not an integer, the number of Control reads for each peak  will be the next highest integer.
}
\item{bCalledMasks}{
logical indicating whether to compute site masks for each peakset indicating which sites were originally identified as peaks (used by \code{\link{dba.report}}). 
}
\item{filterFun}{
function that will be invoked for each interval with a vector of scores for each sample. Returns a score that will be evaluated against the \code{filter} value (only intervals with a score at least as high as \code{filter} will be kept). Default is \code{max}, indicating that at least one sample should have a score of at least \code{filter}; other useful values include \code{sum} (indicating that all the scores added together should be at least \code{filter}) and \code{mean} (setting a minimum mean normalized count level). Users can supply their own function as well.
}
\item{bCorPlot}{
logical indicating whether to plot a correlation heatmap for the counted data
}
\item{bLowMem}{
DEPRECATED. See \code{bUseSummarizeOverlaps}.
}
\item{bUseSummarizeOverlaps}{
logical indicating that \code{\link{summarizeOverlaps}} should be used for counting instead of the built-in counting code. This option is slower but uses the more standard counting function. If \code{TRUE}, all read files must be BAM (.bam extension), with associated index files (.bam.bai extension). The \code{insertLength} parameter must absent. 

See notes for when the \code{bRemoveDuplicates} parameter is set to \code{TRUE}, where the built-in counting code may not correctly handle certain cases and \code{bUseSummarizeOverlaps} should be set to \code{TRUE}.

Four additional parameters for \code{\link{summarizeOverlaps}} may be specified in \code{DBA$config}:
\tabular{ll}{
DBA$config$yieldSize \tab yieldSize indicating how many reads to process at one time; default is 5000000. The lower this value, the less memory will be used, but the more time it will take to complete the count operation.\cr
DBA$config$intersectMode \tab mode indicating which overlap algorithm to use; default is \code{"IntersectionNotEmpty"}\cr 
DBA$config$singleEnd \tab logical indicating if reads are single end; default is \code{TRUE}\cr 
DBA$config$scanbamparam \tab \code{ScanBamParam} object to pass to \code{\link{summarizeOverlaps}}. If present, \code{bRemoveDuplicates} is ignored.\cr 
}
}
\item{readFormat}{
Specify the file type of the read files, over-riding the file extension. Possible values:
\tabular{ll}{
DBA_READS_DEFAULT \tab use file extension (.bam, .bed, .bed.gz) to determine file type\cr
DBA_READS_BAM \tab assume the file type is BAM, regardless of the file extension\cr
DBA_READS_BED \tab assume the file type is BED (or zipped BED), regardless of the file extension.
}
Note that if \code{readFormat} is anything other than \code{DBA_READS_DEFAULT}, all the read files must be of the same file type.
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
 
 \seealso{
 \code{\link{dba.analyze}}
 }

\examples{
# These won't run unless you have the reads available in a BAM or BED file
data(tamoxifen_peaks)
\dontrun{tamoxifen = dba.count(tamoxifen)}


# Count using a peakset made up of only peaks in all responsive MCF7 replicates
data(tamoxifen_peaks)
mcf7Common = dba.overlap(tamoxifen,tamoxifen$masks$MCF7&tamoxifen$masks$Responsive)
\dontrun{tamoxifen = dba.count(tamoxifen,peaks=mcf7Common$inAll)}
tamoxifen

#First make consensus peaksets from each set of replicates, then derive master consensus set for counting from those
data(tamoxifen_peaks)
tamoxifen = dba.peakset(tamoxifen,consensus = -DBA_REPLICATE)
\dontrun{tamoxifen = dba.count(tamoxifen, peaks=tamoxifen$masks$Consensus)}
tamoxifen

# Change binding affinity scores
data(tamoxifen_counts)
tamoxifen = dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_READS)
head(tamoxifen$vectors)
tamoxifen = dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_RPKM_FOLD)
head(tamoxifen$vectors)
tamoxifen = dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_TMM_MINUS_FULL)
head(tamoxifen$vectors)

# Plot effect of a range of filter values and then apply filter 
data(tamoxifen_counts)
rate.max = dba.count(tamoxifen, peaks=NULL, filter=0:250)
rate.sum = dba.count(tamoxifen, peaks=NULL, filter=0:250,filterFun=sum)
plot(0:250,rate.max/rate.max[1],type='l',xlab="Filter Value",ylab="Proportion Retained Sites")
lines(0:250,rate.sum/rate.sum[1],col=2)
tamoxifen = dba.count(tamoxifen,peaks=NULL,filter=125,filterFun=sum)
tamoxifen
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
