\name{dba}
\alias{dba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct a DBA object
}
\description{
Constructs a new DBA object from a sample sheet,  or based on an existing DBA object
}
\usage{
dba(DBA,mask, minOverlap=2,
            sampleSheet="dba_samples.csv", 
            config=data.frame(RunParallel=TRUE, reportInit="DBA", DataType=DBA_DATA_GRANGES, 
            AnalysisMethod=DBA_EDGER, bCorPlot=TRUE, th=.1, bUsePval=FALSE),
            peakCaller="raw", peakFormat, scoreCol, bLowerScoreBetter, filter, skipLines=0, 
            bAddCallerConsensus=FALSE, 
            bRemoveM=TRUE, bRemoveRandom=TRUE, bSummarizedExperiment=FALSE,
            bCorPlot, attributes) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Optional or defaulted arguments}{}

  \item{DBA}{
existing DBA object -- if present, will return a fully-constructed DBA object based on the passed one, using criteria specified in the mask and/or minOverlap parameters. If missing, will create a new DBA object based on the sampleSheet.
}
\item{mask}{
logical or numerical vector indicating which peaksets to include in the resulting model if basing DBA object on an existing one. See  \code{\link{dba.mask}}.
}
\item{minOverlap}{
only include peaks in at least this many peaksets in the main binding matrix if basing DBA object on an existing one.  If minOverlap is between zero and one, peak will be included from at least this proportion of peaksets.
}
\item{sampleSheet}{
data frame containing sample sheet, or file name of sample sheet to load (ignored if DBA is specified). Columns names in sample sheet may include:
\itemize{
\item {SampleID:}	{Identifier string for sample}
\item {Tissue:}	{Identifier string for tissue type}
\item {Factor:}	{Identifier string for factor}
\item {Condition:}	{Identifier string for condition}
\item {Treatment:}	{Identifier string for treatment}
\item {Replicate:}	{Replicate number of sample}
\item {bamReads:}	{file path for bam file containing aligned reads for ChIP sample}
\item {bamControl:}	{file path for bam file containing aligned reads for control sample}
\item {ControlID:}	{Identifier string for control sample}
\item {Peaks:}	{path for file containing peaks for sample. format determined by PeakCaller field or caller parameter}
\item {PeakCaller:}	{Identifier string for peak caller used. If Peaks is not a bed file, this will determine how the Peaks file is parsed. If missing, will use default peak caller specified in caller parameter. Possible values:}
\itemize{
\item {\dQuote{raw}:}	{text file file; peak score is in fourth column}
\item {\dQuote{bed}:}	{.bed file; peak score is in fifth column}
\item {\dQuote{narrow}:}  {default peak.format: narrowPeaks file}
\item {\dQuote{macs}:}	{MACS .xls file}
\item {\dQuote{swembl}:}	{SWEMBL .peaks file}
\item {\dQuote{bayes}:}	{bayesPeak file}
\item {\dQuote{peakset}:}	{peakset written out using pv.writepeakset}
\item {\dQuote{fp4}:}	{FindPeaks v4}
}
\item {PeakFormat:}	{string indicating format for peak files; see PeakCaller and  \code{\link{dba.peakset}}}
\item {ScoreCol:}	{column in peak files that contains peak scores}
\item {LowerBetter:}	{logical indicating that lower scores signify better peaks}
\item{Counts:} {file path for externally computed read counts; see \code{\link{dba.peakset}} (\code{counts} parameter)}
}
}
  \item{config}{
data frame containing configuration options, or file name of config file to load when constructing a new DBA object from a sample sheet. NULL indicates no config file. Relevant fields include:
\itemize{
\item {RunParallel:}	{logical indicating if counting and analysis operations should be run in parallel using multicore by default.}
\item {DataType:}	{default class for peaks and reports (DBA_DATA_GRANGES, DBA_DATA_RANGEDDATA, or DBA_DATA_FRAME).}
\item {ReportInit:}   {string to append to the beginning of saved report file names.}
\item {AnalysisMethod:}	{either DBA_EDGER or DBA_DESEQ.}
\item {bCorPlot:} {logical indicating that a correlation heatmap should be plotted automatically}
\item {th:}   {default threshold for reporting and plotting analysis results.}
\item {bUsePval:}   {logical, default indicating whether to use FDR (\code{FALSE}) or p-values (\code{TRUE}).}

}
}
  \item{peakCaller}{
if a sampleSheet is specified, the default peak caller that will be used if the \code{fPeakCaller} column is absent.
}
 \item{peakFormat}{
if a sampleSheet is specified, the default peak file format that will be used if the \code{fPeakFormat} column is absent.
}
 \item{scoreCol}{
if a sampleSheet is specified, the default column  in the peak files that will be used for scoring if the \code{ScoreCol} column is absent.
}
 \item{bLowerScoreBetter}{
if a sampleSheet is specified, the  sort order for peak scores if the \code{LowerBetter} column is absent.
}
 \item{filter}{
if a sampleSheet is specified, a filter value if the \code{Filter} column is absent. Peaks with scores lower than this value (or higher if \code{ bLowerScoreBetter} or \code{LowerBetter} is TRUE) will be removed.
}
  \item{skipLines}{
if a sampleSheet is specified, the number of lines (ie header lines) at the beginning of each peak file to skip.
}
  \item{bAddCallerConsensus}{
add a consensus peakset for each sample with more than one peakset (i.e. different peak callers) when constructing a new DBA object from a sample sheet.
}
  \item{bRemoveM}{
logical indicating whether to remove peaks on chrM (mitochondria) when constructing a new DBA object from a sample sheet.
}
\item{bRemoveRandom}{
logical indicating whether to remove peaks on chrN_random when constructing a new DBA object from a sample sheet.
}

\item{bSummarizedExperiment}{
logical indicating whether to return resulting object as a \code{\link{SummarizedExperiment}}.
}
\item{bCorPlot}{
logical indicating that a correlation heatmap should be plotted before returning. If \code{DBA} is \code{NULL} (a new DBA object is being created), and \code{bCorPlot} is missing, then this will take the default value (\code{TRUE}). However if \code{DBA} is \code{NULL} (a new DBA object is being created), and \code{bCorPlot} is specified, then the specified value will become the default value of \code{bCorPlot} for the resultant \code{DBA} object.
}
\item{attributes}{
vector of attributes to use subsequently as defaults when generating labels in plotting functions:
\itemize{ 
\item DBA_ID
\item DBA_TISSUE
\item DBA_FACTOR
\item DBA_CONDITION
\item DBA_REPLICATE
\item DBA_CONSENSUS
\item DBA_CALLER
\item DBA_CONTROL
}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

MODE: Construct a new DBA object from a samplesheet:

dba(sampleSheet, config,
    bAddCallerConsensus, bRemoveM, bRemoveRandom,
    attributes)

MODE: Construct a DBA object based on an existing one:

dba(DBA, mask, attributes)

MODE: Convert a DBA object to a SummarizedExperiment object:

dba(DBA, bSummarizedExperiment=TRUE)
}
\value{
DBA object
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
 \code{\link{dba.peakset}},  \code{\link{dba.show}}
}

\examples{
# Create DBA object from a samplesheet
setwd(system.file("extra", package="DiffBind"))
tamoxifen = dba(sampleSheet="tamoxifen.csv")
tamoxifen

tamoxifen = dba(sampleSheet="tamoxifen_allfields.csv")
tamoxifen

tamoxifen = dba(sampleSheet="tamoxifen_allfields.csv",config="config.csv")
tamoxifen

#Create a DBA object with a subset of samples
data(tamoxifen_peaks)
Responsive = dba(tamoxifen,tamoxifen$masks$Responsive)
Responsive

# change peak caller but leave peak format the same
setwd(system.file("extra", package="DiffBind"))
tamoxifen = dba(sampleSheet="tamoxifen.csv", peakCaller="macs", peakFormat="raw")
dba.show(tamoxifen, attributes=c(DBA_TISSUE,DBA_CONDITION,DBA_REPLICATE,DBA_CALLER))

# Convert DBA object to SummarizedExperiment
data(tamoxifen_counts)
sset = dba(tamoxifen,bSummarizedExperiment=TRUE)
sset
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
