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
        config=data.frame(RunParallel=TRUE, reportInit="DBA"),
        caller='raw', skipLines=0, bAddCallerConsensus=FALSE, 
        bRemoveM=TRUE, bRemoveRandom=TRUE, 
        bCorPlot=FALSE, attributes) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Optional or defaulted arguments}{}

  \item{DBA}{
existing DBA object -- if present, will return a fully-constructed DBA object based on the passed one, using criteria specified in the mask and/or minOverlap parameters.
}
\item{mask}{
logical or numerical vetcor indicating which peaksets to include in the resulting model if basing DBA object on an existing one. See dba.mask.
}
\item{minOverlap}{
only include peaks in at least this many peaksets in main matrix if basing DBA object on an existing one.
}
\item{sampleSheet}{
data frame containing sample sheet, or file name of sample sheet to load (ignored if DBA is specified). Columns names in sample sheet should include:
\itemize{
\item {SampleID:}	{Identifier string for sample}
\item {Tissue:}	{Identifier string for tissue type}
\item {Factor:}	{Identifier string for factor}
\item {Condition:}	{Identifier string for condition}
\item {Replicate:}	{Replicate number of sample}
\item {bamReads:}	{file path for bam file containing aligned reads for ChIP sample}
\item {bamControl:}	{file path for bam file containing aligned reads for control sample}
\item {ControlID:}	{Identifier string for control sample}
\item {Peaks:}	{path for file containing peaks for sample. format determined by PeakCaller field or caller parameter}
\item {PeakCaller:}	{Identifier string for peak caller used. If Peaks is not a bed file, this will determine how the Peaks file is parsed. Possible values:}
\itemize{
\item {\dQuote{raw}:}	{text file file; peak score is in fourth column}
\item {\dQuote{bed}:}	{.bed file; peak score is in fifth column}
\item {\dQuote{macs}:}	{MACS .xls file}
\item {\dQuote{swembl}:}	{SWEMBL .peaks file}
\item {\dQuote{bayes}:}	{bayesPeak file}
\item {\dQuote{peakset}:}	{peakset written out using pv.writepeakset}
\item {\dQuote{fp4}:}	{FindPeaks v4}
}
}
}
  \item{config}{
data frame containing sample sheet, or file name of config file to load when constructing a new DBA object from a sample sheet. NULL indicates no config file. Column names used for parallel operations:
\itemize{
\item {RunParallel:}	{logical indicating if counting and analysis operations should be run in parallel using multicore by default.}
\item {RangedData:}	{logical indicating if peaks should be RangedData objects by default.}
}
}
  \item{caller}{
if a sampleSheet is specified, the default peak file format that will be used id the PeakCaller column is absent.
}
  \item{skipLines}{
if a sampleSheet is specified, the number of lines (ie header lines)  at the top of each peak caller to skip.
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

\item{bCorPlot}{
logical indicating that a correlation heatmap should be plotted before returning
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
MODES:

Construct a new DBA object from a samplesheet:
dba(sampleSheet, config,
               bAddCallerConsensus, bRemoveM, bRemoveRandom,
               attributes)

Construct a DBA object based on an existing one:
dba(DBA, mask, attributes)
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# Create DBA object from a samplesheet
setwd(system.file("extra", package="DiffBind"))
tamoxifen = dba(sampleSheet="tamoxifen.csv")
tamoxifen

#Create a DBA object with a subset of samples
data(tamoxifen_peaks)
Responsive = dba(tamoxifen,tamoxifen$masks$Responsive)
Responsive
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
