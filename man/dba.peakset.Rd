\name{dba.peakset}
\alias{dba.peakset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add a peakset to, or retrieve a peakset from, a DBA object
}
\description{
Adds a peakset to, or retrieves a peakset from, a DBA object
}
\usage{
dba.peakset(DBA=NULL, peaks, sampID, tissue, factor, condition, replicate,
            control, peak.caller, reads=0, consensus=FALSE, 
            bamReads, bamControl,
            normCol=4, bRemoveM=TRUE, bRemoveRandom=TRUE,
            minOverlap=2, bMerge=TRUE,
            bRetrieve=FALSE, writeFile, numCols=4,
            DataType=DBA$config$DataType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


\item{DBA}{
DBA object. Required unless creating a new DBA object by adding an initial peakset.
}

\item{peaks}{
When adding a specified peakset: set of peaks, either a GRanges or RangedData object, or a peak dataframe or matrix (chr,start,end,score), or a filename where the peaks are stored.

When adding a consensus peakset: a sample mask or vector of peakset numbers. If missing or NULL, a consensus is derived from all peaksets present in the model. See dba.mask, or dba.show to get peakset numbers.

When adding all the peaks from one DBA object to another: a DBA object. In this case, the only other parameter to have an effect is minOverlap.

When retrieving and/or writing a peakset: either a GRanges or RangedData object, or a peak dataframe or matrix (chr,start,end,score), or a peakset number; if NULL, retrieves/writes the full binding matrix.
}

%\subsection{Optional/defaulted arguments}{}

\item{sampID}{
ID string for the peakset being added; if missing, one is assigned (a serial number for a new peakset, or a concatenation of IDs for a consensus peakset).
}
\item{tissue}{
tissue name for the peakset being added; if missing, one is assigned for a consensus peakset (a concatenation of tissues).
}
  \item{factor}{
factor name for the peakset being added; if missing, one is assigned for a consensus peakset (a concatenation of factors).
}
  \item{condition}{
condition name for the peakset being added; if missing, one is assigned for a consensus peakset (a concatenation of conditions).
}
  \item{replicate}{
replicate number for the peakset being added; if missing, one is assigned for a consensus peakset (a concatenation of replicate numbers).
}
  \item{control}{
control name for the peakset being added; if missing, one is assigned for a consensus peakset (a concatenation of control names).
}
\item{peak.caller}{
peak caller name string. If peaks is specified as a file, this will control how it is interpreted. Supported values:
\itemize{
\item {\dQuote{raw}:}	{text file file; peak score is in fourth column}
\item {\dQuote{bed}:}	{.bed file; peak score is in fifth column}
\item {\dQuote{macs}:}	{MACS .xls file}
\item {\dQuote{swembl}:}	{SWEMBL .peaks file}
\item {\dQuote{bayes}:}	{bayesPeak file}
\item {\dQuote{peakset}:}	{peakset written out using pv.writepeakset}
\item {\dQuote{fp4}:}	{FindPeaks v4}
}
if missing, a name is assigned for a consensus peakset (a concatenation of peak caller names).
}
\item{reads}{
total number of ChIPed library reads for the peakset being added.
}
  \item{consensus}{
TRUE if peakset being added is made from overlap of other peaksets (set automatically when adding a consensus peakset).
}
\item{bamReads}{
file path of the BAM/SAM/BED file containing the aligned reads for the peakset being added.
}
\item{bamControl}{
file path of the BAM/SAM/BED file containing the aligned reads for the control used for the peakset being added.
}
  \item{normCol}{
peak column to normalize to 0...1 scale when adding a peakset; 0 indicates no normalization
}
  \item{bRemoveM}{
logical indicating whether to remove peaks on chrM when adding a peakset
}
  \item{bRemoveRandom}{
logical indicating whether to remove peaks on chrN_random when adding a peakset
}
  \item{minOverlap}{
the minimum number of peaksets a peak must be in to be included when adding a consensus peakset.
}
  \item{bMerge}{
logical indicating whether global binding matrix should be compiled after adding the peakset. When adding several peaksets via successive calls to dba.peakset, it may be more efficient to set this parameter to FALSE and call dba(DBA) after all the peaksets have been added.
}
  \item{bRetrieve}{
logical indicating that a peakset is being retrieved and/or written, not added.
}
 \item{writeFile}{
file to write retrieved peakset.
}
 \item{numCols}{
number of columns to include when writing out peakset. First four columns are chr, start, end, score; the remainder are maintained from the original peakset. Ignored when writing out complete binding matrix.
}
\item{DataType}{
The class of object for returned peaksets:
\itemize{
\item DBA_DATA_GRANGES
\item DBA_DATA_RANGEDDATA
\item DBA_DATA_FRAME
}
Can be set as default behavior by setting DBA$config$DataType.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

MODE: Add a specified peakset:

dba.peakset(DBA=NULL, peaks, sampID, tissue, factor, condition, replicate,
                       control, peak.caller, reads, consensus, bamReads, bamControl,
                       normCol, bRemoveM, bRemoveRandom)

MODE: Add a consensus peakset (derived from overlapping peaks in peaksets already present):

dba.peakset(DBA, peaks, minOverlap)

MODE: Retrieve a peakset:

dba.peakset(DBA, peaks, bRetrieve=T)

MODE: Write a peakset out to a file:

dba.peakset(DBA, peaks, bRetrieve=T,  writeFile, numCols)
}
\value{
DBA object when adding a peakset.
Peakset matrix or RangedData object when retrieving and/or writing a peakset.
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

# create a new DBA object by adding three peaksets
mcf7 = dba.peakset(NULL,
                   peaks=system.file("extra/peaks/MCF7_ER_1.bed.gz", package="DiffBind"),
                   sampID="MCF7.1",tissue="MCF7",factor="ER",condition="Responsive",replicate=1)
mcf7 = dba.peakset(mcf7,
                   peaks=system.file("extra/peaks/MCF7_ER_2.bed.gz", package="DiffBind"),    
                   sampID="MCF7.2",tissue="MCF7",factor="ER",condition="Responsive",replicate=2)
mcf7 = dba.peakset(mcf7,
                   peaks=system.file("extra/peaks/MCF7_ER_3.bed.gz", package="DiffBind"),      
                   sampID="MCF7.3",tissue="MCF7",factor="ER",condition="Responsive",replicate=3)
mcf7

#add a consensus peakset -- peaks in all three replicates
mcf7 = dba.peakset(mcf7, 1:3, minOverlap=3,sampID="MCF7_3of3")
mcf7

#retrieve the consensus peakset as RangedData object
mcf7.consensus = dba.peakset(mcf7,mcf7$masks$Consensus,bRetrieve=TRUE)
mcf7.consensus

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
