\name{dba.show}
\alias{dba.show}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List attributes of peaksets of contrasts associated with a DBA object
}
\description{
Returns attributes of peaksets and/or contrasts associated with a DBA object.
}
\usage{
dba.show(DBA, mask, attributes, bContrasts=FALSE, th=0.1, bUsePval=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object
}
%\subsection{Optional/defaulted arguments}{}
  \item{mask}{
mask of peaksets for which to get attributes (used when obtaining peakset attributes, i.e. bContrasts=FALSE).
}
\item{attributes}{
attribute or vector of attributes to retrieve. Number of intervals is always shown. Used when obtaining peakset attributes, i.e. bContrasts=FALSE. Values:
\itemize{ 
\item DBA_ID
\item DBA_TISSUE
\item DBA_FACTOR
\item DBA_CONDITION
\item DBA_CONDITION
\item DBA_REPLICATE
\item DBA_CONSENSUS
\item DBA_CALLER
\item DBA_CONTROL
\item DBA_INTERVALS
\item DBA_SN_RATIO
}
}

 \item{bContrasts}{
logical indicating whether peaksets or contrast attributes are to be retrieved. TRUE retrieves a dataframe of contrast information instead of peakset attributes. If no contrasts are set, returns possible contrasts. See dba.contrast.
}
  \item{th}{
if bContrasts is TRUE, then th is used as the threshold for determining how many significant sites there are for each contrast. Only relevant when obtaining contrast attributes (bContrasts=TRUE) and dba.analyze has been run.
}
  \item{bUsePval}{
logical indicating that p-values will be used (along with th) to determine how many significant sites there are for each contrast; if FALSE, adjusted p-values (FDR) are used. Only relevant when obtaining contrast attributes (bContrasts=TRUE) and dba.analyze has been run.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
MODE: Return attributes of peaksets associated with a DBA object:

dba.show(DBA, mask, attributes)

MODE: Return contrasts associated with a DBA object:

dba.show(DBA,bContrasts=T, th, bUsePval)
}
\value{
dataframe with peakset attributes. 

If bContrasts == FALSE, each row represents a peakset, and each column is an attributes, with the final column, Intervals, indicating how many sites there are in the peakset.

If bContrasts == TRUE, each row represent a contrast, with the following columns:
\item{Group1}{Label for first group of contrast}
\item{Members1}{Number of samples in first group of contrast}
\item{Group2}{Label for first group of contrast}
\item{Members3}{Number of samples in first group of contrast}

if dba.analyze has been successfully run, there there will be up to two more columns showing the number of significant differentially bound (DB) sites identified for \item{DB.edgeR}{Number of significantly differentially bound sites identified using edgeR}
\item{DB.DESeq}{Number of significantly differentially bound sites identified using DESeq}

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
dba.show(tamoxifen)
dba.show(tamoxifen,tamoxifen$masks$Responsive)
dba.show(tamoxifen,attributes=c(DBA_TISSUE,DBA_REPLICATE,DBA_CONDITION))

data(tamoxifen_counts)
tamoxifen = dba.contrast(tamoxifen)
dba.show(tamoxifen,bContrasts=TRUE)

#alternatively:
tamoxifen
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
