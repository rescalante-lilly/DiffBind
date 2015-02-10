\name{dba.plotVenn}
\alias{dba.plotVenn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw 2-way, 3-way, or 4-way Venn diagrams of overlaps
}
\description{
Draws 2-way, 3-way, or 4-way Venn diagrams of overlaps
}
\usage{
dba.plotVenn(DBA, mask, overlaps, label1, label2, label3, label4, main, sub, 
             contrast, method=DBA$config$AnalysisMethod, 
             th=DBA$config$th, bUsePval=DBA$config$bUsePval,
             bDB=TRUE, bNotDB, bAll=TRUE, bGain=FALSE, bLoss=FALSE,
             labelAttributes, bReturnPeaksets=FALSE, DataType=DBA$config$DataType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}

\item{DBA}{
DBA object; if present, only the mask parameter will apply.
}

\item{mask}{
mask or vector of peakset numbers indicating which peaksets to include in Venn diagram. Only 2 or 3 peaksets should be included.  See  \code{\link{dba.mask}}.
Only one of mask or overlaps is used.
}

\item{overlaps}{
overlap record, as computed by  \code{\link{dba.overlap}(Report=DBA_OLAP_PEAKS)}. 
Only one of mask or overlaps is used.
}
%\subsection{Optional/defaulted arguments}{}

  \item{label1}{
label for first peakset in diagram
}
  \item{label2}{
label for second peakset in diagram
}
  \item{label3}{
label for third peakset in diagram
}
  \item{label4}{
label for fourth peakset in diagram
}
  \item{main}{
main title for plot
}
  \item{sub}{
subtitle for plot
}
  \item{contrast}{
contrast number(s) to use for results-based plots. This can be a vector of contrast numbers.  See \code{\link{dba.show}(DBA, bContrast=T)} to get contrast numbers.
}
 \item{method}{
if \code{contrast} is specified, include results from analyses using this method or methods:
\itemize{ 
\item \link{DBA_EDGER}
\item \link{DBA_DESEQ}
\item \link{DBA_DESEQ2}
\item \link{DBA_EDGER_BLOCK}
\item \link{DBA_DESEQ_BLOCK}
\item \link{DBA_DESEQ2_BLOCK}
\item \link{DBA_ALL_METHODS}
\item \link{DBA_ALL_BLOCK}
\item \link{DBA_ALL_METHODS_BLOCK}
}
}
\item{th}{
if \code{contrast} is specified, use this significance threshold; all sites with FDR (or p-values, see bUsePval) less than or equal to this value will be considered differentially bound (DB). 
}
\item{bUsePval}{
if \code{contrast} is specified, this logical indicates whether to use FDR (FALSE) or p-value (TRUE) for thresholding.
}
\item{bDB}{
if \code{contrast} is specified, this logical indicates that peaksets should include Differentially Bound (DB) sites (respecting the \code{th}, \code{bUsePval}, and \code{fold} parameters).
}
\item{bNotDB}{
if \code{contrast} is specified, this logical indicates that peaksets should include non-Differentially Bound (non-DB) sites (respecting the \code{th}, \code{bUsePval}, and \code{fold} parameters).
}
\item{bAll}{
if \code{contrast} is specified, this logical indicates peaksets combining peaks with both positive and negative fold changes should be included.
}
\item{bGain}{
if \code{contrast} is specified, this logical indicates that peaksets with only positive fold changes should be included.
}
\item{bLoss}{
if \code{contrast} is specified, this logical indicates that peaksets with only negative fold changes should be included.
}

\item{labelAttributes}{
is \code{label}s are not specified, use these attributes to create default labels:
\itemize{ 
\item \link{DBA_ID}
\item \link{DBA_TISSUE}
\item \link{DBA_FACTOR}
\item \link{DBA_CONDITION}
\item \link{DBA_TREATMENT}
\item \link{DBA_REPLICATE}
\item \link{DBA_CONSENSUS}
\item \link{DBA_CALLER}
}
Only specified attributes that differ between peaksets will be used for labels; the ones that have the same value for all peaksets will be used as the default subtitle.
}
\item{bReturnPeaksets}{
logical indicating that the peaksets that make up portions of the Venn diagram should be returned. See \code{\link{dba.overlap}} for the format of the returned object. 
If the \code{DataType} parameter is set to \code{\link{DBA_DATA_DBAOBJECT}}, and a \code{contrast} has been spceified, a result-based DBA object will be returned, as described in \code{\link{dba.report}}.
}
\item{DataType}{
if \code{bReturnPeaksets} is set to \code{TRUE}, the class of object that peaksets should be returned as:
\itemize{
\item \link{DBA_DATA_GRANGES}
\item \link{DBA_DATA_RANGEDDATA}
\item \link{DBA_DATA_FRAME}
}
Can be set as default behavior by setting DBA$config$DataType.

Alternatively, this can be set to:
\itemize{
\item \link{DBA_DATA_DBAOBJECT}
}
to return a results-based DBA object, if a \code{contrast} is specified.
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
None, unless \code{bReturnPeaksets = TRUE}. In this case, either a list of peaksets is returned (as decribed in \code{\link{dba.overlap}}), or, if \code{DataType=DBA_DATA_DBAOBJECT}, a results-based DBA object.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
\note{

When working with results overlaps (a least one \code{contrast} is specified), and results-oriented DBA object is generated internally (as decribed in \code{\link{dba.report}}). In some cases, it may be better to generate the DBA object explicitly (using \code{\link{dba.report}} or setting \code{bReturnPeaksets=TRUE} and \code{DataType=DBA_DATA_DBAOBJECT}). This include the case where mseveral plots are being made of the same results set, and it takes a long time to generate the results-based DBA object, as well as the case where there are more than four results peaksets and a mask needs to be specified. I


This function relies on  \code{\link{vennPlot}} in the  \code{systemPipeR} package, written by Thomas Girke.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{dba.analyze}},  \code{\link{dba.overlap}},  
 \code{\link{dba.report}}, \code{\link{dba.plotPCA}},
 \code{\link{vennPlot}}
}

\examples{
data(tamoxifen_peaks)

par(mfrow=c(2,2))
# 2-way Venn
dba.plotVenn(tamoxifen,6:7)
dba.plotVenn(tamoxifen,tamoxifen$masks$ZR75)

# 3-way Venn (done two different ways)
dba.plotVenn(tamoxifen,tamoxifen$masks$MCF7&tamoxifen$masks$Responsive)
olaps = dba.overlap(tamoxifen,tamoxifen$masks$MCF7&tamoxifen$masks$Responsive)
dba.plotVenn(tamoxifen,overlaps=olaps,
             label1="Rep 1",label2="Rep 2",label3="Rep 3",main="MCF7 (Responsive) Replicates")

#Venn of overlaps
Responsive=dba(tamoxifen,tamoxifen$masks$Responsive)
Responsive
Responsive = dba.peakset(Responsive,1:3,sampID="MCF7")
Responsive = dba.peakset(Responsive,4:5,sampID="T47D")
Responsive = dba.peakset(Responsive,6:7,sampID="ZR75")
par(mfrow=c(1,1))
dba.plotVenn(Responsive,Responsive$masks$Consensus)

#4-way overlap
data(tamoxifen_peaks)
tamoxifen = dba.peakset(tamoxifen, consensus=DBA_TISSUE)
par(mfrow=c(1,1))
dba.plotVenn(tamoxifen,tamoxifen$masks$Consensus,main="Tissue consensus overlaps")

#Venns of differentially bound sites
data(tamoxifen_analysis)
tamoxifen = dba.contrast(tamoxifen,categories=DBA_CONDITION,block=tamoxifen$masks$MCF7)
tamoxifen = dba.analyze(tamoxifen,method=c(DBA_EDGER,DBA_DESEQ2))
dba.plotVenn(tamoxifen,contrast=1,method=DBA_ALL_METHODS_BLOCK)
dba.plotVenn(tamoxifen,contrast=1,method=DBA_ALL_BLOCK,bAll=FALSE,bGain=TRUE,bLoss=TRUE)
par(mfrow=c(2,1))
dba.plotVenn(tamoxifen,contrast=1,method=DBA_ALL_BLOCK,bAll=FALSE,bGain=TRUE,bLoss=FALSE)
dba.plotVenn(tamoxifen,contrast=1,method=DBA_ALL_BLOCK,bAll=FALSE,bGain=FALSE,bLoss=TRUE)
}

