\name{dba.overlap}
\alias{dba.overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute binding site overlaps (occupancy analysis)
}
\description{
Computes binding overlaps and co-occupancy statistics 
}
\usage{
dba.overlap(DBA, mask, mode=DBA_OLAP_PEAKS, minVal=0,
            contrast, method=DBA$config$AnalysisMethod, th=.1, bUsePval=FALSE, report,
            byAttribute, bCorOnly=TRUE, CorMethod="pearson", 
            bRangedData=DBA$config$RangedData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\emph{Required arguments:}


\item{DBA}{
DBA object
}

%\emph{Optional/defaulted arguments:}
 
\item{mask}{
mask or vector of peakset numbers indicating a subset of peaksets to use (see dba.mask). When generating overlapping/unique peaksets, either two or three peaksets must be specified. If the mode type is DBA_OLAP_ALL, and a contrast is specified, a value of TRUE (mask=TRUE) indicates that all samples should be included (otherwise only those present in one of the contrast groups will be included).
}

\item{mode}{
indicates which results should be returned (see MODES below). One of:
\itemize{
\item DBA_OLAP_PEAKS
\item DBA_OLAP_ALL
\item DBA_OLAP_RATE
}
}

\item{minVal}{
minimum score value to be considered a "called" peak.
}
  \item{contrast}{
contrast number to use. Only specified if contrast data is to be used when mode=DBA_OLAP_ALL. See dba.show(DBA, bContrast=T) to get contrast numbers. 
}
  \item{method}{
if contrast is specified and mode=DBA_OLAP_ALL, use data from method used for analysis:
\itemize{ 
\item DBA_EDGER
\item DBA_DESEQ
\item DBA_EDGER_BLOCK
}
}
\item{th}{
if contrast is specified and mode=DBA_OLAP_ALL, significance threshold; all sites with FDR (or p-values, see bUsePval) less than or equal to this value will be included. A value of 1 will include all binding sites, but only the samples included in the contrast.
}
\item{bUsePval}{
if contrast is specified and mode=DBA_OLAP_ALL, logical indicating whether to use FDR (FALSE) or p-value (TRUE) for thresholding.
}
\item{report}{
if contrast is specified and mode=DBA_OLAP_ALL, a report (obtained from dba.report) specifying the data to be used . If this is present, the method, th, and bUsePval parameters are ignored.
}
\item{byAttribute}{
when computing co-occupancy statistics (DBA_OLAP_ALL), limit comparisons to peaksets with the same value for a specific attribute, one of:
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

\item{bCorOnly}{
when computing co-occupancy statistics (DBA_OLAP_ALL), logical indicating that only correlations, and not overlaps, should be computed. This is much faster if only correlations are desired (e.g. to plot the correlations using dba.plotHeatmap).
}

\item{CorMethod}{
when computing co-occupancy statistics (DBA_OLAP_ALL), method to use when computing correlations.
}

\item{bRangedData}{
logical indicating whether, if mode==DBA_OLAP_PEAKS, the peaksets should be returned as RangedData objects. Can be set as default behavior by setting DBA$config$RangedData=TRUE.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~

MODES:

Generate overlapping/unique peaksets:
dba.overlap(DBA, mask, mode=DBA_OLAP_PEAKS, minVal=0)

Compute correlation and co-occupancy statistics (e.g. for dba.plotHeatmap):
dba.overlap(DBA, mask, mode=DBA_OLAP_ALL, byAttribute, minVal=0,
            attributes = DBA$attributes, 
            bCorOnly=T,  
            CorMethod="pearson")

Compute correlation and co-occupancy statistics using significantly differentially bound sites (e.g. for dba.plotHeatmap):
dba.overlap(DBA, mask=T, mode=DBA_OLAP_ALL, byAttribute, minVal=0,
            contrast, method=DBA_EDGER, th=.1, bUsePval=F,
            attributes = DBA$attributes, bCorOnly=T, 
            CorMethod="pearson")

Compute overlap rates at different stringency thresholds:
dba.overlap(DBA, mask, mode=DBA_OLAP_RATE, minVal=0)
}
\value{
Value depends on the mode specified in the mode parameter.

If mode = DBA_OLAP_PEAKS, Value is an overlap record, a list of three peaksets for an A-B overlap, and seven peaksets for a A-B-C overlap:
\item{inAll}{peaks in all peaksets}
\item{onlyA}{peaks unique to peakset A}
\item{onlyB}{peaks unique to peakset B}
\item{onlyC}{peaks unique to peakset C}
\item{notA}{peaks in both peaksets B and C but not peakset A}
\item{notB}{peaks in both peaksets A and C but not peakset B}
\item{notC}{peaks in both peaksets A and B but not peakset C}


If mode = DBA_OLAP_ALL, Value is a correlation record, a matrix with a row for each pair of peaksets and the following columns:
\item{A}{peakset number of first peakset in overlap}
\item{B}{peakset number of second peakset in overlap}
\item{onlyA}{number of sites unique to peakset A}
\item{onlyB}{number of sites unique to peakset B}
\item{inAll}{number of peaks in both peakset A and B (merged)}
\item{R2}{correlation value A vs B}
\item{Overlap}{percentage overlap (number of overlapping sites divided by number of peaks unique to smaller peakset}


If mode = DBA_OLAP_RATE, Value is a vector whose length is the number of peaksets, containing the number of overlapping peaks at the corresponding minOverlaps threshold (i.e., Value[1] is the total number of unique sites, Value[2] is the number of unique sites appearing in at least two peaksets, Value[3] the number of sites overlapping in at least three peaksets, etc.).
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
# default mode: DBA_OLAP_PEAKS -- get overlapping/non overlapping peaksets
mcf7 = dba.overlap(tamoxifen,tamoxifen$masks$MCF7)
names(mcf7)
mcf7$inAll

# mode:  DBA_OLAP_ALL -- get correlation record
mcf7 = dba(tamoxifen,tamoxifen$masks$MCF7)
mcf7.corRec = dba.overlap(mcf7,mode=DBA_OLAP_ALL,bCorOnly=FALSE)
mcf7.corRec

# mode: DBA_OLAP_RATE -- get overlap rate vector
data(tamoxifen_peaks)
rate = dba.overlap(tamoxifen, mode=DBA_OLAP_RATE)
rate
plot(rate,type='b',xlab="# peaksets",ylab="# common peaks",main="Tamoxifen dataset overlap rate")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
