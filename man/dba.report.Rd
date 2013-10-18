\name{dba.report}
\alias{dba.report}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a report for a differential binding affinity analysis
}
\description{
Generates a report for a differential binding affinity analysis
}
\usage{
dba.report(DBA, contrast, method=DBA$config$AnalysisMethod, 
           th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
           fold=0, bNormalized=TRUE,
           bCalled=FALSE, bCounts=FALSE, bCalledDetail=FALSE,
           bDB, bNotDB, bAll=TRUE, bGain=FALSE, bLoss=FALSE,
           file,initString=DBA$config$reportInit,ext='csv',DataType=DBA$config$DataType) 
}        
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


\item{DBA}{
DBA object. A differential binding affinity analysis needs to have been previously carried out (see  \code{\link{dba.analyze}}).
}
%\subsection{Optional/defaulted arguments}{}
  \item{contrast}{
contrast number to report on. When generating a report-based DBA object, this can be a vector of contrast numbers. If missing, defaults to first contrast for reports, and all contrasts when generating a report-based DBA object.  See  \code{\link{dba.show}(DBA, bContrast=T)} to get contrast numbers.
}
  \item{method}{
method used for analysis:
\itemize{ 
\item \link{DBA_EDGER}
\item \link{DBA_DESEQ}
\item \link{DBA_DESEQ2}
\item \link{DBA_EDGER_BLOCK}
\item \link{DBA_DESEQ_BLOCK}
\item \link{DBA_DESEQ2_BLOCK}
}
When generating a report-based DBA object (see \code{bDB} and \code{bNotDB} parameters below), a list of methods may be supplied, including the shortcuts
\itemize{ 
\item \link{DBA_ALL_METHODS}
\item \link{DBA_ALL_BLOCK}
\item \link{DBA_ALL_METHODS_BLOCK}
}
}
\item{th}{
significance threshold; all sites with FDR (or p-values, see bUsePval) less than or equal to this value will be included in the report. A value of 1 will include all binding sites in the report.
}
\item{bUsePval}{
logical indicating whether to use FDR (FALSE) or p-value (TRUE) for thresholding.
}
\item{fold}{
only sites with an absolute Fold value greater than equal to this will be included in the report.
}
\item{bNormalized}{
logical indicating that normalized data (using normalization factors computed by differential analysis method) should be reported. FALSE uses raw count data.
}
\item{bCalled}{
logical indicating that peak caller status should be included (if available from a previous call to  \code{\link{dba.count}(bCalledMasks=TRUE)}). This will add a column for each group, each indicating the number of samples in the group identified as a peak in the original peaksets.
}
\item{bCounts}{
logical indicating that count data for individual samples should be reported as well as group statistics. Columns are added for each sample in the first group, followed by columns for each sample in the second group.
}
\item{bCalledDetail}{
logical indicating that peak caller status should be included for each sample (if available).  Columns are added for each sample in the first group, followed by columns for each sample in the second group.
}
\item{bDB}{
logical indicating that a report-based DBA object should be generated, and that it should include Differentially Bound (DB) sites (respecting the \code{th}, \code{bUsePval}, and \code{fold} parameters).
}
\item{bNotDB}{
logical indicating that a report-based DBA object should be generated, and that it should include non-Differentially Bound (non-DB) sites (respecting the \code{th}, \code{bUsePval}, and \code{fold} parameters).
}
\item{bAll}{
logical indicating that a report-based DBA object should be generated, and that it should include peaksets combining peaks with both positive and negative fold changes.
}
\item{bGain}{
logical indicating that a report-based DBA object should be generated, and that it should include peaksets with only positive fold changes.
}
\item{bLoss}{
logical indicating that a report-based DBA object should be generated, and that it should include peaksets with only negative fold changes.
}
  \item{file}{
if present, also save the report to a comma separated value (csv) file, using this filename.
}
  \item{initString}{
if saving to a file, pre-pend this string to the filename.
}
  \item{ext}{
if saving to a file, append this extension to the filename.
}
\item{DataType}{
The class of object for returned report:
\itemize{
\item \link{DBA_DATA_GRANGES}
\item \link{DBA_DATA_RANGEDDATA}
\item \link{DBA_DATA_FRAME}
}
If set to \code{\link{DBA_DATA_SUMMARIZED_EXPERIMENT}}, the result will be a \code{\link{SummarizedExperiment}} object, with all the count data and sample metadata for the experiment. The contrast statistics will be included as metadata columns in the \code{rowData} of the object.

Can be set as default behavior by setting DBA$config$DataType.
}
}
%\details{
%}
\value{
if neither \code{bDB} or \code{bNotDB} is set to \code{TRUE}, a report dataframe, GRanges, or RangedData object is returned, with a row for each binding site within the thresholding parameters, and the following columns:
\item{Chr}{Chromosome of binding site}
\item{Start}{Starting base position of binding site}
\item{End}{End base position of binding site}
\item{Conc}{Concentration -- mean (log) reads across all samples in both groups}
\item{Conc_group1}{Group 1 Concentration -- mean (log) reads across all samples first group}
\item{Conc_group2}{Group 2 Concentration -- mean (log) reads across all samples in second group}
\item{Fold}{Fold difference -- mean fold difference of binding affinity of group 1 over group 2 (Conc1 - Conc2). Absolute value indicates magnitude of the difference, and sign indicates which one is bound with higher affinity, with a positive value indicating higher affinity in the first group }
\item{p-value}{p-value calculation -- statistic indicating significance of difference (likelihood difference is not attributable to chance)}
\item{FDR}{adjusted p-value calculation -- p-value subjected to multiple-testing correction}

If bCalled is TRUE and caller status is available, two more columns will follow:
\item{Called1}{Number of samples in group 1 that identified this binding site as a peak}
\item{Called2}{Number of samples in group 2 that identified this binding site as a peak}

If bCounts is TRUE, a column will be present for each sample in group 1, followed by each sample in group 2. The Sample ID will be used as the column header. This column contains the read counts for the sample.

If bCalledDetail is TRUE, a column will be present for each sample in group 1, followed by each sample in group 2. The Sample ID will be used as the column header. This column contains a "+" to indicate for which sites the sample was called as a peak, and a "-" if it was not so identified.

If  \code{bDB} or \code{bNotDB} is set to \code{TRUE}, a special \code{DBA object} is returned,containing peaksets based on sites determined to be differentiall bound (or not) as spcified using the \code{bDB}, \code{bNotDB}, \code{bGain}, \code{bLoss}, and \code{bAll} parameters. In this \code{DBA object}, the \code{Tissue} value will specify the direction of the change (\code{Gain} for positive fold changes, \code{Loss} for negative fold changes, and \code{All} for any fold change). The \code{Factor} value specifies if the peaks are differnetially bound (\code{DB}) or not (\code{!DB}). The \code{Condition} value specifies the analysis method (e.g. \code{edgeR}), and the \code{Treatment} value is blank for unblocked analyses and set to \code{block} for blocked analyses.
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

\seealso{
 \code{\link{dba.analyze}}
}

\examples{
data(tamoxifen_analysis)

#Retrieve DB sites with FDR <  0.1
tamoxifen.DB = dba.report(tamoxifen)
tamoxifen.DB

#Retrieve DB sites with p-value <  0.05 and Fold > 2
tamoxifen.DB = dba.report(tamoxifen, th=.05, bUsePval=TRUE, fold=2)
tamoxifen.DB

#Retrieve all sites with confidence stats
# and how many times each site was identified as a peak
tamoxifen.DB = dba.report(tamoxifen, th=1, bCalled=TRUE)
tamoxifen.DB

#Retrieve all sites with confidence stats and normalized counts
tamoxifen.DB = dba.report(tamoxifen, th=1, bCounts=TRUE)
tamoxifen.DB

#Retrieve all sites with confidence stats and raw counts
tamoxifen.DB = dba.report(tamoxifen, th=1, bCounts=TRUE,bNormalized=FALSE)
tamoxifen.DB

#Retrieve report as a SummarizedObject
tamoxifen.sset = dba.report(tamoxifen, DataType=DBA_DATA_SUMMARIZED_EXPERIMENT)
tamoxifen.sset

#Retrieve report-based DBA object
data(tamoxifen_counts)
tamoxifen = dba.contrast(tamoxifen, categories=DBA_CONDITION, block=tamoxifen$masks$MCF7)
tamoxifen = dba.analyze(tamoxifen,bCorPlot=FALSE)
tamoxifen.DB = dba.report(tamoxifen,method=c(DBA_EDGER,DBA_EDGER_BLOCK),
                          bDB=TRUE, bGain=TRUE, bLoss=TRUE)
dba.plotVenn(tamoxifen.DB,1:4,label1="Single Factor GAIN",label2="Single Factor LOSS",
                             label3="Blocking Factor GAIN",label4="Blocking Factor LOSS")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
