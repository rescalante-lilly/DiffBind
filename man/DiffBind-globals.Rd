\name{DiffBind -- DBA global constant variables}
\alias{DBA_ID}
\alias{DBA_FACTOR}
\alias{DBA_TISSUE}
\alias{DBA_CONDITION}
\alias{DBA_TREATMENT}
\alias{DBA_REPLICATE}
\alias{DBA_CALLER}
\alias{DBA_CONSENSUS}
\alias{DBA_CONTROL}
\alias{DBA_INTERVALS}
\alias{DBA_SN_RATIO}
\alias{DBA_GROUP}
\alias{DBA_OLAP_PEAKS}
\alias{DBA_OLAP_ALL}
\alias{DBA_OLAP_RATE}
\alias{DBA_SCORE_READS}
\alias{DBA_SCORE_READS_MINUS}
\alias{DBA_SCORE_READS_FOLD}
\alias{DBA_SCORE_RPKM}
\alias{DBA_SCORE_RPKM_FOLD}
\alias{DBA_SCORE_TMM_READS_FULL}
\alias{DBA_SCORE_TMM_READS_EFFECTIVE}
\alias{DBA_SCORE_TMM_MINUS_FULL}
\alias{DBA_SCORE_TMM_MINUS_EFFECTIVE}
\alias{DBA_EDGER}
\alias{DBA_DESEQ}
\alias{DBA_EDGER_BLOCK}
\alias{DBA_DESEQ_BLOCK}
\alias{DBA_EDGER_CLASSIC}
\alias{DBA_DESEQ_CLASSIC}
\alias{DBA_EDGER_GLM}
\alias{DBA_DESEQ_GLM}
\alias{DBA_DATA_FRAME}
\alias{DBA_DATA_GRANGES}
\alias{DBA_DATA_RANGEDDATA}
\alias{DBA_DATA_SUMMARIZED_EXPERIMENT}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constant variables used in DiffBind package
}
\description{
Constant variables used in DiffBind package
}
\usage{
DBA_ID
DBA_FACTOR
DBA_TISSUE
DBA_CONDITION
DBA_TREATMENT
DBA_REPLICATE
DBA_CALLER
DBA_CONSENSUS
DBA_CONTROL

DBA_INTERVALS
DBA_SN_RATIO

DBA_GROUP

DBA_OLAP_PEAKS
DBA_OLAP_ALL
DBA_OLAP_RATE

DBA_SCORE_READS
DBA_SCORE_READS_MINUS
DBA_SCORE_READS_FOLD
DBA_SCORE_RPKM
DBA_SCORE_RPKM_FOLD
DBA_SCORE_TMM_READS_FULL
DBA_SCORE_TMM_READS_EFFECTIVE
DBA_SCORE_TMM_MINUS_FULL
DBA_SCORE_TMM_MINUS_EFFECTIVE

DBA_EDGER
DBA_DESEQ
DBA_EDGER_BLOCK
DBA_DESEQ_BLOCK
DBA_EDGER_CLASSIC
DBA_DESEQ_CLASSIC
DBA_EDGER_GLM
DBA_DESEQ_GLM

DBA_DATA_FRAME
DBA_DATA_GRANGES
DBA_DATA_RANGEDDATA
DBA_DATA_SUMMARIZED_EXPERIMENT
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{DBA_ID}{
DBA peakset metadata: Peakset ID
}
\item{DBA_FACTOR}{
DBA peakset metadata: Factor
}
\item{DBA_TISSUE}{
DBA peakset metadata: Tissue
}
\item{DBA_CONDITION}{
DBA peakset metadata: Condition
}
\item{DBA_TREATMENT}{
DBA peakset metadata: Treatment
}
\item{DBA_REPLICATE}{
DBA peakset metadata: Replicate
}
\item{DBA_CALLER}{
DBA peakset metadata: Peak Caller
}
\item{DBA_CONSENSUS}{
DBA peakset metadata: Is this a consensus peakset?
}
\item{DBA_CONTROL}{
DBA peakset metadata: ID of Control sample
}

\item{DBA_INTERVALS}{
DBA peakset metadata: Number of intervals in peakset
}
\item{DBA_SN_RATIO}{
DBA peakset metadata: Signal to Noise ratio (number of reads in intervals divided by total number of reads in library)
}

\item{DBA_GROUP}{
DBA peakset metadata: color PCA plot using contras groups
}

\item{DBA_OLAP_PEAKS}{
dba.overlap mode: return overlapping/unique peaksets
}
\item{DBA_OLAP_ALL}{
dba.overlap mode: return report of correlations/overlaps for each pair of samples
}
\item{DBA_OLAP_RATE}{
dba.overlap mode: return overlap rates
}

\item{DBA_SCORE_READS}{
dba.count score is number of reads in ChIP
}
\item{DBA_SCORE_READS_FOLD}{
dba.count score is number of reads in ChIP divided by number of reads in Control
}
\item{DBA_SCORE_READS_MINUS}{
dba.count score is number of reads in ChIP minus  number of reads in Control
}
\item{DBA_SCORE_RPKM}{
dba.count score is RPKM of ChIP 
}
\item{DBA_SCORE_RPKM_FOLD}{
dba.count score is RPKM of ChIP divided by RPKM of Control
}
\item{DBA_SCORE_TMM_READS_FULL}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts and Full Library size
}
\item{DBA_SCORE_TMM_READS_EFFECTIVE}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts and Effective Library size
}
\item{DBA_SCORE_TMM_MINUS_FULL}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size
}
\item{DBA_SCORE_TMM_MINUS_EFFECTIVE}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Effective Library size
}

\item{DBA_EDGER}{
differential analysis method: edgeR (default: DBA_EDGER_GLM)
}
\item{DBA_DESEQ}{
differential analysis method: DESeq (default: DBA_DESEQ_CLASSIC)
}

\item{DBA_EDGER_CLASSIC}{
differential analysis method: "classic" edgeR for two-group comparisons
}
\item{DBA_DESEQ_CLASSIC}{
differential analysis method: "classic" DESeq for two-group comparisons
}
\item{DBA_EDGER_GLM}{
differential analysis method: use GLM in edgeR for two-group comparisons
}
\item{DBA_DESEQ_GLM}{
differential analysis method: use GLM in DESeq for two-group comparisons
}
\item{DBA_EDGER_BLOCK}{
differential analysis method: edgeR with blocking factors (GLM)
}
\item{DBA_DESEQ_BLOCK}{
differential analysis method: DESeq with blocking factors (GLM)
}

\item{DBA_DATA_GRANGES}{
Use GRanges class for peaksets and reports. This is the  default (DBA$config$DataType = DBA_DATA_GRANGES).
}
\item{DBA_DATA_RANGEDDATA}{
Use RangedData class for peaksets and reports. Can be set as default (DBA$config$DataType = DBA_DATA_RANGEDDATA).
}
\item{DBA_DATA_FRAME}{
Use data.frame class for peaksets and reports. Can be set as default (DBA$config$DataType = DBA_DATA_FRAME).
}
\item{DBA_DATA_SUMMARIZED_EXPERIMENT}{
Used in \code{\link{dba.peakset}} to convert a \code{DBA} object to a \code{\link{SummarizedExperiment}} object.
}
}
%\details{
%%
%}
%\value{
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
\note{
Variables with ALL CAP names are used as constants within DiffBind.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
