\name{DiffBind -- DBA global constant variables}
\alias{DBA_ID}
\alias{DBA_FACTOR}
\alias{DBA_TISSUE}
\alias{DBA_CONDITION}
\alias{DBA_REPLICATE}
\alias{DBA_CALLER}
\alias{DBA_CONSENSUS}
\alias{DBA_CONTROL}
\alias{DBA_GROUP}
\alias{DBA_OLAP_PEAKS}
\alias{DBA_OLAP_ALL}
\alias{DBA_OLAP_RATE}
\alias{DBA_SCORE_READS}
\alias{DBA_SCORE_READS_MINUS}
\alias{DBA_SCORE_RPKM}
\alias{DBA_SCORE_READS_FOLD}
\alias{DBA_SCORE_RPKM_FOLD}
\alias{DBA_SCORE_RPKM_MINUS}
\alias{DBA_EDGER}
\alias{DBA_DESEQ}
\alias{DBA_EDGER_BLOCK}
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
DBA_REPLICATE
DBA_CALLER
DBA_CONSENSUS
DBA_CONTROL

DBA_GROUP

DBA_OLAP_PEAKS
DBA_OLAP_ALL
DBA_OLAP_RATE

DBA_SCORE_READS
DBA_SCORE_RPKM
DBA_SCORE_READS_FOLD
DBA_SCORE_READS_MINUS
DBA_SCORE_RPKM_FOLD
DBA_EDGER
DBA_DESEQ
DBA_EDGER_BLOCK
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
dba.count score should be number of reads
}
\item{DBA_SCORE_RPKM}{
dba.count score should be RPKM
}
\item{DBA_SCORE_READS_FOLD}{
dba.count score should be number of reads divided by number of reads in control
}
\item{DBA_SCORE_READS_MINUS}{
dba.count score should be number of reads minus  number of reads in control
}
\item{DBA_SCORE_RPKM_FOLD}{
dba.count score should be RPKM divided by RPKM of control
}
\item{DBA_EDGER}{
differential analysis method: edgeR
}
\item{DBA_DESEQ}{
differential analysis method: DESeq
}
\item{DBA_EDGER_BLOCK}{
differential analysis method: edgeR with blocking factor
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
