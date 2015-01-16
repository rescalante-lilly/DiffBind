\name{DiffBind-package}
\alias{DiffBind-package}
\alias{DiffBind}
\docType{package}
\title{
Differential Binding Analysis of ChIP-seq peaksets
}
\description{
Differential binding analysis of ChIP-seq peaksets

}
\details{
Computes differentially bound sites from multiple ChIP-seq experiments using affinity (quantitative) data.
Also enables occupancy (overlap) analysis and plotting functions.


Entry Points:
\tabular{ll}{
\code{\link{dba}}:		\tab 	Construct a dba object\cr
\code{\link{dba.peakset}}:	\tab	Add a peakset to, or retrieve a peakset from, a dba object\cr
\cr

\code{\link{dba.overlap}}:	\tab	Compute binding site overlaps and/or correlations\cr
\code{\link{dba.count}}:	\tab	Count reads in binding sites\cr
\cr
	
\code{\link{dba.contrast}}:	\tab	Establish contrast(s) for analysis\cr
\code{\link{dba.analyze}}:	\tab	Execute affinity analysis\cr
\code{\link{dba.report}}:	\tab	Generate report for a contrast analysis\cr
\cr

\code{\link{dba.plotHeatmap}}:\tab	Heatmap plot\cr
\code{\link{dba.plotPCA}}:	\tab	Principal Components plot\cr
\code{\link{dba.plotBox}}:	\tab	Boxplots\cr
\code{\link{dba.plotMA}}:	\tab	MA/scatter plot\cr
\code{\link{dba.plotVenn}}:	\tab 	Venn diagram plot\cr
\cr
	
\code{\link{dba.show}}:	\tab	Show dba metadata\cr
\code{\link{dba.mask}}:	\tab	Mask samples or sites \cr
\cr
	
\code{\link{dba.save}}:	\tab	Save dba object\cr
\code{\link{dba.load}}:	\tab	Load dba object\cr
}
}
\author{
Rory Stark <rory.stark @at@ cruk.cam.ac.uk> and Gord Brown <gdbzork @at@ gmail.com>
}
%\references{
%}
\keyword{ package }
%\seealso{
%}
%\examples{
%
%}
