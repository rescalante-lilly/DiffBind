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
dba:		\tab 	Construct a dba object\cr
dba.peakset:	\tab	Add a peakset to, or retrieve a peakset from, a dba object\cr
\cr

dba.overlap:	\tab	Compute binding site overlaps and/or correlations\cr
dba.count:	\tab	Count reads in binding sites\cr
\cr
	
dba.contrast:	\tab	Establish contrast(s) for analysis\cr
dba.analyze:	\tab	Execute affinity analysis\cr
dba.report:	\tab	Generate report for a contrast analysis\cr
\cr

dba.plotHeatmap:\tab	Heatmap plot\cr
dba.plotPCA:	\tab	Principal Components plot\cr
dba.plotBox:	\tab	Boxplots\cr
dba.plotMA:	\tab	MA/scatter plot\cr
dba.plotVenn:	\tab 	Venn diagram plot\cr
\cr
	
dba.show:	\tab	Show dba metadata\cr
dba.mask:	\tab	Mask samples or sites \cr
\cr
	
dba.save:	\tab	Save dba object\cr
dba.load:	\tab	Load dba object\cr
}
}
\author{
Rory Stark <rory.stark@cancer.org.uk> and Gordon Brown <gordon.brown@cancer.org.uk>
}
%\references{
%}
\keyword{ package }
%\seealso{
%}
%\examples{
%
%}
