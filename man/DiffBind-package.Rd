\name{DiffBind-package}
\alias{DiffBind-package}
\alias{DiffBind}
\docType{package}
\title{
DiffBind: Differential Binding Analysis
}
\description{
DiffBind -  differential binding analysis of ChIP-seq peaksets

}
\details{
\tabular{ll}{
Package: \tab DiffBind\cr
Type: \tab Package\cr
Version: \tab 0.99.7\cr
Date: \tab 2011-10-20\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}

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
Rory Stark and Gordon Brown
}
%\references{
%}
\keyword{ package }
%\seealso{
%}
%\examples{
%
%}
