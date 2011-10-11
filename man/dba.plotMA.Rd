\name{dba.plotMA}
\alias{dba.plotMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate MA and scatter plots of differential binding analysis results
}
\description{
Generates MA and scatter plots of differential binding analysis results.
}
\usage{
dba.plotMA(DBA, contrast=1, method=DBA_EDGER, th=.1, bUsePval=FALSE, 
                      bNormalized=TRUE, factor="", bXY=FALSE, dotSize=.33, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object, on which dba.analyze should have been successfully run.
}
%\subsection{Optional/defaulted arguments}{}
  \item{contrast}{
number of contrast to report on, or vector of contrasts for multiple plots. See dba.show(DBA, bContrast=T) to get contrast numbers.
}
\item{method}{
 method or vector of methods to plot results for:
\itemize{ 
\item DBA_EDGER
\item DBA_DESEQ
\item DBA_EDGER_BLOCK
}
}


\item{th}{
significance threshold; all sites with FDR (or p-values, see bUsePval) less than or equal to this value will be colored red in the plot
}
  \item{bUsePval}{
logical indicating whether to use FDR (FALSE) or p-value (TRUE) for thresholding.
}
  \item{bNormalized}{
logical indicating whether to plot normalized data using normalization factors computed by differential analysis method (TRUE) or raw read counts (FALSE).
}
 \item{factor}{
string to be prepended to plot main title; e.g. factor name.
}
\item{bXY}{
logical indicating whether to draw MA plot (FALSE) or XY scatter plot (TRUE).
}

  \item{dotSize}{
size of points on plot (cex).
}
  \item{\dots}{
passed to plot.
}
}
%\details{
%}
%\value{
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
%\note{
%plots log2 values of read counts.
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(tamoxifen_analysis)
tamoxifen = dba.contrast(tamoxifen,categories=DBA_CONDITION)
tamoxifen = dba.analyze(tamoxifen)

# default MA plot
dba.plotMA(tamoxifen)

#XY plots (with raw and normalized data)
par(mfrow=c(1,2))
dba.plotMA(tamoxifen,bXY=TRUE,bNormalized=FALSE)
dba.plotMA(tamoxifen,bXY=TRUE,bNormalized=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
