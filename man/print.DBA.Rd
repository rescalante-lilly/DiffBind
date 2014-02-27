\name{DBA object methods}
\alias{DBA}
\alias{print.DBA}
\alias{summary.DBA}
\alias{plot.DBA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard S3 methods for DBA object
}
\description{
Standard S3 methods for DBA object.
}
\usage{
\method{print}{DBA}(x, \dots)
\method{summary}{DBA}(object, \dots)
\method{plot}{DBA}(x, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
DBA object
}
\item{object}{
DBA object
}

\item{\dots}{
Arguments passed on to parent methods
}
}
\details{
S3 methods for DBA object from the \code{\link{DiffBind}} package.

DBA objects are usually constructed using the \code{\link{dba}} function.
}
%\value{
%}
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
tamoxifen
data(tamoxifen_counts)
tamoxifen
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
