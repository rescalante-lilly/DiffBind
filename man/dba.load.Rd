\name{dba.load}
\alias{dba.load}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
load DBA object
}
\description{
Reads in saved DBA object
}
\usage{
dba.load(file='DBA', dir='.', pre='dba_', ext='RData')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%\subsection{Optional/defaulted arguments}{}
  \item{file}{
main filename
}
  \item{dir}{
directory in which to save model
}
  \item{pre}{
string to pre-pend to filename
}
  \item{ext}{
file extension to use
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
loaded DBA object
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
\code{\link{dba.save}}
}

\examples{
data(tamoxifen_peaks)
dba.save(tamoxifen,'tamoxifenPeaks')
tamoxifen = dba.load('tamoxifenPeaks')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
