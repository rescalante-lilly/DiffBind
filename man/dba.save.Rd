\name{dba.save}
\alias{dba.save}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
save DBA object
}
\description{
Writes out DBA object
}
\usage{
dba.save(DBA, file='DBA', dir='.', pre='dba_', ext='RData', bMinimize=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object
}
%\subsection{Optional/defaulted arguments}{}
  \item{file}{
main filename
}
  \item{dir}{
directory to save model in
}
  \item{pre}{
string to pre-pend to filename
}
 \item{ext}{
extensions to use
}
\item{bMinimize}{
logical indicating saved DBA object should be compressed as much as possible.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
string containing full path and filename.
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
savefile = dba.save(tamoxifen,'tamoxifenPeaks')
savefile
tamoxifen = dba.load('tamoxifenPeaks')
unlink(savefile)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
