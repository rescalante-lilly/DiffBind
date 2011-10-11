\name{dba.plotVenn}
\alias{dba.plotVenn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw 2-way or 3-way Venn diagrams of overlaps
}
\description{
Draws 2-way or 3-way Venn diagrams of overlaps
}
\usage{
dba.plotVenn(DBA, mask, overlaps, label1, label2, label3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}

\item{DBA}{
DBA object; if present, only the mask parameter will apply.
}

\item{mask}{
mask or vector of peakset numbers indicating which peaksets to include in Venn diagram. Only 2 or 3 peaksets should be included.  See dba.mask.
Only one of mask or overlaps is used.
}

\item{overlaps}{
overlap record, as computed by dba.overlap(Report=DBA_OLAP_PEAKS). 
Only one of mask or overlaps is used.
}
%\subsection{Optional/defaulted arguments}{}

  \item{label1}{
label for first peakset in diagram
}
  \item{label2}{
label for second peakset in diagram
}
  \item{label3}{
label for third peakset in diagram
}
  \item{\dots}{
arguments passed on to vennDiagram\{limma\}
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%NULL
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

par(mfrow=c(2,2))
# 2-way Venn
dba.plotVenn(tamoxifen,6:7)
dba.plotVenn(tamoxifen,tamoxifen$masks$ZR75)

# 3-way Venn (done two different ways)
dba.plotVenn(tamoxifen,tamoxifen$masks$MCF7)
olaps = dba.overlap(tamoxifen,tamoxifen$masks$MCF7)
dba.plotVenn(tamoxifen,overlaps=olaps,label1="Rep 1",label2="Rep 2",label3="Rep 3",main="MCF7 Replicates")

#Venn of overlaps
Responsive=dba(tamoxifen,tamoxifen$masks$Responsive)
Responsive
Responsive = dba.peakset(Responsive,1:3,sampID="MCF7")
Responsive = dba.peakset(Responsive,4:5,sampID="T47D")
Responsive = dba.peakset(Responsive,6:7,sampID="ZR75")
dba.plotVenn(Responsive,Responsive$masks$Consensus)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
