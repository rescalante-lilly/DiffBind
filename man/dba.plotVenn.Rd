\name{dba.plotVenn}
\alias{dba.plotVenn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw 2-way, 3-way, or 4-way Venn diagrams of overlaps
}
\description{
Draws 2-way, 3-way, or 4-way Venn diagrams of overlaps
}
\usage{
dba.plotVenn(DBA, mask, overlaps, label1, label2, label3, label4, main="", sub="")
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
  \item{label4}{
label for fourth peakset in diagram
}
  \item{main}{
main title for plot
}
  \item{sub}{
subtitle for plot
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
Rory Stark and Thomas Girke
}
\note{
Venn plotting code written by Thomas Girke as part of overLapper code:

http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#R_graphics_venn
}

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
dba.plotVenn(tamoxifen,tamoxifen$masks$MCF7&tamoxifen$masks$Responsive)
olaps = dba.overlap(tamoxifen,tamoxifen$masks$MCF7&tamoxifen$masks$Responsive)
dba.plotVenn(tamoxifen,overlaps=olaps,
             label1="Rep 1",label2="Rep 2",label3="Rep 3",main="MCF7 (Responsive) Replicates")

#Venn of overlaps
Responsive=dba(tamoxifen,tamoxifen$masks$Responsive)
Responsive
Responsive = dba.peakset(Responsive,1:3,sampID="MCF7")
Responsive = dba.peakset(Responsive,4:5,sampID="T47D")
Responsive = dba.peakset(Responsive,6:7,sampID="ZR75")
par(mfrow=c(1,1))
dba.plotVenn(Responsive,Responsive$masks$Consensus)

#4-way overlap
data(tamoxifen_peaks)
tamoxifen = dba.peakset(tamoxifen, consensus=DBA_TISSUE)
par(mfrow=c(1,1))
dba.plotVenn(tamoxifen,tamoxifen$masks$Consensus,main="Tissue consensus overlaps")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
