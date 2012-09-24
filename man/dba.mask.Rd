\name{dba.mask}
\alias{dba.mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Derive a mask to define a subset of peaksets or sites for a DBA object
}
\description{
Derives a mask to define a subset of peaksets or sites for a DBA object.
}
\usage{
dba.mask(DBA, attribute, value, combine='or', mask, merge='or', bApply=FALSE,
         peakset, minValue=-1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object
}
  \item{attribute}{
when deriving a peakset mask, attribute to base mask on:
\itemize{  	
\item DBA_ID
\item DBA_TISSUE
\item DBA_FACTOR
\item DBA_CONDITION
\item DBA_TREATMENT
\item DBA_REPLICATE
\item DBA_CONSENSUS
\item DBA_CALLER
\item DBA_CONTROL
}
}
%\subsection{Optional/defaulted arguments}{}
  \item{value}{
when deriving a peakset/sample mask, attribute value (or vector of attribute values) to match.
}

  \item{combine}{
when deriving a peakset/sample mask, if value is a vector, OR when deriving a site mask, and peaksets is a vector, this is method for combining result of each value:
\itemize{ 
\item \dQuote{or}
\item \dQuote{and}
\item \dQuote{nor}
\item \dQuote{nand}  
}
}
\item{mask}{
when deriving a peakset/sample mask, this specifies an existing mask to merge with; if missing, create new mask
}
  \item{merge}{
when deriving a peakset/sample mask, and an existing mask is supplied, this speficies the method for combining new mask with supplied mask:
\itemize{ 
\item \dQuote{or}
\item \dQuote{and}
\item \dQuote{nor}
\item \dQuote{nand}  note: if mask is missing, \dQuote{nand} results in negative of mask
}
}
\item{bApply}{
when deriving a peakset/sample mask, a logical indicating that a new DBA object with the mask applied will be returned.
}

\item{peakset}{
when deriving a peak/site mask, this specifies a peakset number, or a vector of peakset numbers. The resulting mask will indicate which of the overall sites were called as peaks in this peakset or set of peaksets. If a vector, the masks for each of the peaksets will be combined using the method specified in the combine parameter.
}

\item{minValue}{
when deriving a peak/site mask, scores greater than this value will be considered as indicating that the site corresponds to a called peakset.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~

MODE: Derive a a mask of peaksets/samples:

dba.mask(DBA, attribute, value, combine, mask, merge, bApply)

MODE: Derive a mask of peaks/sites:

dba.mask(DBA, combine, mask, merge,bApply, peakset, minValue)
}
\value{
either a logical mask, or new DBA object if bApply is TRUE.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
\note{
%%  ~~further notes~~
dba automatically generates masks for each unique value of DBA_TISSUE, DBA_FACTOR, DBA_CONDITION, DBA_TREATMENT,  DBA_CALLER, and DBA_REPLICATE. These are accessible using masks field of the DBA object (DBA$masks), and can be viewed using names(DBA$masks).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dba.show}}
}
\examples{
data(tamoxifen_peaks)

# Pre-made masks
names(tamoxifen$masks)
dba.show(tamoxifen,tamoxifen$masks$MCF7)

# New masks
mcf7Mask = dba.mask(tamoxifen,DBA_TISSUE, "MCF7")
mcf7DerivedMask = dba.mask(tamoxifen,DBA_TISSUE,"TAMR",mask=mcf7Mask)
mcf7Derived = dba(tamoxifen,mcf7DerivedMask)
mcf7Derived
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
