\name{dba.contrast}
\alias{dba.contrast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set up contrasts for differential binding affinity analysis
}
\description{
Sets up contrasts for differential binding affinity analysis
}
\usage{
dba.contrast(DBA, group1, group2=!group1, name1="group1", name2="group2",
             minMembers=3, block,
             categories = c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION,DBA_TREATMENT))

}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object with count data
}
%\subsection{Optional/defaulted arguments}{}
  \item{group1}{
mask of samples in first group (when adding a specific contrast). See \code{\link{dba.mask}}.
}
  \item{group2}{
mask of samples in second group (when adding a specific contrast). See \code{\link{dba.mask}}.
}
  \item{name1}{
label for samples in first group (when adding a specific contrast). 
}
  \item{name2}{
label for samples in second group (when adding a specific contrast). 
}
  \item{minMembers}{
when automatically generating contrasts, minimum number of unique samples in a group. Must be at least 2, as replicates are strongly advised. If you wish to do an analysis with no replicates, you can set the group1 and group2 parameters explicitly.
}
  \item{categories}{
when automatically generating contrasts, attribute or vector of attributes to base contrasts on:
\itemize{ 
\item DBA_ID
\item DBA_TISSUE
\item DBA_FACTOR
\item DBA_CONDITION
\item DBA_TREATMENT
\item DBA_REPLICATE
\item DBA_CALLER
}
}
\item{block}{
blocking attribute for multi-factor analysis.  This may be specified as either a value, a vector, or a list. 

If block is a value, the specified metadata field is used to derive the blocking factor. One of:
\itemize{ 
\item DBA_TISSUE
\item DBA_FACTOR
\item DBA_CONDITION
\item DBA_TREATMENT
\item DBA_REPLICATE
\item DBA_CALLER
}

If block is a vector, it can either be a mask (logical vector) or a vector of peakset numbers. In this case, the peaksets indicated in the blocking vector are all given the same value (true), while any peaksets not included in the vector take the alternative value (false).

If block is a list, it should be a list of vectors (either logical masks or vectors of peakset numbers), with each indicating a set of peaksets that should share the same value. Each peasket should appear at most once, and any peaksets not specified will be given an default value (other).

}
}
\details{

MODE: Set up all possible contrasts:

dba.contrast(DBA, minMembers, categories)

MODE: Set up a specific contrast:

dba.contrast(DBA, group1, group2, name1, name2, block)
}
\value{
DBA object with contrast(s) set as DBA$contrasts. Contrast list can be retrieved using dba.show(DBA, bContrasts=T).
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
\note{
Contrasts will only be set up for peaksets where DBA_CALLER == "counts".

Contrasts can be cleared by DBA$contrasts=NULL.
}

\seealso{

\code{\link{dba.analyze}}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(tamoxifen_counts)
tamoxifen = dba.contrast(tamoxifen, categories=DBA_CONDITION)
tamoxifen

# Another way to do the same thing
tamoxifen$contrasts=NULL
tamoxifen = dba.contrast(tamoxifen, tamoxifen$masks$Responsive, tamoxifen$masks$Resistant,
                                              "Responsive", "Resistant")
tamoxifen

# Add add default contrasts
tamoxifen$contrasts=NULL
tamoxifen = dba.contrast(tamoxifen)
tamoxifen

# Specify a blocking factor
tamoxifen$contrasts=NULL
tamoxifen = dba.contrast(tamoxifen, categories=DBA_CONDITION, block=DBA_TISSUE)
tamoxifen

tamoxifen$contrasts=NULL
tamoxifen = dba.contrast(tamoxifen, categories=DBA_CONDITION, block=list(c(3,4,5,8,9),c(1,2,10,11)))
tamoxifen

tamoxifen$contrasts=NULL
tamoxifen = dba.contrast(tamoxifen, categories=DBA_CONDITION, block=tamoxifen$masks$MCF7)
tamoxifen = dba.analyze(tamoxifen)
tamoxifen

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
