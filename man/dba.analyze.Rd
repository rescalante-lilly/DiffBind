\name{dba.analyze}
\alias{dba.analyze}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform differential binding affinity analysis
}
\description{
Performs differential binding affinity analysis
}
\usage{
dba.analyze(DBA, method=DBA$config$AnalysisMethod, 
            bSubControl=TRUE, bFullLibrarySize=FALSE, bTagwise=TRUE,
            bCorPlot=TRUE,  bParallel=DBA$config$RunParallel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}

  \item{DBA}{
DBA object. If no contrasts are specified (DBA$contrast is NULL), all deault contrasts will be added via a call to dba.contrast(DBA).
}
%\subsection{Optional/defaulted arguments}{}
  \item{method}{
method, or vector of methods, by which to analyze differential binding affinity. Supported methods:
\itemize{ 
\item DBA_EDGER
\item DBA_DESEQ
}
}

\item{bSubControl}{
logical indicating whether Control read counts are subtracted for each site in each sample before performing analysis.
}
\item{bFullLibrarySize}{
logical indicating if the full library size (total number of reads in BAM/SAM/BED file) for each sample is used for scaling normalization. If FALSE, the total number of reads present in the peaks for each sample is used (generally preferable).
}
\item{bTagwise}{
logical indicating if dispersion should be calculated on a tagwise (or per-condition) basis. If there are only a very few members of each group in a contrast (e.g. no replicates), this should be set to FALSE.
}
\item{bCorPlot}{
logical indicating whether to plot a correlation heatmap for the analyzed data (first contrast only). 
If no sites are significantly differentially bound using the default threholds, no heatmap will be plotted.
}
\item{bParallel}{
logical indicating that the analyses is to be done in parallel using multicore (one process for each contrast for each method, plus an additional process per method).
}
}
\details{
See the DBA User Guide for more details on how the edgeR and DESeq analyses are carried out.  
}
\value{
DBA object with results of analysis added to DBA$contrasts. 
}

\author{
Rory Stark
}
\note{
If the "edgeR" method is specified, and there is a blocking factor for the contrast(s) specified using a previous call to dba.contrast, a multi-factor analysis will automatically be carried out in addition to a single factor analysis.
}



\examples{
data(tamoxifen_counts)

tamoxifen = dba.analyze(tamoxifen)

}


