\name{DiffBindNews}
\title{DiffBind News}
\encoding{UTF-8}

\section{Version 1.5.10}{
   \itemize{
      \item{Add alias for DBA in man page for print.DBA}
      \item{Add tamoxifen_GEO.R and tamoxifen_GEO.csv to extras; update vignette to explain how to use them to obtain all of the experimental data for the vignette from GEO}
   }
}

\section{Version 1.5.9}{
   \itemize{
      \item{Add bSummarizedExperiment option to dba() to convert DBA object}
      \item{Add DataType = DBA_DATA_SUMMARIZED_REPORT option to dba.report() to return SummarizedExperiment}
   }
}

\section{Version 1.5.5}{
   \itemize{
      \item{Add bLowMem option to dba.count for low memory counting using summarizeOverlaps}
   }
}

\section{Version 1.5.4}{
   \itemize{
      \item{Add ability to add externally derived counts via dba.peakset and a sample sheet}
   }
}

\section{Version 1.5.2}{
   \itemize{
      \item{add features to dba.count to deal with filtering}
      \itemize{
         \item{Change parameter name: maxFilter -> filter}
         \item{Allow maxFilter to be a numerical vector to retrieve filtering rate}
         \item{Add parameter: filterFun to control filtering method}
         \item{Fix bug so SN gets recomputed after filtering}
         \item{Update documentation with example}
      }
   }
}

\section{Version 1.5.1}{
   \itemize{
      \item{Roll up bugfixes for 1.4.0}
   }
}

\section{Version 1.4.0}{
   \itemize{
      \item{Plotting}
      \itemize{
          \item{dba.plotMA}
          \itemize{
             \item{Smooth plots now default}
             \item{Added fold parameter in addition to th (threshold)}
          }
          \item{dba.plotHeatmap}
          \itemize{
             \item{Side colorbars added}
                \item{Add support for specifying sample mask to include any subset of samples in a contrast plot, including samples that were not in the original contrast}
           }
           \item{dba.plotVenn}
           \itemize{
                \item{Changed plotter from limma to T. Girke's overLapper}
                \item{Added support for 4-way Venns (also in dba/overlap)}
	  }
	  \item{dba.plotPCA}
	  \itemize{
               \item{Add support for specifying sample mask to include any subset of samples in a contrast plot, including samples that were not in the original contrast}  
            }
      }
      \item{Peaksets (dba and dba.peakset)}
      \itemize{
         \item{Peakset formats}
         \itemize{
            \item{narrowPeaks format supported}
            \item{Can override file format, score column, and score orientation defaults for supported peak callers}
         }
         \item{Consensus peaksets}
         \itemize{
            \item{Added ability to generate sets of consensus peaksets based on metadata attributes: for example create consensus peaksets for each tissue type and/or condition, or for all unique samples by taking the consensus of their replicate peaksets}
         }
      }
      \item{Read counting (dba.count)}
      \itemize{
         \item{Compute Signal-to-Noise ratio when counting}
         \item{Added bScaleControl to down-scale control reads by default}
         \item{Add option to specify a mask in peak parameter to limit which peaksets are used to for a consensus by overlap. Works with new consensus peakset options in dba.peakset}
         \item{Remove references to support for SAM files}
      }
      \item{Analysis (dba.analyze)}
      \itemize{
         \item{edgeR: updated calls to math change sin edgeR; updated vignette and references}
         \item{DESeq: updated to work with current DESeq; use pooled-CR dispersion estimation method for blocking analysis; update vignette}
      }
      \item{Various bug fixes; more informative warnings; update documentation including vignette, new examples and cross-referencing in man pages}
   }
 }

\section{version 1.2.3:(2012-09-01)}{
 \itemize{
   \item{more informative warnings and minor bug fixes.}
 }
}

\section{version 1.2.0:(2012-03-30)}{
 \itemize{
   \item{GRanges is default class for peaksets and reports instead of RangedData, controlled by DataType parameter.}

   \item{Both analysis methods (edgeR and DESeq) use generalized linear models (GLMs) for two-group contrasts by default.}

   \item{Blocking factors (for two-factor analysis) can be specified flexibly such that arbitrary blocking factors can be used. }

   \item{Section added to vignette showing an ananalysis using a blocking factor.}

   \item{Added new metadata type, DBA_TREATMENT.}

   \item{New DBA_SCORE_ options for specifying scoring method, including TMM normalized counts, and ability to change scoring 
             method on the fly in dba.plotHeatmap and dba.plotPCA when plotting global binding matrix.}

   \item{bRemoveDuplicates parameter in dba.count allows duplicate reads to be discarded when computing counts}

   \item{More efficient use of memory when analyzing (controlled by bReduceObjects parameter in dba.analyze).}

   \item{various bugs fixed, man pages updated, and warning messages added.}
  }
}
