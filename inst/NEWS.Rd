\name{DiffBindNews}
\title{DiffBind News}
\encoding{UTF-8}

\section{version 1.9.16}{
 \itemize{
  \item{Vignette updated to use BiocStyles and dynamically generated Figures.}
  \item{Plotting colors updated to CRUK colors.}  
  \item{Data objects updaed to use hg19 alignments/peaks and bam instead of bed.}
  }
}


\section{version 1.9.15}{
 \itemize{
  \item{Feature: added legends to dba.plotPCA}
  \item{Feature: added label parameter to dba.plotPCA to lable point in plot}  
  \item{Feature: added mapQCth and fragmentSize to config, and made config a list}
  }
}

\section{version 1.9.14}{
 \itemize{
  \item{Feature: added mapQCth filter to dba.count.}
  }
}


\section{version 1.9.9}{
 \itemize{
  \item{Feature: added support to pass in "fragments" when using summarizeOverlaps to count.}
  \item{Fix: bug in dba when samplesheet passed in as data frame.}   
  \item{Change: now using FRIP (fractions of reads in peaks) instead of SN (signal-to-noise).}
  }
}

\section{version 1.9.7}{
 \itemize{
   \item{New set of features relating to computing peak summits}
  \itemize{
     \item{Feature: dba.count computes summit heights and locations for each peak.}
     \item{Feature: new scores: DBA_SCORE_SUMMIT and DBA_SCORE_SUMMIT_ADJ}
     \item{Feature: dba.count re-centers consensus peakset on peak summits, with defined width}
   }
   \item{Interface changes to dba.count}
   \itemize{
     \item{removed parameter bCalledMasks (now always TRUE)}
     \item{changed insertLength to fragmentSize}
     \item{fragmentSize default is now 125 bp}
     \item{fragmentSize can also be a vector with separate values for each library}
   }
   \item{Feature: new TMM scores reporting values in CPM}
   \item{Feature: can config custom names for ID, Tissue, Factor, Condition, Treatment, Caller}
   \item{Fix: disabled tamoxifen_GEO.R as SRA changed data formats and script no longer works}
   \item{Fix: shortened length of some very long lines in man pages and examples}
 }
}

\section{version 1.9.5}{
 \itemize{
  \item{Fix: bug when generating GRanges overlaps and one peakset is empty.}
  \item{Fix: need to be dependant on limma because of edgeR call to nonEstimable}
  \item{Fix: summarizeOverlaps moved to GenomicAlignments package; add dependancy/import}
 }
}

\section{version 1.9.3}{
 \itemize{
  \item{Fix: Missing header file introduced in version 1.9.1.}
 }
}

\section{version 1.9.2}{
\itemize{
  \item{Finalize deprecation of bLowMem in dba.count.}
  \item{add new constants: DBA_ALL_ATTRIBUTES, DBA_ALL_METHODS_BLOCK}
  \item{add config options th and bUsePval to set default thresholds on a per-DBA object basis.}
  \item{redo report-based DBA object interface}
  \itemize{
    \item{bDown is now bLoss}
    \item{bUp is now bGain}
    \item{default values changed}
    \item{in resulting object, all is now All, down us now Loss, and up is now Gain}
  }
  \item{pv.plotVenn changes/improvements}
  \itemize{
     \item{labelAttributes parameter added for generaiton of default labels}
     \item{more checks for consistency of masks}
     \item{default main and sub titles added}
     \item{added plotting of DB sites by specifying contrasts and other parameters similar to report-dbased DBA objects in dba.report}
     \item{added bReturnPeaksets to enable all peaksets in the Venn diagram ot be returned}
     \item{added DataType to control clas of returned peaskets, including DBA_DATA_DBAOBJECT to retrieve a DBa object made up of DB sites}
  }
  \item{prevent dba.plotHeatmap from plotting column vector for attributes where every sample has a different value}
  \item{warning and errors; fix non-expoerted constants}
 }
}

\section{version 1.9.1}{
 \itemize{
  \item{Fix: longstanding bug in dba.count that could case reads mapped to chrN_random to be counted aginst intervals on chrN.}
 }
}

\section{version 1.8.0}{

Changes since last release:

\itemize{

     \item{Add support for DESeq2:}{
     \itemize{
      \item{New: Add DBA_DESEQ2, DBA_ALL_METHODS and DBA_ALL_BLOCK method constants}     
      \item{Change: dba.analyze can analyze using DESeq2}
      \item{Change: all reporting and plotting functions support DESeq2 results}
      \item{Change: vignette includes comparison of edgeR, DESeq, and DESeq2}
     }
   }
    
     \item{Changes to counting using dba.count:}{
       \itemize{
        \item{Change: optimize built-in counting code to use much less memory and run faster}
        \item{Change: deprecate bLowMem, replaced with bUseSummarizeOverlaps}
        \item{New: add readFormat parameter to specify read file type (instead of using file suffix)}
       }
     }
     
   \item{New: generation of result-based DBA object using dba.report (makes it easier to work with differentially bound peaksets)}{}
     
    \item{Changes to defaults:}{
      \itemize{
        \item{Change: default score is now DBA_SCORE_TMM_MINUS_FULL instead of DBA_SCORE_TMM_MINUS_EFFECTIVE in dba.count}
       \item{Change: default value for bFullLibrarySize is now TRUE in dba.analyze}
        \item{New: add bCorPlot option to DBA$config to turn off CorPlot by default}
      }
    }
    
   \item{Various bugfixes, improved warnings, updated documentation}{}
  }
}

\section{Version 1.6.0}{
   \itemize{
      \item{New: Low memory counting of bam files using Rsamtools and summarizeOverlaps (bLowMem in dba.count)}
      \item{New: Ability to read in externally derived counts (e.g. from htSeq) (dba.count)}
      \item{Improved: Features to deal with filtering intervals based on read scores (dba.count)}
         \itemize{
            \item{Change parameter name: maxFilter -> filter}
            \item{Allow maxFilter to be a numerical vector to retrieve filtering rate}
            \item{Add parameter: filterFun to control filtering method}
        }
      \item{New: Support for SummarizedExperiment objects (dba and dba.report)} 
          \itemize{
            \item{Add bSummarizedExperiment option to dba() to convert DBA object}
            \item{Add DataType = DBA_DATA_SUMMARIZED_REPORT option to dba.report() to return SummarizedExperiment}
         }
      \item{Documentation: Add section to vignette showing how to obtain full tamoxifen resistance dataset}  
         \itemize{
            \item{Add section to vignette showing how to obtains full tamoxifen dataset}
            \item{Add script (tamoxifen_GEO.R) and sample sheet (tamoxifen_GEO.csv) to extras for full tamoxifen dataset}
            \item{Add examples to man page for dba.count to show filtering}
            \item{Add examples to man pages for dba and dba.report to show retrieval of SummarizedExperiment objects}
            \item{Update and cleanup vignette and man pages}
         }
      \item{Various bugfixes and improved warnings}      
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


\section{version 1.2.0:(2012-03-30)}{
 \itemize{
   \item{GRanges is default class for peaksets and reports instead of RangedData, controlled by DataType parameter.}

   \item{Both analysis methods (edgeR and DESeq) use generalized linear models (GLMs) for two-group contrasts by default.}

   \item{Blocking factors (for two-factor analysis) can be specified flexibly such that arbitrary blocking factors can be used. }

   \item{Section added to vignette showing an an analysis using a blocking factor.}

   \item{Added new metadata type, DBA_TREATMENT.}

   \item{New DBA_SCORE_ options for specifying scoring method, including TMM normalized counts, and ability to change scoring 
             method on the fly in dba.plotHeatmap and dba.plotPCA when plotting global binding matrix.}

   \item{bRemoveDuplicates parameter in dba.count allows duplicate reads to be discarded when computing counts}

   \item{More efficient use of memory when analyzing (controlled by bReduceObjects parameter in dba.analyze).}

   \item{various bugs fixed, man pages updated, and warning messages added.}
  }
}

\section{version 1.0.0}{
   \itemize{
     \item{Initial release}
  }
}
