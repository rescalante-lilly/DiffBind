\name{DiffBindNews}
\title{DiffBind News}
\encoding{UTF-8}

\section{Version 1.14.4}{
\itemize{
   \item{Fix bug when score=DBA_SCORE_TMM_READS_FULL_CPM}
 }
}


\section{Version 1.14.2}{
\itemize{
   \item{Fix bug when not all chromosome names in first peakset}
 }
}

\section{Version 1.14.1}{
\itemize{
   \item{Fix bug when summit heights are zero}
   \item{Update version number in object}
   \item{dba.plotHeatmap returns binding sites in row order}
 }
}

\section{Version 1.14.0}{
\itemize{
   \item{Add support for reading Excel-format sample sheets (.xls, .xlsx extensions}
   \item{Update DESeq2 reference in vignette; fix vignette samplesheet}
   \item{use vennPlot from systemPiper}
   \item{Fix Makevars to avoid gnu-specific extensions}
   \item{Replace 'require' with 'requireNamespace' to eliminate NOTEs regarding misuse of 'require'}
   \item{Remove non-ASCII characters from a couple of comments}
   \item{Change Gord's email address}
   \item{New: color vector lists for dba.plotHeatmap and colos for dba.plotPCA labels}
   \item{Fix: bug causing two plots when changing score in dba.plotHeatmap and dba.plotPCA}
 }
}

\section{Version 1.12.0}{
\itemize{
   \item{Mostly bug fixes!}
 }
}

\section{Version 1.10.0}{

Changes since last release:

\itemize{
   \item{Counting}
      \itemize{
         \item{New: option to compute summits}
         \item{New: option to center peaks with fixed width around summits}
         \item{New: scores for summits (height, position) and CPM for TMM values}
         \item{New: filter reads by mapping quality (mapQCth)}
         \item{New: support for PE bam data using summarizeOverlaps}
         \item{Remove: bCalledMask (now always TRUE)}
         \item{Change: insertLength to fragmentSize}
         \item{Add: fragmentSize can be a vector with a size for each sample}
         \item{Change: fragmentSize default is 125 bp}
      }
   \item{Plotting}
      \itemize{
        \item{Change: colors based on CRUK color scheme}
        \item{PCA plots}
           \itemize{
              \item{New: legend}
              \item{New: label parameter for adding text labels of points in 2D plot}
           }
         \item{Venn diagrams}
         \itemize{
            \item{New: plot overlaps of differentially bound sites by specifying contrasts, thresholds etc.}
            \item{New: able to return overlapping peaksets as GRanges directly}
            \item{New: able to generate new DBA object consisting of overlapping peaks}
            \item{New: labelAttributes for controlling default labels}
            \item{New: default main and sub titles}
         }
        \item{Heatmaps}
         \itemize{
            \item{Fix: don’t plot column vector for attributes where every sample has a different value}
         }
      }

   \item{General}
   \itemize{
     \item{New: add attribute value: DBA_ALL_ATTRIBUTES}
     \item{Change: SN  (signal/noise) to FRIP (fraction of reads in peaks)}
     \item{Change: “Down” to “Loss” and Up” to “Gain”}
     \item{Vignette}
       \itemize{
          \item{Change: vignette uses BiocStyles and dynamically generated figures}
           \item{Change: example data based on hg19 instead of hg18}
           \item{Change: example reads from bam files instead of bed files}
            \item{New: section on using DiffBind and ChIPQC together}
     }
     \item{New configuration defaults options (DBA$config):}
        \itemize{
          \item{Metadata name strings: ID, Tissue, Factor, Condition, Treatment, Caller}
          \item{th: significance threshold}
          \item{bUsePval}
          \item{fragmentSize}
          \item{mapQCth: filter reads by mapping quality}
          \item{fragments (for summarizeOverlaps)}
      }
   }
   \item{Bugs/Issues}
   \itemize{
      \item{Fix: bRemoveDuplicates had some unpredictable behaviour}
      \item{Fix: chrN_random were being counted against chrN}
      \item{Disable: tamoxifen_GEO.R doesn’t work after SRA changed format of archived data}
   }
}
}


\section{Version 1.8.0}{

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
