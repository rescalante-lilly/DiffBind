################################################
## cpp_wrapper.R -- wrap C/C++ functions in R ##
## first created: 2014 03 12                  ##
## Gord Brown                                 ##
## CR-UK Cambridge Institute                  ##
################################################

################################################
## Entry Points to C/C++ Routines             ##
################################################

## cpp_count_reads -- count reads on intervals
## cpp_mergeOne -- merge overlapping intervals into extended intervals
## cpp_mergeTwo -- merge two sets of intervals, calculating the union of them

cpp_count_reads <- function(bamfile,insertLength,fileType,bufferSize,
                            intervals,bWithoutDupes,summits,minMappingQual=0) {
  icount <- length(intervals[[1]])
  counts <- vector(mode="integer",length=icount)
  if (!missing(summits)) {
    summits.vec <- vector(mode="integer",length=icount)
    heights.vec <- vector(mode="integer",length=icount)
    bSummits = TRUE
  } else {
    summits.vec <- vector()
    heights.vec <- vector()
    bSummits = FALSE
  }
  libsize <- .Call("croi_count_reads",bamfile,
                   as.integer(insertLength),
                   as.integer(fileType),
                   as.integer(bufferSize),
                   as.integer(minMappingQual),
                   as.character(intervals[[1]]),
                   as.integer(intervals[[2]]),
                   as.integer(intervals[[3]]),
                   as.integer(icount),
                   as.logical(bWithoutDupes),
                   as.logical(bSummits),
                   counts,
                   summits.vec,
                   heights.vec)
  counts[counts==0]=1

  widths = intervals[,3] - intervals[,2]
  rpkm = (counts/(widths/1000))/(libsize/1E6)

  result <- list(counts=counts,rpkm=rpkm,libsize=libsize)
  if (bSummits==T) {
    result$summits <- summits.vec;
    result$heights <- heights.vec;
  }
  return(result)
}


cpp_mergeOne <- function(peaks,bKeepAll,minOverlap) {
  return(.Call("mo_mergeOne",peaks,bKeepAll,minOverlap))
}


cpp_mergeTwo <- function(peaksA,peaksB,bKeepAll,minOverlap,zero) {
  return(.Call("mo_mergeTwo",peaksA,peaksB,bKeepAll,minOverlap,zero))
}
