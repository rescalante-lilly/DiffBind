## NOTE: THIS SCRIPT NO LONGER WORKS DUE TO THE SRA CHANGING THE DATA FORMAT ON THE ARCHIVE

cat("NOTE: THIS SCRIPT NO LONGER WORKS DUE TO THE SRA CHANGING THE DATA FORMAT ON THE ARCHIVE\n")

if(FALSE){
   
library(GEOquery)

file.copy(file.path(system.file("extra", package="DiffBind"),"tamoxifen_GEO.csv"),getwd())

dir.create("GEO_DATA")

getGEOSuppFiles("GSM798430",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798431",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798443",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798440",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798423",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798424",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798425",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798428",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798429",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798442",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798432",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798433",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798444",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798426",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798427",makeDirectory=FALSE,baseDir="GEO_DATA")
getGEOSuppFiles("GSM798441",makeDirectory=FALSE,baseDir="GEO_DATA")

tamoxifen = dba(sampleSheet='tamoxifen_GEO.csv')

}

              
              
              