library(DiffBind)

pdf("DBA_vignette_plots.pdf")
savewd = getwd()
setwd(sprintf("%s/extra",system.file(package='DiffBind')))

### SECTION 3 ###

tamoxifen = dba(sampleSheet="tamoxifen.csv")
tamoxifen

#tamoxifen = dba.count(tamoxifen, minOverlap=3) 
data("tamoxifen_counts")
tamoxifen = dba.contrast(tamoxifen,categories=DBA_CONDITION)
tamoxifen = dba.analyze(tamoxifen)
tamoxifen

tamoxifen.DB = dba.report(tamoxifen)

### SECTION 4 ###

#4.1         
dba.plotMA(tamoxifen)

#4.2
dba.plotPCA(tamoxifen)
dba.plotPCA(tamoxifen, contrast=1)

#4.3
sum(tamoxifen.DB$Fold<0)
sum(tamoxifen.DB$Fold>0)
pvals = dba.plotBox(tamoxifen)
pvals

#4.4
dba.plotHeatmap(tamoxifen,contrast=1)
dba.plotHeatmap(tamoxifen,contrast=1, correlations=FALSE)


## SECTION 5 ##

data(tamoxifen_peaks)

#5.1
olap.rate = dba.overlap(tamoxifen,mode=DBA_OLAP_RATE)
olap.rate
plot(olap.rate,type='b') 

#5.2
names(tamoxifen$masks)
dba.overlap(tamoxifen,tamoxifen$masks$MCF7,mode=DBA_OLAP_RATE)
dba.plotVenn(tamoxifen,tamoxifen$masks$MCF7)

tamoxifen = dba.peakset(tamoxifen,tamoxifen$masks$MCF7,sampID="MCF7+") 
tamoxifen = dba(tamoxifen,!(!tamoxifen$masks$Consensus&tamoxifen$masks$MCF7))
tamoxifen

#5.3
data(tamoxifen_peaks)

dba.overlap(tamoxifen,tamoxifen$masks$Resistant,mode=DBA_OLAP_RATE)
dba.overlap(tamoxifen,tamoxifen$masks$Responsive,mode=DBA_OLAP_RATE)

tamoxifen = dba.peakset(tamoxifen,tamoxifen$masks$Resistant,sampID="Resistant",minOverlap=2)
tamoxifen = dba.peakset(tamoxifen,tamoxifen$mask$Responsive,sampID="Responsive",minOverlap=3)
dba.plotVenn(tamoxifen,tamoxifen$masks$Consensus)

tamoxifen.OL = dba.overlap(tamoxifen, tamoxifen$masks$Consensus)
tamoxifen.OL$onlyA
tamoxifen.OL$onlyB

#5.4
tamoxifen = dba.peakset(tamoxifen,tamoxifen$masks$Consensus,minOverlap=1,sampID="Consensus_OL")
tamoxifen = dba.peakset(tamoxifen,!tamoxifen$masks$Consensus,minOverlap=3,sampID="Consensus_3") 
dba.plotVenn(tamoxifen,14:15)

data(tamoxifen_counts)
tamoxifen = dba.contrast(tamoxifen,categories=DBA_CONDITION)
tamoxifen = dba.analyze(tamoxifen)

tamoxifen.rep = dba.report(tamoxifen,bCalled=T,th=1)
onlyResistant  = tamoxifen.rep$Called1>=2 & tamoxifen.rep$Called2<3
sum(onlyResistant)
onlyResponsive = tamoxifen.rep $Called2>=3 & tamoxifen.rep$Called1<2
sum(onlyResponsive)
bothGroups     = tamoxifen.rep$Called1>=2 & tamoxifen.rep$Called2>=3
sum(bothGroups)

tamoxifen.DB = dba.report(tamoxifen,bCalled=T,th=.1)
onlyResistant  = tamoxifen.DB$Called1>=2 & tamoxifen.SDB$Called2<3
sum(onlyResistant)
onlyResponsive = tamoxifen.SDB$Called2>=3 & tamoxifen.SDB$Called1<2
sum(onlyResponsive)
bothGroups = tamoxifen.SDB$Called1>=2 & tamoxifen.SDB$Called2>=3
sum(bothGroups)

## RESTORE WORKING DIRECTORY ###
setwd(savewd)
dev.off()
