plotPDFs=F
padjVal=0.05
lfcVal=0
aligner="BWA"
geneset="rptOnlyNoX" #rpt, rptOnly, rptFam, rptFamNorm. "NoX" indicates excluding X chromosome genes
#rpt is for analying protein coding genes and repeats together
#rptOnly filters out all the protein coding genes before doing statistics in DESeq2 (increases power)
#rptFam aggregates counts for whole repeat families before doing stats in DESeq2
#rptFamNorm uses aggregate family counts normalised by the size of the family
multimappers="random" # none or random

fileNamePrefix=paste0("p",padjVal,"_lfc",lfcVal,"/",aligner,"_",multimappers, "_")
filterPrefix=paste0("p",padjVal,"_lfc",lfcVal,"_",aligner,"_",
                    geneset,"_",multimappers,"/",aligner,"_",
                    geneset,"_",multimappers,"_")
if(grepl("rptFam",geneset)){
  dataset=paste0("_",multimappers,"_",geneset)
} else {
  dataset=paste0("_",multimappers)
}


outPath="."
minReads=1
genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)
dfamVer="Dfam_3.3"

filterData=T
