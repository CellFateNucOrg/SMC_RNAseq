library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)

df<-data.frame(sample=c("382_B","775_B","784_B"),
control=c(rep("366_B",3)))

ws235toCe11<-function(gr){
	seqlevels(gr)<-c("I","II","III","IV","V","X","MtDNA")
	seqlevels(gr)<-seqlevels(BSgenome.Celegans.UCSC.ce11::Celegans)
	seqinfo(gr)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
	return(gr)
}

makeRatio<-function(sampleName, controlName,outName){
	ref1<-import.bw(controlName)
	cov1<-coverage(ref1,weight="score")
	treat1<-import.bw(sampleName)
	cov2<-coverage(treat1,weight="score")
	ratio<-log2((cov2+mean(unlist(cov2)))/(cov1+mean(unlist(cov1))))
	export(GRanges(ratio),outName,format="bw")
}


for(i in 1:nrow(df)){
	controlName<-paste0("./tracks/avr_",df$control[i],"_F_ribo0_UniqueMultiple_RPM.bw")
	sampleName<-paste0("./tracks/avr_",df$sample[i],"_F_ribo0_UniqueMultiple_RPM.bw")
	outName<-paste0("./tracks/lfc_",df$sample[i],"-",df$control[i],"_F_ribo0_STAR.bw")
	makeRatio(sampleName,controlName,outName)
	controlName<-paste0("./tracks/avr_",df$control[i],"_R_ribo0_UniqueMultiple_RPM.bw")
	sampleName<-paste0("./tracks/avr_",df$sample[i],"_R_ribo0_UniqueMultiple_RPM.bw")
	outName<-paste0("./tracks/lfc_",df$sample[i],"-",df$control[i],"_R_ribo0_STAR.bw")
	makeRatio(sampleName,controlName,outName)
}


