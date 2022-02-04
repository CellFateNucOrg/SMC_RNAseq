library(rtracklayer)
library(GenomicRanges)


df<-data.frame(sample=c("382_B","775_B","784_B","828_C","844_C"),
control=c(rep("366_B",3),rep("366_C",2)))


makeRatio<-function(sampleName, controlName,outName){
	ref1<-import.bw(controlName)
	cov1<-coverage(ref1,weight="score")
	treat1<-import.bw(sampleName)
	cov2<-coverage(treat1,weight="score")
	ratio<-log2((cov2+mean(unlist(cov2)))/(cov1+mean(unlist(cov1))))
	export(GRanges(ratio),outName,format="bw")
}


for(i in 1:nrow(df)){
	controlName<-paste0("./tracks/avr_",df$control[i],"_F_UniqueMultiple_RPM.bw")
	sampleName<-paste0("./tracks/avr_",df$sample[i],"_F_UniqueMultiple_RPM.bw")
	outName<-paste0("./tracks/lfc_",df$sample[i],"-",df$control[i],"_F_STAR.bw")
	makeRatio(sampleName,controlName,outName)
	controlName<-paste0("./tracks/avr_",df$control[i],"_R_UniqueMultiple_RPM.bw")
	sampleName<-paste0("./tracks/avr_",df$sample[i],"_R_UniqueMultiple_RPM.bw")
	outName<-paste0("./tracks/lfc_",df$sample[i],"-",df$control[i],"_R_STAR.bw")
	makeRatio(sampleName,controlName,outName)
}

