
workDir="."
#tbl<-read.csv(paste0(workDir,"/McMurchy2017_PRJNA345014.txt"))

rowIdx<-(grepl("WT",tbl$genotype.variation) |
           grepl("met-2",tbl$genotype.variation)) &
          tbl$growth_condition=="condition1"

setmet<-tbl[rowIdx,c("Run","Age","Assay.Type","AvgSpotLen","BioProject","BioSample",
               "Experiment","genotype.variation","GEO_Accession..exp.",
               "growth_condition","LibraryLayout","Sample.Name",
               "strain.background","Tissue")]



setmet

SRRs<-data.frame(SRRnumber=setmet$Run,
                 dataset=setmet$BioProject,
                 bioType=paste0(gsub("\\([[:alnum:]]*\\)","",
                                     setmet$genotype.variation),
                              "_ribo0"),
                 replicate=paste0(1:2))

write.table(SRRs,paste0(workDir,"/SRR_McMurchy2017_setmet.tsv"),quote=F,
            sep="\t",row.names=F)



mcmurchyAll<-tbl[,c("Run","Age","Assay.Type","AvgSpotLen","BioProject","BioSample",
                     "Experiment","genotype.variation","GEO_Accession..exp.",
                     "growth_condition","LibraryLayout","Sample.Name",
                     "strain.background","Tissue")]

SRRs<-data.frame(SRRnumber=mcmurchyAll$Run,
                 dataset=mcmurchyAll$BioProject,
                 bioType=paste0(gsub("\\([[:alnum:]]*\\)","",
                                     mcmurchyAll$genotype.variation),
                                "_ribo0"),
                 replicate=paste0("r",1:2),
                 condition=gsub("condition","c",mcmurchyAll$growth_condition))


write.table(SRRs,paste0(workDir,"/SRR_McMurchy2017_all.tsv"),quote=F,
            sep="\t",row.names=F)
