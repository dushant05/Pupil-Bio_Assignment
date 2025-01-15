
###################################################CNV CALLS#############################################

#source ("http://bioconductor.org/biocLite.R")
#biocLite ("CNVPanelizer")
library(CNVPanelizer)
#library(devtools)
#install_github("biostuff/CNVPanelizer")


args = commandArgs(TRUE)

#reading targets
#genomicRangesFromBed <-  BedToGenomicRanges("~/Reference_Files/Human_reference/targetRegions/350.bed",ampliconColumn = 4,split = "_")
genomicRangesFromBed <-  BedToGenomicRanges(args[1],ampliconColumn = 4,split = "_")
metadataFromGenomicRanges <-  elementMetadata(genomicRangesFromBed)
geneNames = metadataFromGenomicRanges["geneNames"][, 1]
ampliconNames = metadataFromGenomicRanges["ampliconNames"][, 1]
# reading sample & reference reads
# referenceReadCounts <-  ReadCountsFromBam("~/Reference_Files/Human_reference/healthy/TGS350/pooledNormalFemale.bam",
#                                           gr = genomicRangesFromBed,ampliconNames = ampliconNames,
#                                           sampleNames = "Normal",removeDup = FALSE)
#referenceReadCounts <-  ReadCountsFromBam(paste0(args[3],".bam"),
#                                          gr = genomicRangesFromBed,ampliconNames = ampliconNames,
#                                          sampleNames = "Normal",removeDup = FALSE)


referenceDirectory <- "~/app/Normal_data/"
referenceFilenames <- list.files(path = referenceDirectory,pattern = "bqsr.bam$", full.names = TRUE)
removePcrDuplicates <- FALSE
referenceReadCounts <- ReadCountsFromBam(referenceFilenames,genomicRangesFromBed,sampleNames = referenceFilenames,ampliconNames = ampliconNames,removeDup = removePcrDuplicates)



# sampleReadCounts <-  ReadCountsFromBam("~/Desktop/Javitri.realigned.bam",
#                                        gr = genomicRangesFromBed,ampliconNames = ampliconNames,
#                                        sampleNames = "Tumor",removeDup = FALSE)

sampleReadCounts <-  ReadCountsFromBam(paste0(args[2],".bqsr.bam"),
                                       gr = genomicRangesFromBed,ampliconNames = ampliconNames,
                                       sampleNames = "Tumor",removeDup = FALSE)

normalizedReadCounts <-  CombinedNormalizedCounts(sampleReadCounts,referenceReadCounts,ampliconNames = ampliconNames)
# After normalization data sets need to be splitted again to perform bootstrap
samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]
print(dim(samplesNormalizedReadCounts))
print(dim(referenceNormalizedReadCounts))

bootList <- BootList(geneNames,samplesNormalizedReadCounts,referenceNormalizedReadCounts,replicates = 10000)
print(summary(bootList))
backgroundNoise <-  Background(geneNames,samplesNormalizedReadCounts,referenceNormalizedReadCounts,bootList,
                               replicates = 10000,significanceLevel = 0.001)
print(summary(backgroundNoise))
reportTables <-  ReportTables(geneNames,samplesNormalizedReadCounts,referenceNormalizedReadCounts,bootList,
                              backgroundNoise)
print(reportTables)
summary(reportTables)
cnvs =as.data.frame(reportTables)
print(dim(cnvs))
names(cnvs)
cnvs = cnvs[cnvs$Tumor.Passed != 0,]
cnvs$Tumor.Gene = rownames(cnvs)
print("good till now")
#cnvs = cnvs
#cnvs = cnvs[,c(13,1,11,12)]
cnvs = cnvs[,c(14,1,11,12)]


#writing table to file
file=args[2]
fname=gsub("_T","",file)
write.csv(cnvs, file = paste0(fname,".copyNumberVariations.csv"),row.names=FALSE)
#write.csv(cnvs, file = "/home/shruti/samples/OMPRAKASH/cnv.csv",row.names=FALSE)
