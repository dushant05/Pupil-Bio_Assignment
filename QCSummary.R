
###################################################COVERAGE##############################################

#source("https://bioconductor.org/biocLite.R")
#biocLite("TEQC")
library(TEQC)
#install.packages("gridExtra","data.table","splitstackshape")
#install.packages("data.frame")
library(gridExtra)
library(data.table)
#library(data.frame)
library(splitstackshape)
#reading the arguments

args = commandArgs(TRUE)

targets = get.targets(targetsfile =args[1] , chrcol = 1, startcol = 2,
                      endcol = 3, zerobased = TRUE, sep = "\t", skip = 0, header = FALSE)

reads = get.reads(readsfile = paste0(args[2],".bqsr.bam"),filetype = "bam" , chrcol = 1, startcol = 2, endcol = 3,
                  idcol, zerobased = TRUE)
#creating the target summary coverage
targetCoverageSummary = coverage.target(reads, targets, Offset = 0, perTarget = TRUE, perBase = TRUE)
onTarget = fraction.reads.target(reads, targets, Offset = 0, mappingReads = FALSE)
coverage = as.data.frame(targetCoverageSummary$targetCoverages)
coverage$avgCoverage = round(coverage$avgCoverage, digits = 0)
Ychromosome = subset(coverage, coverage$space == "chrY")
avgCoverage = mean(coverage$avgCoverage)
medianCoverage = median(coverage$avgCoverage)
summaryCov = as.data.frame(quantile(coverage$avgCoverage,seq(0,1,0.1)))
colnames(summaryCov) = "Coverage"
summaryCov$Percentile = rownames(summaryCov)
summaryCov =summaryCov[c("Percentile", "Coverage")]

coverage$key = paste(coverage$seqnames,coverage$start-1,coverage$end, sep = "-")
#gene = read.delim(file = "/home/shruti/refGenome/100.bed",sep = "\t", header = F)
gene = read.delim(file = args[1],sep = "\t", header = F)
gene$key = paste(gene$V1,gene$V2,gene$V3, sep = "-")

coverage = merge(coverage , gene[,4:5], by = "key")
genecoverage = by(coverage[,c(7,8)], coverage$V4, colMeans)
genecoverage <- as.data.frame(do.call(rbind,genecoverage))
genecoverage$avgCoverage = round(genecoverage$avgCoverage, digits = 0)
genecoverage$gene = rownames(genecoverage)
#writing gene coverage
write.csv(genecoverage[,c(3,1)], file = paste0(args[2],".coverageGene.csv"),row.names=FALSE)

#writing exon coverage
write.csv(coverage[,c(2,3,4,8,6)], file = paste0(args[2],".coverageExon.csv"),row.names=FALSE)

#writing Y chromosome coverage
write.csv(Ychromosome, file = paste0(args[2],".coverageChrY.csv"),row.names=FALSE)


#################################################NGS QC Output##############################################
#ngsqc = read.delim("/home/shruti/samples/PB_CG_WE_2016_67/IlluQC_Filtered_files/Raghu_S4_R1_001.fastq.gz_Raghu_S4_R2_001.fastq.gz_stat",sep = "\t")
ngsqc = read.delim(file = args[3],sep = "\t")

###info 
info = data.table(ngsqc[c(1,3,4,5),])
colnames(info) = "Basic Information"
####for read pairs quality
reads = data.table(ngsqc[c(10,12),])
reads = cSplit(reads, "V1",sep = " ")
reads$Parameter = paste(reads$V1_1,reads$V1_2,reads$V1_3,reads$V1_4,sep = " ")
reads$R1 = reads$V1_5
reads$R2 = reads$V1_6
reads = reads[,7:9]

readBase = data.table(ngsqc[c(28,30,32),])
readBase = cSplit(readBase, "V1",sep = " ")
readBase$Parameter = paste(readBase$V1_1,readBase$V1_2,readBase$V1_3,readBase$V1_4,sep = " ")
readBase$R1 = readBase$V1_5
readBase$R2 = readBase$V1_6
readBase = readBase[,7:9]

##read length
readLength = data.table(ngsqc[22:24,])
readLength = cSplit(readLength, "V1",sep = " ")
readLength$Parameter = paste(readLength$V1_1,readLength$V1_2,readLength$V1_3,sep = " ")
readLength$R1 = readLength$V1_4
readLength$R2 = readLength$V1_5
readLength = readLength[,6:8]

readSummary = rbind(reads,readBase,readLength)
readSummary$R1 = as.integer(sub("%","",readSummary$R1))
readSummary$R2 = as.integer(sub("%","",readSummary$R2))

readSummary$Remark = apply(readSummary[,2:3],1,mean)

#preparing sample QC summary report
coverage = if (medianCoverage >=700) {"PASS"} else {"FAIL"}
ontargetReads = if (onTarget >= .8) {"PASS"} else {"FAIL"}
HQBase = if (readSummary[4,4] >= 90) { "PASS"} else {"FAIL"}
HQReads = if (readSummary[4,4] >= 85) { "PASS"} else {"FAIL"}
avgLength = if (readSummary[8,4] >= 100) { "PASS"} else {"FAIL"}
NBase = if (readSummary[5,4] <= 2) { "PASS"} else {"FAIL"}



sampleQCSummary = as.data.frame(matrix(c("Median Coverage","On Target Reads","High Quality Bases", "High Quality Reads", "Non ATGC Bases", "Average Read Length", 
                            "700X", "90%","90%","85%","2%","100 bp",
                            coverage, ontargetReads,HQBase, HQReads, NBase,avgLength),nrow = 6,ncol = 3))
colnames(sampleQCSummary)<- c("Parameter","Criterion","Remark")

readSummary$Remark = NULL
#preparing a QC report
#filenamePdf = paste0(args[2],".QCSummary.pdf")
#pdf("/home/shruti/samples/Hrangthan/Hrangthan.QCSummary1.pdf")
#pdf("/home/shruti/samples/kartik/kartik3/kartik.QCSummary.pdf")
pdf(paste0(args[2],".QCSummary.pdf"),height=11, width=8.5)
g1 = tableGrob(info,rows=NULL,theme = ttheme_default(base_size =9))
g2 = tableGrob(sampleQCSummary,rows=NULL,theme = ttheme_default(base_size =9))
g3 = tableGrob(readSummary,rows=NULL,theme = ttheme_default(base_size =9))
g4 = tableGrob(summaryCov,rows=NULL,theme = ttheme_default(base_size =9))

grid.arrange(g1,g2,g3,g4,nrow =2, ncol =2)
coverage.hist(targetCoverageSummary$coverageTarget, col.hist = "lightblue", col.line = "orange", 50, breaks = "Sturges")
#insert.size.hist(readpairs, returnInserts = FALSE, legendpos="topleft", main, xlab, ylab)
#grid.arrange(g4,g1,nrow =2)
dev.off()
