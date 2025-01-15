####################################################################################################################################
#					TEST PIPELINE
####################################################################################################################################


####################################################################################################################################
#					FILES & DIRECTORIES
###################################################################################################################################

ref=~~/app/hg38/hg38.fasta
bedfile=~/app/hg38/XT-HS_coverage_hg38.bed
genesBed=~/app/hg38/XT-HS_sorted_hg38.bed
annovar=~/app/Reference_hg38/annovar/
snpEff=~/app/snpEff/snpEff.jar


echo "ANALYSIS STARTED" >>output.log
date +"%D %H:%M:%S" >>output.log

echo "patient data initialisation" >>output.log
patientID=$1
workDir=$PWD
cd ${workDir}

mkdir raw
cd raw
ln -s ../*gz .
cd ../

#tumor Sample
rawFqTR1=$2
rawFqTR2=$3

#normal matched
rawFqNR1=$4
rawFqNR2=$5
####################################################################################################################################
#					PARAMETERS FOR CUSTOMISING PROCESSING AND ANALYSIS
####################################################################################################################################

echo "initialising parameters for analysis">>output.log
threads=24

# parameters for qc

phredScore=20 #min phred score while trimming the reads.
library=2 #Primer/Adaptor libraries for NGS QC: 1 = Genomic DNA/Chip-seq Library;2 = Paired End DNA Library;3 = DpnII gene expression Library;4 = NlaIII gene;expression Library;5 = Small RNA Library;6 = Multiplexing DNA Library;N = Do not filter for Primer/Adaptor
fastQcVar=A #FASTQ variants: 1 = Sanger (Phred+33, 33 to 73); 2 = Solexa (Phred+64, 59 to 104); 3 = Illumina (1.3+) (Phred+64, 64 to 104); 4 = Illumina (1.5+) (Phred+64, 66 to 104); 5 = Illumina (1.8+) (Phred+33, 33 to 74); A = Automatic detection of FASTQ variant


# parameters for alignment
sortOrder=coordinate #parameter of picard tool to mark duplicates. Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate}
rdgplLib=library ##parameter of picard tool to mark duplicates.Read Group library 
rdGpPtf=illumina #parameter of picard tool to mark duplicates. Read Group platform (e.g. illumina, solid) Required. 
rdGpPtU=HaloPlex #Read Group platform unit (eg. run barcode) Required. 
compLev=0 #exclusive for all picard tools. Compression level for all compressed files created (e.g. BAM and GELI). Default value: 5. This option can be set to null to clear the default value. 
createIndex=true #exclusive for all picard tools.Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. This option can be set to null to clear the default value. Possible values: {true, false} 
valStrategy=LENIENT # exclusive for all picard tools. Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT, LENIENT, SILENT} 
#insertSize=500

####################################################################################################################################
#					QUALITY CHECK - DONE
####################################################################################################################################
echo "Trimming adapters and removing adapters" >>output.log
date +"%D %H:%M:%S" >>output.log
mkdir clean


fastp -h clean/${patientID}_T.fastp.html -j clean/${patientID}_T.fastp.json -w 16 -l 30 --cut_mean_quality 20 -a auto -f 1 -F 1 -p -z 9 -i ${rawFqTR1} -I ${rawFqTR2} -o clean/${patientID}_T_1.clean.fastq.gz -O clean/${patientID}_T_2.clean.fastq.gz &

fastp -h clean/${patientID}_N.fastp.html -j clean/${patientID}_N.fastp.json -w 16 -l 30 --cut_mean_quality 20 -a auto -f 1 -F 1 -p -z 9 -i ${rawFqNR1} -I ${rawFqNR2} -o clean/${patientID}_N_1.clean.fastq.gz -O clean/${patientID}_N_2.clean.fastq.gz &
wait


echo "Checking fastQ data quality"
cd clean
for i in *gz; do echo "fastqc "$i" "; done| parallel -j2

echo "Trimming done using fastp" >>../output.log
date +"%D %H:%M:%S" >>../output.log
#multiqc -p .

echo “Starting Analysis of QC” >>../output.log
date +"%D %H:%M" >>../output.log

for i in *.clean.fastq.gz; do echo "pigz -ckdp24 "$i">$(basename "$i" .clean.fastq.gz).clean.fastq";done| parallel -j4

wait

~/app/IlluQC_PRLL.pl -pe ${patientID}_T_1.clean.fastq  ${patientID}_T_2.clean.fastq N A -c 24 -s 20 -onlyStat
~/app/IlluQC_PRLL.pl -pe ${patientID}_N_1.clean.fastq  ${patientID}_N_2.clean.fastq N A -c 24 -s 20 -onlyStat
wait	
echo “End Analysis of QC” >>../output.log
date +"%D %H:%M" >>../output.log

echo “End Analysis of QC” >>../output.log
date +"%D %H:%M" >>../output.log


####################################################################################################################################
#					ALIGNMENT AND RECALIBRATION - DONE
####################################################################################################################################


echo "Aligning the files to Ref (hg38.fasta) genome" >>../output.log
date +"%D %H:%M" >>../output.log

for i in *_1.clean.fastq.gz; do echo "bwa mem ${ref} "$i" $(basename "$i" _1.clean.fastq.gz)_2.clean.fastq.gz -M -t 24 | samtools view -Sbh /dev/stdin | samtools "sort" -@ 15 /dev/stdin> ../$(basename "$i" _1.clean.fastq.gz).alnpe.sort.bam";done | parallel -j2

cd ../

echo "Alignment done" >>../output.log
date +"%D %H:%M" >>../output.log


echo "Indexing bam files" >>output.log
date +"%D %H:%M" >>output.log

echo "indexing bam file"
for i in *sort.bam; do echo "samtools index "$i"" ; done| parallel -j2

echo "Indexing bam files done" >>output.log
date +"%D %H:%M" >>output.log

echo "Adding RG to bam files" >>output.log
date +"%D %H:%M" >>output.log

echo "adding read groups" 

for i in *sort.bam; do echo " picard AddOrReplaceReadGroups I="$i" O=$(basename "$i" sort.bam)sortrg.bam RGLB=library RGPL=illumina RGPU=HaloPlex RGSM=$(basename "$i" .alnpe.sort.bam) COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT SO=coordinate CREATE_INDEX=true"; done| parallel -j 2


echo "Adding RG to bam files done" >>output.log
date +"%D %H:%M" >>output.log

echo "Marking duplacates in bam files" >>output.log
date +"%D %H:%M" >>output.log


echo "Marking duplicates"

for i in *.sortrg.bam; do echo "picard MarkDuplicates I="$i"  O=$(basename "$i" .sortrg.bam).sort.markdup.bam  REMOVE_DUPLICATES=true AS=true METRICS_FILE=$(basename "$i" .sortrg.bam).sort.markdup.metrics" ; done| parallel -j2

for i in *.sort.markdup.bam; do echo "samtools index -b "$i""; done | parallel -j2
echo "Marking duplacates in bam files done" >>output.log
date +"%D %H:%M" >>output.log

echo "BQSR of bam files" >>output.log
date +"%D %H:%M" >>output.log


for i in *.sort.markdup.bam; do echo "gatk BaseRecalibrator -I "$i" -R ${ref} --known-sites ~/app/Mills_and_1000G_gold_standard.indels.hg38.vcf --known-sites  ~/app/dbsnp_146.hg38.vcf.gz -O $(basename "$i" .sort.markdup.bam).recal_data.table"; done | parallel -j2



for i in *.sort.markdup.bam; do echo " gatk ApplyBQSR -R ${ref} -I "$i" --bqsr-recal-file $(basename "$i" .sort.markdup.bam).recal_data.table -O $(basename "$i" .alnpe.sort.markdup.bam).bqsr.bam" ; done| parallel -j2

echo "BQSR of bam files done" >>output.log
date +"%D %H:%M" >>output.log


echo "running QC post alignment and generating final QC report" >>output.log
date +"%D %H:%M" >>output.log

Rscript ~/app/QCSummary.R  ${genesBed} ${patientID}_T clean/IlluQC_Filtered_files/${patientID}_T_1.clean.fastq_${patientID}_T_2.clean.fastq_stat
Rscript ~/app/QCSummary.R  ${genesBed} ${patientID}_N clean/IlluQC_Filtered_files/${patientID}_N_1.clean.fastq_${patientID}_N_2.clean.fastq_stat

echo "QUALITY CHECK SUMMARY REPORT GENERATED" >>output.log
date +"%D %H:%M" >>output.log


echo "calling variants with GATK haplotype caller starts " >>output.log
date +"%D %H:%M" >>output.log
~/app/gatk-4.6.0.0/gatk HaplotypeCaller -L ${bedfile}  -R ${ref}  -I ${patientID}_T.bqsr.bam -O ${patientID}.haplotypeCaller_raw.vcf  --dbsnp ~/app/dbsnp_146.hg38.vcf.gz

~/app/gatk-4.6.0.0/gatk SelectVariants -R ${ref} -V ${patientID}.haplotypeCaller_raw.vcf --select-type-to-include SNP  -O ${patientID}.haplotypeCaller_raw_snps.vcf
:'
~/app/gatk-4.6.0.0/gatk VariantFiltration -R ${ref} -V ${patientID}.haplotypeCaller_raw_snps.vcf --filter-expression "QD <2.0 || FS > 60.0 || MQ<40.0 ||MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0 || HaplotypeScore > 13.0" --filter-name ${patientID}.snp.filter -O ${patientID}.haplotypeCaller_filtered_snps.vcf

~/app/gatk-4.6.0.0/gatk SelectVariants  -R ${ref} -V ${patientID}.haplotypeCaller_raw.vcf --select-type-to-include INDEL -O ${patientID}.haplotypeCaller_raw_indels.vcf


~/app/gatk-4.6.0.0/gatk VariantFiltration -R ${ref} -V ${patientID}.haplotypeCaller_raw_indels.vcf --filter-expression "QD <2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name ${patientID}.indel.filter -O ${patientID}.haplotypeCaller_filtered_indels.vcf

'

~/app/gatk-4.6.0.0/gatk VariantFiltration -R ${ref} -V ${patientID}.haplotypeCaller_raw_snps.vcf --filter-expression "(QD <2.0) || (FS > 60.0) || (MQ<40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 4.0) || (HaplotypeScore > 13.0)" --filter-name ${patientID}.snp.filter -O ${patientID}.haplotypeCaller_filtered_snps.vcf

~/app/gatk-4.6.0.0/gatk SelectVariants  -R ${ref} -V ${patientID}.haplotypeCaller_raw.vcf --select-type-to-include INDEL -O ${patientID}.haplotypeCaller_raw_indels.vcf


~/app/gatk-4.6.0.0/gatk VariantFiltration -R ${ref} -V ${patientID}.haplotypeCaller_raw_indels.vcf --filter-expression "(QD <2.0 )|| (FS > 200.0) ||( ReadPosRankSum < -8.0) || (SOR > 10.0)" --filter-name ${patientID}.indel.filter -O ${patientID}.haplotypeCaller_filtered_indels.vcf
#gatk SortVcf -I ${patientID}.haplotypeCaller_filtered_snps.vcf -I ${patientID}.haplotypeCaller_filtered_indels.vcf -O ${patientID}.haplotypeCaller_filtered.vcf

~/app/gatk-4.6.0.0/gatk SortVcf -I ${patientID}.haplotypeCaller_filtered_snps.vcf -I ${patientID}.haplotypeCaller_filtered_indels.vcf -O ${patientID}.haplotypeCaller_filter_marked.vcf

grep -i "^#\|PASS" ${patientID}.haplotypeCaller_filter_marked.vcf > ${patientID}.haplotypeCaller_filtered.vcf



echo "calling variants with GATK haplotype caller done" >>output.log
date +"%D %H:%M" >>output.log



echo "calling variants with GATK Mutect caller starts " >>output.log
date +"%D %H:%M" >>output.log

~/app/gatk-4.6.0.0/gatk Mutect2 -R ${ref} -I ${patientID}_T.bqsr.bam   --germline-resource ~/app/af-only-gnomad.hg38.vcf.gz -L ${bedfile} -pon ~/app/somatic-hg38_1000g_pon.hg38.vcf.gz --f1r2-tar-gz ${patientID}_f1r2.tar.gz -O ${patientID}_mutect_unfiltered.vcf

# pass this raw data to LearnReadOrientationModel
~/app/gatk-4.6.0.0/gatk LearnReadOrientationModel -I ${patientID}_f1r2.tar.gz -O ${patientID}_f1r2_read-orientation-model.tar.gz

# Generate pileup summaries on tumor sample:
~/app/gatk-4.6.0.0/gatk GetPileupSummaries -I ${patientID}_T.bqsr.bam -O ${patientID}.targeted_sequencing.table -V ~/app/af-only-gnomad.hg38.vcf.gz  -L ${bedfile} -R ${ref}

# Calculate contamination on tumor sample
~/app/gatk-4.6.0.0/gatk CalculateContamination -I ${patientID}.targeted_sequencing.table -tumor-segmentation segments.table -O ${patientID}.targeted_sequencing.contamination.table

# Filter variant calls from MuTect
~/app/gatk-4.6.0.0/gatk FilterMutectCalls -R ${ref}  -V ${patientID}_mutect_unfiltered.vcf --tumor-segmentation segments.table --contamination-table ${patientID}.targeted_sequencing.contamination.table --ob-priors ${patientID}_f1r2_read-orientation-model.tar.gz --min-allele-fraction 0.05 -O ${patientID}_mutect_pre_filtered.vcf

grep -i "^#\|PASS" ${patientID}_mutect_pre_filtered.vcf > ${patientID}_T_mutect_filtered.vcf


echo "calling variants with GATK mutect caller done" >>output.log
date +"%D %H:%M" >>output.log

echo "Annotation of variants starts" >>output.log
date +"%D %H:%M" >>output.log
${annovar}convert2annovar.pl -format vcf4  --includeinfo ${patientID}_T_mutect_filtered.vcf >${patientID}_haplotypeCaller_annovar_input.vcf
${annovar}table_annovar.pl ${patientID}_haplotypeCaller_annovar_input.vcf ${annovar}humandb/ --otherinfo -out ${patientID}.allannotated -remove -protocol refGene,cosmic90,clinvar_20191111,dbnsfp35c,ljb26_sift,ljb26_pp2hdiv,ljb26_lrt,ljb26_fathmm,ljb26_mt,snp147,exac03nontcga,ALL.sites.2015_08 -operation g,f,f,f,f,f,f,f,f,f,f,f -build hg38 -nastring . -csvout &
echo "Annotation of variants done" >>output.log
date +"%D %H:%M" >>output.log

echo "Call Somatic Mutation using MuTect2 starts " >>output.log
date +"%D %H:%M" >>output.log

sh ~/Desktop/krypton_share/scripts/mutect2_tumor_paired_normal.sh ${patientID}

date +"%D %H:%M" >>output.log


echo "Annotation of variants starts" >>output.log
date +"%D %H:%M" >>output.log

date +"%D %H:%M" >>output.log

echo "calling variants with SNVer starts" >>output.log
date +"%D %H:%M" >>output.log

java -jar ~/apps/SNVerIndividual.jar -i ${patientID}_T.bqsr.bam -o ${patientID}_T.SNVer -n 2 -r ${ref} -b 0.25 -bq 20 -l ${bedfile}

bgzip ${patientID}_T.SNVer.filter.vcf

bgzip ${patientID}_T.SNVer.indel.filter.vcf

tabix -p vcf  ${patientID}_T.SNVer.filter.vcf.gz

tabix -p vcf ${patientID}_T.SNVer.indel.filter.vcf.gz

vcf-concat ${patientID}_T.SNVer.filter.vcf.gz ${patientID}_T.SNVer.indel.filter.vcf.gz | gzip -c > ${patientID}_T.SNVer.concat.vcf.gz
vcf-sort ${patientID}_T.SNVer.concat.vcf.gz > ${patientID}_T.SNVer.concat.sorted.vcf


echo "calling variants with SNVer done" >>output.log
date +"%D %H:%M" >>output.log
${annovar}convert2annovar.pl --format vcf4 --includeinfo  ${patientID}_T.SNVer.concat.sorted.vcf > ${patientID}_T.SNVer.annovarInput.vcf

${annovar}table_annovar.pl ${patientID}_T.SNVer.annovarInput.vcf ${annovar}humandb/ --otherinfo -out ${patientID}_snver.allannotated -remove -protocol refGene,cosmic90,clinvar_20191111,dbnsfp35c,ljb26_sift,ljb26_pp2hdiv,ljb26_lrt,ljb26_fathmm,ljb26_mt,snp147,exac03nontcga,ALL.sites.2015_08 -operation g,f,f,f,f,f,f,f,f,f,f,f -build hg38 -nastring . -csvout &

echo "calling variants with Vardict"

echo "Extract Somatic and Germline variants"


#Structural variants
echo "SVs calling by breakdancer starts"
~/apps/breakdancer/perl/bam2cfg.pl -g -h ${patientID}_T.bqsr.bam ${patientID}_N.bqsr.bam > ${patientID}.Breakdancer.cfg
breakdancer-max  -q 10 ${patientID}.Breakdancer.cfg > ${patientID}.ctx
echo "SVs calling by breakdancer finished"


echo "SVs calling by LUMPY started"
for i in *bqsr.bam; do echo " samtools view -b -F 1294 "$i"| samtools sort -@ 21 /dev/stdin >$(basename "$i" .bqsr.bam).discordants.sorted.bam" ; done| parallel -j2
for i in *bqsr.bam; do echo " samtools view -h "$i" | extractSplitReads_BwaMem  -i stdin  | samtools view -Sb - > $(basename "$i" .bqsr.bam).splitters.unsorted.bam" ; done | parallel -j2
for i in *splitters.unsorted.bam; do echo " samtools sort -@ 21 "$i" >$(basename "$i" .splitters.unsorted.bam).splitters.sorted.bam"; done | parallel -j2
~/apps/lumpy-sv/bin/lumpyexpress -B ${patientID}_T.bqsr.bam,${patientID}_N.bqsr.bam -S  ${patientID}_T.splitters.sorted.bam,${patientID}_N.splitters.sorted.bam -D ${patientID}_T.discordants.sorted.bam,${patientID}_N.discordants.sorted.bam -o ${patientID}_T_lumpy.vcf

echo "Annotating SV file"
java -Xmx4G -jar ${snpEff} eff -v  -i vcf -o vcf hg38 ${patientID}_T_lumpy.vcf > ${patientID}_T.lumpy.SV.annot.vcf
echo "SVs calling by LUMPY finished"
grep -w "SVTYPE=BND" ${patientID}_T.lumpy.SV.annot.vcf | grep -Fwf ~/Reference_Files/Human_reference/560_gene_lumpy_manta.csv /dev/stdin> ${patientID}.lumpy.SV.annot_BND.vcf

echo "Mutation calling by MANTA"
python2.7 ~/Desktop/krypton_share/apps/manta/build/bin/configManta.py --tumorBam ${patientID}_T.bqsr.bam  --exome --referenceFasta /home/thunder/Tapan/ref/hg38/hg38.fasta --runDir ${workDir}
python2.7 runWorkflow.py  -m local
gunzip  results/variants/tumorSV.vcf.gz

echo "Annotating SV file"
java -Xmx4G -jar ${snpEff} eff -v  -i vcf -o vcf hg38 results/variants/tumorSV.vcf  > ${patientID}T.Manta.SV.annot.vcf

grep -w "SVTYPE=BND" ${patientID}T.Manta.SV.annot.vcf | grep -Fwf ~/Reference_Files/Human_reference/560_gene_lumpy_manta.csv /dev/stdin> ${patientID}.Manta.SV.annot_BND.vcf


p=ID
echo "${patientID}_N control">>file.txt
echo "${patientID}_T tumor">>file.txt
sed 's/ /\t/g' file.txt>samples.tsv
rm file.txt




echo "calling copy number variations"
Rscript   ~/Desktop/krypton_share/scripts/callCNV.R ${genesBed} ${patientID}_T ${patientID}_N

Rscript  ~/Desktop/krypton_share/scripts/clinicallyActionable.R ~/Desktop/krypton_share/Reference_hg38/hg38/drugRecoMatrices.xlsx ${patientID} ~/Desktop/krypton_share/Reference_hg38/hg38/transcript_560.gtf "One"


export VEP_PATH=~/Desktop/krypton_share/apps/ensembl-vep/
export VEP_DATA=~/Desktop/krypton_share/.vep
export PERL5LIB=$VEP_PATH:$PERL5LIB


~/Desktop/krypton_share/apps/ensembl-vep/vep --species homo_sapiens --assembly GRCh38 --offline  --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir ~/Desktop/krypton_share/.vep/ --fasta ~/Desktop/krypton_share/.vep/homo_sapiens/98_GRCh38/Homo_sapiens.GRCh38.dna.toplevel_header_changed.fa.gz --input_file ${patientID}.haplotypeCaller_filtered.vcf  --output_file ${patientID}_T.haplotypeCaller.vep.vcf --polyphen b --af --af_1kg --af_esp --regulatory --custom ~/Desktop/krypton_share/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz,ExAC,vcf,exact,1,AC,AN

perl ~/Desktop/krypton_share/apps/ensembl-vep/mskcc-vcf2maf-5453f80/vcf2maf.pl --vep-path ~/Desktop/krypton_share/apps/ensembl-vep/  --vep-data ~/Desktop/krypton_share/.vep/  --ref-fasta ~/Desktop/krypton_share/.vep/homo_sapiens/98_GRCh38/Homo_sapiens.GRCh38.dna.toplevel_header_changed.fa.gz  --ncbi-build GRCh38 --input-vcf ${patientID}.haplotypeCaller_filtered.vcf  --tumor-id ${patientID}_T --output-maf ${patientID}T.haplotypeCaller.vep.maf


rm clean/*fastq
mkdir Output
mv error.log Output
mv ${patientID}_T.bqsr.bam Output
mv ${patientID}_T.bqsr.bai Output
mv ${patientID}_N.bqsr.bam Output
mv ${patientID}_N.bqsr.bai Output
mv ${patientID}_T.all.cleaned.vcf Output
mv ${patientID}.allannotated.hg38_multianno.csv Output
mv ${patientID}.clinicallyRelevantMutations.csv Output
mv ${patientID}.DrugRecomedations.xlsx Output
mv ${patientID}.TMB.pdf Output
mv ${patientID}.copyNumberVariations.csv Output
mv ${patientID}.structuralVariants.csv Output
mv ${patientID}_MSI_Analysis.txt Output
mv ${patientID}T.haplotypeCaller.vep.maf Output
mv ${patientID}_T.haplotypeCaller.vep.vcf Output
mv ${patientID}.haplotypeCaller_filtered.vcf* Output
#mv ${patientID}T.prefix
mv ${patientID}_msi* Output
mv ${patientID}_T.SNVer.concat.sorted.vcf* Output
mv ${patientID}_snver.allannotated.hg38_multianno.csv Output
mv ${patientID}_T.lumpy.SV.annot.vcf Output
mv ${patientID}_T_lumpy.vcf Output
mv ${patientID}T.Manta.SV.annot.vcf Output
mv ${patientID}T_D.Deletion_pindel.hg38_multianno.csv Output
mv *utect* Output
mv *_BND.vcf *QC_Metrics.txt Output
echo "Creating QC folder and moving files required for analysis"
mkdir QC
mv ${patientID}_T.coverageExon.csv QC
mv ${patientID}_T.coverageGene.csv QC
mv ${patientID}_T.QCSummary.pdf QC
mv ${patientID}_N.coverageExon.csv QC
mv ${patientID}_N.coverageGene.csv QC
mv ${patientID}_N.QCSummary.pdf QC


echo "SUCCESSFULLY COPIED TO POSITIVE_SHARE!!"
echo "ANALYSIS FINISHED" >>output.log
date +"%D %H:%M:%S" >>output.log
