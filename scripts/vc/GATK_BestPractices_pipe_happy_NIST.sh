#!/bin/bash
####### Script for performing GATK best practices Variant Calling ########

### input args
tempRoot=$1
num_threads=$2
fastq1=$3
fastq2=$4
###


##################################
# Modify paths
##################################

## Reference FASTA
ref="/home/human_g1k_v37.fasta"

### Reference VCF files 
hapmapFile="/home/bundle_2.8/hapmap_3.3.b37.vcf"
omniFile="/home/bundle_2.8/1000G_omni2.5.b37.vcf"
KGFile="/home/bundle_2.8/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnpsFile="/home/bundle_2.8/dbsnp_138.b37.vcf"
Mills1KGIndels="/home/bundle_2.8/Mills_and_1000G_gold_standard.indels.b37.vcf"
KGFileIndels="/home/bundle_2.8/1000G_phase1.indels.b37.vcf"

### Ground Truth
gt_vcf="/home/NIST_dataset/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.vcf.gz"
gt_bed="/home/NIST_dataset/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed"

### Software tools
picard_program="/home/picard-tools-2.4.1/picard.jar"
genomeanalysistk_program="/home/GenomeAnalysisTK.jar"
happy_program="/home/hap.py/bin/hap.py"
reppy_program="/home/bin/rep.py"

##################################

### Haplotype Caller Sensitivity
SEC=10
SCC=30
###

### VQSR Parameters
resource1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmapFile"
resource2="omni,known=false,training=true,truth=true,prior=12.0 $omniFile"
resource3="1000G,known=false,training=true,truth=false,prior=10.0 $KGFile"
resource4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnpsFile"
resourceIndels="mills,known=true,training=true,truth=true,prior=12.0 $Mills1KGIndels" 
recalParams="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum"
recalParamsIndels="-an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum"
filter_level="99.0"
###

### temp files
aln_SAM=$tempRoot".aln.sam"
aln_BAM=$tempRoot".aln.bam"
aln_sorted_BAM=$tempRoot".aln.sorted.bam"
aln_sorted_dedup_BAM=$tempRoot".aln.sorted.dedup.bam"
aln_sorted_dedup_rg_BAM=$tempRoot".aln.sorted.dedup.rg.bam"
aln_sorted_dedup_rg_recal_BAM=$tempRoot".aln.sorted.dedup.rg.recal.bam"
dedupMetrics=$tempRoot".dedupMetrics.txt"
bqsr_data=$tempRoot".bqsr.table"
raw_VCF=$tempRoot".raw.vcf"
recalFile=$tempRoot".vqsr.recal"
tranchesFile=$tempRoot".vqsr.tranches"
rscriptFile=$tempRoot".vqsr.R"
recalFileIndels=$tempRoot".vqsr.indels.recal"
tranchesFileIndels=$tempRoot".vqsr.indels.tranches"
rscriptFileIndels=$tempRoot".vqsr.indels.R"
recal_snps_VCF=$tempRoot".vqsr.recal.snps.vcf"
recal_snps_indels_VCF=$tempRoot".vqsr.recal.snps.indels.vcf"
###


# BWA MEM alignment
bwa mem -t $num_threads -M $ref $fastq1 $fastq2 > $aln_SAM
# Convert SAM to BAM
samtools view -@ $num_threads -b -h $aln_SAM > $aln_BAM
# Sort File using samtools
samtools sort -T ./tmp -@ $num_threads -O bam $aln_BAM > $aln_sorted_BAM
# Mark Duplicates in the input file (using Picard tools)
java -jar -Djava.io.tmpdir=./tmp $picard_program MarkDuplicates I=$aln_sorted_BAM O=$aln_sorted_dedup_BAM M=$dedupMetrics ASSUME_SORTED=true
# Add read group name
java -jar -Djava.io.tmpdir=./tmp $picard_program AddOrReplaceReadGroups I=$aln_sorted_dedup_BAM O=$aln_sorted_dedup_rg_BAM RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=NA12878
# We need to index the BAM file
java -jar -Djava.io.tmpdir=./tmp $picard_program BuildBamIndex I=$aln_sorted_dedup_rg_BAM

# Recalibrate Quality Scores
java -jar -Djava.io.tmpdir=./tmp $genomeanalysistk_program -nct $num_threads -T BaseRecalibrator -R $ref -I $aln_sorted_dedup_rg_BAM -knownSites $dbsnpsFile -knownSites $Mills1KGIndels -knownSites $KGFileIndels -o $bqsr_data
java -jar -Djava.io.tmpdir=./tmp $genomeanalysistk_program -nct $num_threads -T PrintReads -R $ref -I $aln_sorted_dedup_rg_BAM -BQSR $bqsr_data -o $aln_sorted_dedup_rg_recal_BAM

# Call variants using Haplotype Caller
java -jar -Djava.io.tmpdir=./tmp $genomeanalysistk_program -T HaplotypeCaller -R $ref -I $aln_sorted_dedup_rg_recal_BAM --genotyping_mode DISCOVERY -stand_emit_conf $SEC -stand_call_conf $SCC -o $raw_VCF

# Fiter the variants using VQSR
java -jar -Djava.io.tmpdir=./tmp $genomeanalysistk_program -R $ref -T VariantRecalibrator -input $raw_VCF -resource:$resource1 -resource:$resource2 -resource:$resource3 -resource:$resource4 $recalParams -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $recalFile -tranchesFile $tranchesFile -rscriptFile $rscriptFile
java -jar -Djava.io.tmpdir=./tmp $genomeanalysistk_program -R $ref -T ApplyRecalibration -input $raw_VCF -mode SNP -recalFile $recalFile -tranchesFile $tranchesFile --ts_filter_level $filter_level -o $recal_snps_VCF

java -jar -Djava.io.tmpdir=./tmp $genomeanalysistk_program -R $ref -T VariantRecalibrator -input $recal_snps_VCF -resource:$resourceIndels $recalParamsIndels -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $recalFileIndels -tranchesFile $tranchesFileIndels -rscriptFile $rscriptFileIndels
java -jar -Djava.io.tmpdir=./tmp $genomeanalysistk_program -R $ref -T ApplyRecalibration -input $recal_snps_VCF -mode INDEL -recalFile $recalFileIndels -tranchesFile $tranchesFileIndels --ts_filter_level $filter_level -o $recal_snps_indels_VCF

######################
happy_root=$tempRoot".happy.nist"
out_html=$tempRoot".happy.nist.html"
# run hap.py
python $happy_program $gt_vcf $recal_snps_indels_VCF -f $gt_bed -o $happy_root -r $ref --roc VQLSOD

#run the benchmarking util rep.py
rm -r "$happy_root".rep.tsv
printf "method\tcomparisonmethod\tfiles\n" >> "$happy_root".rep.tsv
printf "gatk3\tNIST\t" >> "$happy_root".rep.tsv
printf "$happy_root.extended.csv," >> "$happy_root".rep.tsv
for i in "$happy_root".roc.Locations.*; do printf "$i,"; done >> "$happy_root".rep.tsv
sed -i '$ s/.$//' "$happy_root".rep.tsv
python $reppy_program -o $out_html -l "$happy_root".rep.tsv
