## Step 0. Build your DB
#### Database was already build and it is found on /home/cesar/TEST/hg38.fa
#### Database was built with BWA

## Step 1.- Align reads against DB

bwa mem -t 4 -R "@RG\tID:SON\tPL:ILLUMINA\tSM:SON"  /home/cesar/TEST/hg38.fa /data/gene5150/SON/SON_R1.fastq.gz /data/gene5150/SON/SON_R2.fastq.gz > /data/gene5150/Aguirre.Son/Aguirre.Son.bam


### Step 2.- Sort reads

samtools sort -o /data/gene5150/Aguirre.Son/son.sorted.paired.sam /data/gene5150/Aguirre.Son/Aguirre.Son.bam

### Step 3.- Mark Duplicates

/home/cesar/TEST/gatk-4.3.0.0/gatk MarkDuplicates -I /data/gene5150/Aguirre.Son/son.sorted.paired.sam -O /data/gene5150/Aguirre.Son/son.sorted.dedup.paired.sam -M metrics

### Step 4 .- Base calibrator 

/home/cesar/TEST/gatk-4.3.0.0/gatk BaseRecalibrator -I /data/gene5150/Aguirre.Son/son.sorted.dedup.paired.sam  \
-R /home/cesar/TEST/hg38.fa --known-sites /home/cesar/TEST/Homo_sapiens_assembly38.dbsnp138.vcf \
-O recal_data.table

### Step 5.- Apply recalibration.

/home/cesar/TEST/gatk-4.3.0.0/gatk ApplyBQSR -I /data/gene5150/Aguirre.Son/son.sorted.dedup.paired.sam \
-R /home/cesar/TEST/hg38.fa --bqsr-recal-file recal_data.table \
-O /data/gene5150/Aguirre.Son/son_sorted_bqsr_dedup_reads.bam

### Step 6.- Collect metrics.

/home/cesar/TEST/gatk-4.3.0.0/gatk CollectAlignmentSummaryMetrics \
R=/home/cesar/TEST/hg38.fa I=/data/gene5150/Aguirre.Son/son_sorted_bqsr_dedup_reads.bam \	
O=alignment_metrics.txt

/home/cesar/TEST/gatk-4.3.0.0/gatk CollectInsertSizeMetrics \
INPUT=/data/gene5150/Aguirre.Son/son_sorted_bqsr_dedup_reads.bam \
OUTPUT=insert_size_metrics.txt \
HISTOGRAM_FILE=histogram.pdf

### Step 7.- MAIN EVENT! Variant Calling 

home/cesar/TEST/gatk-4.3.0.0/gatk HaplotypeCaller -R /home/cesar/TEST/hg38.fa \
-I /data/gene5150/Aguirre.Son/son_sorted_bqsr_dedup_reads.bam -O raw_variants.vcf

###########
######################################################

### Step 8.- 
### VCF processing :

/home/cesar/TEST/gatk-4.3.0.0/gatk SelectVariants -R /home/cesar/TEST/hg38.fa \
-V /home/cesar/DATASETS/ALL/son_raw_variants.vcf \
--select-type SNP -O raw_snps.vcf

##########

### Step 9.-
### Filtering 

/home/cesar/TEST/gatk-4.3.0.0/gatk VariantFiltration \
-R /home/cesar/TEST/hg38.fa \
-V raw_snps.vcf \
-O filtered_snps.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 60.0" \
-filter-name "MQ_filter" -filter "MQ < 40.0" \
-filter-name "SOR_filter" -filter "SOR > 4.0" \
-genotype-filter-expression "DP < 10" \
-genotype-filter-name "DP_filter" \
-genotype-filter-expression "GQ < 30" \
-genotype-filter-name "GQ_filter"

/home/cesar/TEST/gatk-4.3.0.0/gatk VariantFiltration \
-R /home/cesar/TEST/hg38.fa \
-V raw_snps.vcf \
-O filtered_indels.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 60.0" \
-filter-name "SOR_filter" -filter "SOR > 4.0" \
-genotype-filter-expression "DP < 10" \
-genotype-filter-name "DP_filter" \
-genotype-filter-expression "GQ < 30" \
-genotype-filter-name "GQ_filter"

#### Step 10.- Select variants based on QC

/home/cesar/TEST/gatk-4.3.0.0/gatk SelectVariants \
--exclude-filtered \
-V filtered_snps.vcf \
-O analysis-ready-snps.vcf

/home/cesar/TEST/gatk-4.3.0.0/gatk SelectVariants \
--exclude-filtered \
-V filtered_indels.vcf \
-O analysis-ready-indels.vcf

#### Step 11. Annotation 

/home/cesar/TEST/gatk-4.3.0.0/gatk Funcotator --variant analysis-ready-snps.vcf \
--reference /home/cesar/TEST/hg38.fa --ref-version hg38 --data-sources-path  /home/cesar/TEST/funcotator_dataSources.v1.7.20200521g \
--output analysis-ready-snps-filteredGT-functotated.vcf --output-file-format VCF

/home/cesar/TEST/gatk-4.3.0.0/gatk Funcotator --variant analysis-ready-indels.vcf \
--reference /home/cesar/TEST/hg38.fa --ref-version hg38 --data-sources-path  /home/cesar/TEST/funcotator_dataSources.v1.7.20200521g \
--output analysis-ready-snps-filteredGT-functotated.vcf --output-file-format VCF

#### Step 11.5 Alternative Annotation 

grep '^chr' analysis-ready-snps-filteredGT-functotated.vcf | awk '{print $1":"$2"\t"$8}' | sed 's/FUNCOTATION=\[/\t/g' | sed 's/|/\t/g' | awk '{print $1"\t"$3"\t"$8}' > son_annotation.tsv
 
grep '^chr' analysis-ready-snps-filteredGT-functotated.vcf | awk '{print $1,$2,$4,$5,$6}' > First_Cols

grep '^chr' analysis-ready-snps-filteredGT-functotated.vcf | awk '{print $1,$2,$4,$5,$6,$7,$8,$9}' > First_Cols_2
