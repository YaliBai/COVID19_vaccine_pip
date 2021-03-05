#!/usr/bin/bash
bcftools='~/./tools/anaconda2/bin/bcftools'
samtools='~/./tools/anaconda2/bin/samtools'
ref_fa='~/./database/human/hg38/WholeGenomeFasta/hg38.fa
ref_snps='~/./hg38_snp153_rmi.bed'

vcftools=~/./tools/anaconda2/bin/vcftools
demuxlet=/share/soft/demuxlet/bin/demuxlet

bamdir=$1
MalbacBamdir=$2
outdir=$3


set -x


function preprocess {
  mkdir -p $outdir/logs
  mkdir -p $outdir/malbac

  ls $MalbacBamdir/*/*.bam|while read bamin
  do

  id=`echo $bamin|sed 's/.*\///g'|sed 's/_human_Aligned.sortedByCoord.out.bam//g'`


sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J $id
#SBATCH --partition=compute_new
#SBATCH --exclude=node10,node17
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH -o $outdir/logs/${id}.callsnp.log
#SBATCH -e $outdir/logs/${id}.callsnp.log
set -x
echo $outdir/malbac


$samtools sort -@ 6 -n -o $outdir/malbac/${id}.sortn.bam $bamin
$samtools fixmate -@ 6 -m $outdir/malbac/${id}.sortn.bam $outdir/malbac/${id}.fixed.bam
$samtools sort -@ 6 -o $outdir/malbac/${id}.sorted.bam $outdir/malbac/${id}.fixed.bam
$samtools markdup -@ 6 -r  $outdir/malbac/${id}.sorted.bam $outdir/malbac/${id}.rmdup.bam
$samtools index $outdir/malbac/${id}.rmdup.bam

RUN
done

}


function ref_snps {
# 切分跑

less  -S SNPs/hg38_snp153_rmi.new.bed |cut -f -3 > SNPs/hg38_snp153_rmi.clean.bed
bedtools split -i SNPs/hg38_snp153_rmi.clean.bed --prefix SNPs/hg38_snp153_rmi.cut -a simple  -n 20

ls SNPs/hg38_snp153_rmi.cut.*bed |while read cut_ref
do
  cut_name=`echo $cut_ref|sed 's/.*hg38_snp153_rmi.//g'|sed 's/.bed//g'`

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J $cut_name
#SBATCH --partition=compute_new
#SBATCH --exclude=node10,node17
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH -o $outdir/logs/${prefix}.${cut_name}.log
#SBATCH -e $outdir/logs/${prefix}.${cut_name}.log
set -x

$bcftools mpileup \
  --regions-file $cut_ref \
  ${MalbacBamdir}/*.rmdup.bam \
  --fasta-ref $ref_fa \
  --output ${outdir}/ref_snps/${prefix}.${cut_name}.snp.vcf

$bcftools call \
  -c \
  ${outdir}/ref_snps/${prefix}.${cut_name}.snp.vcf \
  -o ${outdir}/ref_snps/${prefix}.${cut_name}.snp_call.vcf

$bcftools filter -i 'ALT!="." && INFO/DP>10 && QUAL >30' ${outdir}/ref_snps/${prefix}.${cut_name}.snp_call.vcf -O v --output ${outdir}/ref_snps/${prefix}.${cut_name}.snp_call_filter.vcf

bgzip ${outdir}/ref_snps/${prefix}.${cut_name}.snp_call_filter.vcf
$bcftools index ${outdir}/ref_snps/${prefix}.${cut_name}.snp_call_filter.vcf.gz


set +x
RUN

done

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J chrM
#SBATCH --partition=compute_new
#SBATCH --exclude=node10,node17
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH -o $outdir/logs/${prefix}.chrM.log
#SBATCH -e $outdir/logs/${prefix}.chrM.log
set -x

$bcftools mpileup \
  --regions chrM \
  ${MalbacBamdir}/*.rmdup.bam \
  --fasta-ref $ref_fa \
  --output ${outdir}/ref_snps/${prefix}.chrM.snp.vcf

$bcftools call \
  -c \
  ${outdir}/ref_snps/${prefix}.chrM.snp.vcf \
  -o ${outdir}/ref_snps/${prefix}.chrM.snp_call.vcf

  $bcftools filter -i 'ALT!="." && INFO/DP>10 && QUAL >30' ${outdir}/ref_snps/${prefix}.chrM.snp_call.vcf -O v --output ${outdir}/ref_snps/${prefix}.chrM.snp_call_filter.vcf

  bgzip ${outdir}/ref_snps/${prefix}.chrM.snp_call_filter.vcf
  $bcftools index ${outdir}/ref_snps/${prefix}.chrM.snp_call_filter.vcf.gz


set +x
RUN

}

function merge_snp {
new_order="chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr2,chr20,chr21,chr22,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chrM,chrX,chrY"


zcat  ${outdir}/ref_snps/${prefix}.cut.000*.snp_call_filter.vcf.gz ${outdir}/ref_snps/${prefix}.chrM.snp_call_filter.vcf.gz |grep -v "#" > ${outdir}/ref_snps/${prefix}.snp_call_filter.vcf.tmp
less ${outdir}/ref_snps/${prefix}.chrM.snp_call_filter.vcf.gz|grep "#" > ${outdir}/ref_snps/${prefix}.snp_call_filter.vcf.header
cat ${outdir}/ref_snps/${prefix}.snp_call_filter.vcf.header ${outdir}/ref_snps/${prefix}.snp_call_filter.vcf.tmp > ${outdir}/ref_snps/${prefix}.snp_call_filter.vcf

$bcftools sort ${outdir}/ref_snps/${prefix}.snp_call_filter.vcf -O v -o ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted.vcf


bgzip ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted.vcf
bcftools index ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted.vcf.gz

bcftools view \
  -r $new_order \
  ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted.vcf.gz \
  -O v \
  -o ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted2.vcf

new_cc=`echo $new_order|sed 's/,/\n/g'|while read cc; do less ${outdir}/ref_snps/${prefix}.snp_call_filter.vcf|grep "##contig=<ID=${cc},"; done`# > ${outdir}/ref_snps/${prefix}.snp_call_filter.contig.vcf

##less ${outdir}/ref_snps/${prefix}.snp_call_filter.vcf|grep -v "#contig=<ID=[EKG]"| awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted2.vcf
less ${outdir}/ref_snps/${prefix}.snp_call_filter.vcf|grep -v "#contig=<ID=[EKG]"| awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' |sed "s#${MalbacBamdir}/##g"|sed 's/.rmdup.bam//g'  > ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted2.vcf
grep -v "##contig=" ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted2.vcf > ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted.vcf

}

function sep_snps {
set_prefix=$1

/share/soft/demuxlet/bin/demuxlet --sam $bamdir/possorted_genome_bam.bam \
  --vcf ${outdir}/ref_snps/${prefix}.snp_call_filter.sorted.vcf \
  --field GT \
  --out ${outdir}/sep_files/${set_prefix}
}

bamdir=$1
MalbacBamdir=$2
outdir=$3


#preprocess
MalbacBamdir=malbac

bamdir=01-cellranger/PR10_XXL_10X-RBD_1/outs/count
outdir=Seurat/PR10_XXL_10X-RBD_1/PhaseSNPs
prefix=SF_refSnps

mkdir -p $outdir
mkdir -p $outdir/logs
mkdir -p $outdir/ref_snps


ref_snps

merge_snp

phased_snps





