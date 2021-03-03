#!/usr/bin/bash

trim_galore='~/./TrimGalore-0.4.5/trim_galore'
STAR='~/./STAR-2.5.3a/STAR/bin/Linux_x86_64/STAR'
genome='~/./database/human/hg38/Index/STAR'
samtools='~/./tools/anaconda2/bin/samtools'

function run_pip {
indir=$1
id=$2
outdir=$3
mkdir -p $outdir/trim
mkdir -p $outdir/star
mkdir -p $outdir/logs

Read1=$indir/${id}.R1.fastq.gz
Read2=$indir/${id}.R2.fastq.gz

trim_r1=$outdir/trim/${id}.R1_val_1.fq.gz
trim_r2=$outdir/trim/${id}.R2_val_2.fq.gz

prefix=$outdir/star/${id}
bamout=$outdir/star/${id}Aligned.sortedByCoord.out.bam

SLURM_CPUS_PER_TASK=10

sbatch << RUN
#!/usr/bin/bash
#SBATCH -J $id
#SBATCH --partition compute_new
#SBATCH -x node10
#SBATCH -c 10
#SBATCH --mem 10g
#SBATCH -o $outdir/logs/${id}.malbac.log
#SBATCH -e $outdir/logs/${id}.malbac.log
set -x

$trim_galore --fastqc --paired \
        --phred33 \
        --retain_unpaired \
        --output_dir $outdir/trim \
        $Read1 $Read2



$STAR --runThreadN $SLURM_CPUS_PER_TASK \
	--genomeDir $genome \
	--readFilesIn $trim_r1 $trim_r2 \
	--readFilesCommand zcat \
	--outFileNamePrefix $prefix \
	--outSAMtype BAM \
	SortedByCoordinate --quantMode TranscriptomeSAM



$samtools sort -@ $SLURM_CPUS_PER_TASK -n -o $outdir/${id}.sortn.bam $bamout
$samtools fixmate -@ $SLURM_CPUS_PER_TASK -m $outdir/${id}.sortn.bam $outdir/${id}.fixed.bam
$samtools sort -@ $SLURM_CPUS_PER_TASK -o $outdir/${id}.sorted.bam $outdir/${id}.fixed.bam
$samtools markdup -@ $SLURM_CPUS_PER_TASK -r  $outdir/${id}.sorted.bam $outdir/${id}.rmdup.bam
$samtools index $outdir/${id}.rmdup.bam



set +x
RUN

}

raw_dir=malbac_data
out_dir=output

mkdir -p $out_dir
mkdir -p $out_dir/malbac

ls $raw_dir|grep R1|sed 's/.R1.*//g'|while read id
do

	echo $id
	run_pip $raw_dir $id $out_dir/malbac


done

