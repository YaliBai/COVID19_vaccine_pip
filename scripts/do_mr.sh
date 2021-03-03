#!/usr/bin/bash


function do_mr {
## VDJ_SHMrate -- mutation rate
Dir='~/./ncbi-igblast-1.15.0
Rscript='~/./R/R-3.6.3/bin/Rscript'
sp=human

id=$1
query=$2   ##vdj_b/outs/filtered_contig.fasta
OUT=$3

out=$OUT/${id}.igblast_output.txt
reformat=$OUT/${id}.igblast_output_reformat.txt
mr=$OUT/${id}.mutation_rate.txt
mkdir -p $OUT/logs

sbatch << RUN
#!/usr/bin/bash
#SBATCH -J $id
#SBATCH --partition=compute_new
#SBATCH -o $OUT/logs/${id}.o.txt
#SBATCH -e $OUT/logs/${id}.e.txt
#SBATCH --mem=$mem 
#SBATCH -c 1

set -x

cd $Dir

cd $Dir

bin/igblastn \\
-germline_db_V database/${sp}_IG_V.fasta \\
-germline_db_D database/${sp}_IG_D.fasta \\
-germline_db_J database/${sp}_IG_J.fasta \\
-organism $sp \\
-query $query \\
-out $out \\
-auxiliary_data optional_file/${sp}_gl.aux \\
-show_translation -outfmt 7

grep ^[VDJ] $out |cut -f 1-6,8-10 > $reformat
$Rscript mutation_rate.R $reformat $mr

set +x
RUN

}



outdir=01-cellranger/PR10_XXL_10X-RBD_1/outs/vdj_b



do_mr RBD_1 ${outdir}/filtered_contig.fasta ${outdir}





