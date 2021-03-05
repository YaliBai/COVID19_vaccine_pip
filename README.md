# COVID19_vaccine_pip


---------------------------------------
### > MALBAC-DT RNA-sequencing data processing
> scripts/do_star.sh




### > Call donors' SNPs and 10x cell assignment
>scripts/Snp_phased.v03.sh




#### >> Get ref SNPs
total SNPs:15,160,915
>>not in repeats
>> bedtools intersect -a hg38/dbsnp/dbSnp153Common.bed -b hg38/repeatMask/hg38.repeats.bed -v > SNPs/hg38_snp153_repeatmasker.bed
6942367 SNPs left

>>not in imprinted genes (from geneimprint)
>>
bedtools intersect -a SNPs/hg38_snp153_repeatmasker.bed -b hg38/imprinted_Genes/Impri.pos.bed -v > SNPs/hg38_snp153_rmi.bed
6909945 SNPs left





### > the Somatic hypermutation rate
>> main software: ncbi-igblast-1.15.0

> scripts/do_mr.sh





