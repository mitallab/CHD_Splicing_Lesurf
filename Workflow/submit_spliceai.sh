#!/bin/bash

#################################################################
# Developer: Robert Lesurf
# Arguments: config file containing various directories and paths
# Run SpliceAI, various annotations, and filters
#################################################################


##### Read in the config file #####
while getopts ":c:" opt; do
  case $opt in
    c) CONFIG=$OPTARG;;
  esac
done
. ${CONFIG}


##### Setup #####
mkdir -p ${RESULTS_DIR}
cp ${CONFIG} ${RESULTS_DIR}/
mkdir ${RESULTS_DIR}/Scratch
mkdir ${RESULTS_DIR}/Scripts
mkdir ${RESULTS_DIR}/Scripts/Logs
mkdir ${RESULTS_DIR}/Norm
mkdir ${RESULTS_DIR}/SpliceAI
mkdir ${RESULTS_DIR}/Annotate
mkdir ${RESULTS_DIR}/Summary
cd ${RESULTS_DIR}


##### Create per chr script #####
for contig in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do

prefix="chr"
contig_nc=${contig#"$prefix"}

cat << EOF > ${RESULTS_DIR}/Scripts/sbatch_${contig}.sh
#!/bin/bash
#
#SBATCH -J ${SCRIPT_NAME}_${contig}
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --tmp=70G
#SBATCH --mem=30g
#SBATCH --time=100:00:00
#SBATCH -e ${RESULTS_DIR}/Scripts/Logs/${SCRIPT_NAME}_${contig}.err.txt
#SBATCH -o ${RESULTS_DIR}/Scripts/Logs/${SCRIPT_NAME}_${contig}.out.txt


##### Load modules #####
module load bcftools/1.9
module load python/3.7.6
module load vep/102


##### Split the variants #####
bcftools view -f PASS -R ${GENES_BED_BASEPATH}${contig}.bed ${INPUT_VCF} | \
bcftools norm -m - | \
bcftools annotate -x ID,INFO,^FORMAT/GT,FORMAT/DP,FORMAT/DPI -O z -o ${RESULTS_DIR}/Norm/norm_${contig}.vcf.gz
bcftools index -t ${RESULTS_DIR}/Norm/norm_${contig}.vcf.gz


############### SpliceAI: Use pre-computed SpliceAI ###############
##### Split out SpliceAI-precomputed SNVs/indels, and non-precomputed SNVs/indels #####
mkdir ${RESULTS_DIR}/SpliceAI/${contig}

bcftools isec -O z -o ${RESULTS_DIR}/SpliceAI/${contig}/intersect.vcf.gz -n=2 -w1 ${RESULTS_DIR}/Norm/norm_${contig}.vcf.gz ${SPLICEAI_PRECOMPUTED_ALL_BASEPATH}${contig}.vcf.gz
bcftools index -t ${RESULTS_DIR}/SpliceAI/${contig}/intersect.vcf.gz
bcftools isec -O z -o ${RESULTS_DIR}/SpliceAI/${contig}/setdiff.vcf.gz -n~10 -w1 ${RESULTS_DIR}/Norm/norm_${contig}.vcf.gz ${SPLICEAI_PRECOMPUTED_ALL_BASEPATH}${contig}.vcf.gz
bcftools index -t ${RESULTS_DIR}/SpliceAI/${contig}/setdiff.vcf.gz

##### Annotate with precomputed SpliceAI calls, and run new SpliceAI calls #####
bcftools annotate -a ${SPLICEAI_PRECOMPUTED_ALL_BASEPATH}${contig}.vcf.gz -c INFO/SpliceAI -O z -o ${RESULTS_DIR}/SpliceAI/${contig}/precomputed.vcf.gz ${RESULTS_DIR}/SpliceAI/${contig}/intersect.vcf.gz
bcftools index -t ${RESULTS_DIR}/SpliceAI/${contig}/precomputed.vcf.gz
spliceai -M 1 -D 100 -R ${REFERENCE_FASTA} -A grch38 -I ${RESULTS_DIR}/SpliceAI/${contig}/setdiff.vcf.gz -O ${RESULTS_DIR}/SpliceAI/${contig}/postcomputed.vcf
bgzip ${RESULTS_DIR}/SpliceAI/${contig}/postcomputed.vcf
bcftools index -t ${RESULTS_DIR}/SpliceAI/${contig}/postcomputed.vcf.gz

##### Concatenate the SpliceAI outputs #####
bcftools concat ${RESULTS_DIR}/SpliceAI/${contig}/precomputed.vcf.gz ${RESULTS_DIR}/SpliceAI/${contig}/postcomputed.vcf.gz | \
bcftools sort --temp-dir ${RESULTS_DIR}/SpliceAI/${contig} | \
bcftools norm -m - -O z -o ${RESULTS_DIR}/SpliceAI/spliceai_${contig}.vcf.gz
bcftools index -t ${RESULTS_DIR}/SpliceAI/spliceai_${contig}.vcf.gz

############### End SpliceAI Annotation ###############


##### Filter SpliceAI to keep only variants with Delta score >= 0.2 #####
filter_vep \
-i ${RESULTS_DIR}/SpliceAI/spliceai_${contig}.vcf.gz -o ${RESULTS_DIR}/SpliceAI/spliceai_${contig}_DS02.vcf \
--format vcf --vcf_info_field SpliceAI --filter 'DS_AG >= 0.2 or DS_AL >= 0.2 or DS_DG >= 0.2 or DS_DL >= 0.2'
bgzip ${RESULTS_DIR}/SpliceAI/spliceai_${contig}_DS02.vcf
bcftools index -t ${RESULTS_DIR}/SpliceAI/spliceai_${contig}_DS02.vcf.gz

##### Annotate with lots of other variant information #####
vep \
-i ${RESULTS_DIR}/SpliceAI/spliceai_${contig}_DS02.vcf.gz -o ${RESULTS_DIR}/Annotate/annotate_${contig}_DS02_vep.vcf.gz -fasta ${REFERENCE_FASTA} --assembly GRCh38 \
--fork 4 --offline --cache --dir_cache /hpf/tools/centos6/vep/cache102 --format vcf --vcf --compress_output bgzip \
--variant_class --hgvs --protein --symbol --ccds --numbers --canonical --mane --biotype --af --af_1kg --af_esp \
--custom ${CLINVAR_VCF},ClinVar,vcf,exact,0,ALLELEID,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNSIG,CLNSIGCONF \
--custom ${HGMD_VCF},HGMD,vcf,exact,0,CLASS,PHEN,RANKSCORE \
--custom ${GNOMAD2E_DIR}/gnomad.exomes.r2.1.1.sites.${contig_nc}.liftover_grch38.vcf.bgz,gnomAD2e,vcf,exact,0,AC,AN,AF,popmax,faf95,AC_popmax,AN_popmax,AF_popmax,AF_afr,AF_asj,AF_nfe,AF_amr,AF_fin,AF_sas,AF_eas,AF_oth,AF_female,AF_male,lcr,nonpar,segdup \
--custom ${GNOMAD2G_DIR}/gnomad.genomes.r2.1.1.sites.${contig_nc}.liftover_grch38.vcf.bgz,gnomAD2g,vcf,exact,0,AC,AN,AF,popmax,faf95,AC_popmax,AN_popmax,AF_popmax,AF_afr,AF_asj,AF_nfe,AF_amr,AF_fin,AF_sas,AF_eas,AF_oth,AF_female,AF_male,lcr,nonpar,segdup \
--custom ${GNOMAD3_DIR}/gnomad.genomes.v3.1.2.sites.${contig}.vcf.bgz,gnomAD3,vcf,exact,0,AC,AN,AF,popmax,faf95_popmax,AC_popmax,AN_popmax,AF_popmax,AF_afr,AF_ami,AF_asj,AF_mid,AF_nfe,AF_amr,AF_fin,AF_sas,AF_eas,AF_oth,AF_XX,AF_XY,lcr,nonpar,segdup

bcftools index -t ${RESULTS_DIR}/Annotate/annotate_${contig}_DS02_vep.vcf.gz

##### end per chr execution and submit to cluster #####
EOF
chmod +x ${RESULTS_DIR}/Scripts/sbatch_${contig}.sh
sbatch ${RESULTS_DIR}/Scripts/sbatch_${contig}.sh
done

