#!/bin/bash

#################################################################
# Developer: Robert Lesurf
# Arguments: config file containing various directories and paths
# Delete the old VEP annotations and generate new ones
#################################################################


##### Read in the config file #####
while getopts ":c:" opt; do
  case $opt in
    c) CONFIG=$OPTARG;;
  esac
done
. ${CONFIG}


##### Move the old VEP annotations to a backup directory #####
CURRENTDATE=`date +"%Y-%m-%d-%H-%M-%S"`
mkdir ${RESULTS_DIR}/Annotate/Backup_${CURRENTDATE}
mv ${RESULTS_DIR}/Annotate/annotate_* ${RESULTS_DIR}/Annotate/Backup_${CURRENTDATE}
cd ${RESULTS_DIR}


##### Create per chr script #####
for contig in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do

prefix="chr"
contig_nc=${contig#"$prefix"}

cat << EOF > ${RESULTS_DIR}/Scripts/sbatch_reAnnotateVep_${contig}.sh
#!/bin/bash
#
#SBATCH -J ${SCRIPT_NAME}_${contig}_reAnnotateVep
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --tmp=50G
#SBATCH --mem=20g
#SBATCH --time=30:00:00
#SBATCH -e ${RESULTS_DIR}/Scripts/Logs/${SCRIPT_NAME}_${contig}_reAnnotateVep.err.txt
#SBATCH -o ${RESULTS_DIR}/Scripts/Logs/${SCRIPT_NAME}_${contig}_reAnnotateVep.out.txt


##### Load modules #####
module load bcftools/1.9
module load vep/102


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
chmod +x ${RESULTS_DIR}/Scripts/sbatch_reAnnotateVep_${contig}.sh
sbatch ${RESULTS_DIR}/Scripts/sbatch_reAnnotateVep_${contig}.sh
done

