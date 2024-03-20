#!/bin/bash

#################################################################
# Developer: Robert Lesurf
# Arguments: config file containing various directories and paths
# Create and submit script to cluster for variant summarization
#################################################################


##### Read in the config file #####
while getopts ":c:" opt; do
  case $opt in
    c) CONFIG=$OPTARG;;
  esac
done
echo -e "$CONFIG"
. ${CONFIG}


##### Create script to qsub #####
cat << EOF > ${RESULTS_DIR}/Scripts/sbatch_Summarize.sh
#!/bin/bash
#
#SBATCH -J ${SCRIPT_NAME}_Summarize
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --tmp=50G
#SBATCH --mem=100G
#SBATCH --time=30:00:00
#SBATCH -e ${RESULTS_DIR}/Scripts/Logs/${SCRIPT_NAME}_Summarize.err.txt
#SBATCH -o ${RESULTS_DIR}/Scripts/Logs/${SCRIPT_NAME}_Summarize.out.txt


##### Load modules #####
module load R/4.1.2
module load bcftools/1.9
module load gatk/4.2.2.0


##### Launch the Rscript #####
Rscript --vanilla $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/summarize_annotations.R ${CONFIG}


##### End the script #####
EOF

##### Submit script to the cluster #####
chmod +x ${RESULTS_DIR}/Scripts/sbatch_Summarize.sh
sbatch ${RESULTS_DIR}/Scripts/sbatch_Summarize.sh

