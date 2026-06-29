#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --partition io,cpu
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL,COMPLETED
#SBATCH --cpus-per-task 8
#SBATCH -t 50:00:00
#SBATCH --mem 20G
 
donor=$1
mybatch=$SLURM_ARRAY_TASK_ID
 
 
# Folder to upload (pre-encrypted .gpg files should live here)
UPLOAD_DIR="/gpfs/commons/groups/landau_lab/ResolveOME/EGA_upload/${donor}_batch${mybatch}/"
mkdir $UPLOAD_DIR

mylist=$(awk '{print $2}' data_paths/${donor}_batches/bam.list.${donor}.batch${mybatch}.txt | paste -sd, -)
echo $mylist

echo "Starting encryptation for donor "${donor}" batch "$mybatch
module load Java/21.0.2
java -Xmx16g -jar EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar \
	-f \
	-i "$(awk '{print $2}' data_paths/${donor}_batches/bam.list.${donor}.batch${mybatch}.txt | paste -sd, -)" \
	-o $UPLOAD_DIR

ls ${UPLOAD_DIR}*ba* > ${UPLOAD_DIR}/ListFiles.txt

#source ~/.ega_credentials
#module purge
#module load aspera/4.19.0-GCCcore-12.3.0
#export PATH="$HOME/.aspera/sdk:$PATH"   # adjust to the dir from `find` above
#which ascp
#ascp --version
#
#export ASPERA_SCP_PASS=$ASPERA_SCP_PASS
#ascp -P33001 -O33001 -QT -l300M \
#	-L- \
#	--mode=send \
#    --user="${EGA_BOX}" \
#    --host=fasp.ega.ebi.ac.uk \
#    --file-list="${UPLOAD_DIR}/ListFiles.txt" \
#    /

