#!/bin/bash
#
#SBATCH --account=biodept
#SBATCH -p biodept,common                # Partition to submit to (comma separated)
#SBATCH --array=45-47,56-57,91-94,96-99,127-130,001
#SBATCH -J er_array              # Job name
#SBATCH -n 8                     # Number of cores
#SBATCH -N 1                     #Ensure that all cores are on one machine
#SBATCH --mem 10GB               # Memory in MB
#SBATCH -o er_%A_%a.out   # File for STDOUT (with jobid = %j)
#SBATCH -e er_%A_%a.err	 # File for STDERR (with jobid = %j)
#SBATCH --mail-type=END,FAIL          # Type of email notification: BEGIN,END,FAIL,A$
#SBATCH --mail-user=amelia.lawrence.nn@gmail.com  #Email where notifications will be sent




genome="/datacommons/willislab/Mguttatus_genome_5/MguttatusTOL_551_v5.0.fa"
path_to_picard="picard.jar"
module load samtools

#copy over TruSeq.fa from datacommons

#trim low quality reads
trimmomatic PE -threads 6 -phred33 -trimlog sample_${SLURM_ARRAY_TASK_ID}.trimlog -quiet -validatePairs sample_${SLURM_ARRAY_TASK_ID}.R1.fastq.gz sample_${SLURM_ARRAY_TASK_ID}.R2.fastq.gz sample_${SLURM_ARRAY_TASK_ID}.PE.R1.fq.gz sample_${SLURM_ARRAY_TASK_ID}.U.R1.fq.gz sample_${SLURM_ARRAY_TASK_ID}.PE.R2.fq.gz sample_${SLURM_ARRAY_TASK_ID}.U.R2.fq.gz ILLUMINACLIP:TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 AVGQUAL:25 
fastqc sample_${SLURM_ARRAY_TASK_ID}.PE.R1.fq.gz
fastqc sample_${SLURM_ARRAY_TASK_ID}.PE.R2.fq.gz

#Align reads
bwa mem -t 8 $genome sample_${SLURM_ARRAY_TASK_ID}.PE.R1.fq.gz  sample_${SLURM_ARRAY_TASK_ID}.PE.R2.fq.gz > sample_${SLURM_ARRAY_TASK_ID}.aligned.bam


#Sort reads
samtools sort sample_${SLURM_ARRAY_TASK_ID}.aligned.bam -T sample_${SLURM_ARRAY_TASK_ID}.sort -o sample_${SLURM_ARRAY_TASK_ID}.sort.bam
samtools index sample_${SLURM_ARRAY_TASK_ID}.sort.bam

#Build Index
java -Xmx2g -jar $path_to_picard BuildBamIndex INPUT=sample_${SLURM_ARRAY_TASK_ID}.sort.bam VALIDATION_STRINGENCY=LENIENT

#Fix Mate
java -Xmx14g -jar $path_to_picard FixMateInformation INPUT=sample_${SLURM_ARRAY_TASK_ID}.sort.bam OUTPUT=sample_${SLURM_ARRAY_TASK_ID}.FM.bam SORT_ORDER=coordinate TMP_DIR=$temp VALIDATION_STRINGENCY=LENIENT

#Mark Duplicates 
java -Xmx14g -jar $path_to_picard MarkDuplicates INPUT=sample_${SLURM_ARRAY_TASK_ID}.FM.bam OUTPUT=sample_${SLURM_ARRAY_TASK_ID}.MD.bam M=sample_${SLURM_ARRAY_TASK_ID}.metrics_file VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

#Read Groups 
java -Xmx14g -jar $path_to_picard AddOrReplaceReadGroups RGLB=sample_${SLURM_ARRAY_TASK_ID}.MD.bam RGPL=illumina RGPU=run RGSM=sample_${SLURM_ARRAY_TASK_ID} I=sample_${SLURM_ARRAY_TASK_ID}.MD.bam O=sample_${SLURM_ARRAY_TASK_ID}.RG.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT TMP_DIR=$temp
java -Xmx2g -jar $path_to_picard BuildBamIndex INPUT=sample_${SLURM_ARRAY_TASK_ID}.RG.bam VALIDATION_STRINGENCY=LENIENT

