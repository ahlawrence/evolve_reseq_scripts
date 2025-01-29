ls *.RG.bam > bamlist.txt

module load bcftools
genome="/hpc/home/ahl28/Mguttatus_genome_5/MguttatusTOL_551_v5.0.fa"
for f in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do 
	echo '#!/bin/bash' > chr_$f.sh
  	echo '#$ -S /bin/bash' >> chr_$f.sh 
  	echo '#SBATCH --get-user-env' >> chr_$f.sh
  	echo '#SBATCH --job-name='$c >> chr_$f.sh
  	echo '#SBATCH --output'=$f.out >> chr_$f.sh
  	echo '#SBATCH --error='$f.err >> chr_$f.sh
  	echo '#SBATCH --cpus-per-task=1' >> chr_$f.sh
  	echo '#SBATCH --account=biodept' >> chr_$f.sh
  	echo '#SBATCH -p common,biodept' >> chr_$f.sh
  	echo '#SBATCH --mail-type=END' >> chr_$f.sh
  	echo '#SBATCH --mail-user=amelia.lawrence.nn@gmail.com' >> chr_$f.sh
  	echo '#SBATCH --mem=30G' >> chr_$f.sh
  	echo -en '\n' >> chr_$f.sh
	echo bcftools mpileup --redo-BAQ --min-BQ 30 --per-sample-mF --annotate FORMAT/AD,FORMAT/DP,INFO/AD -f $genome -b bamlist.txt -r Chr_$f '|' bcftools call --multiallelic-caller -Oz '>' all_indel.$f.bcftools.vcf.gz >> chr_$f.sh
done

for file in chr_*.sh; do sbatch $file; done

