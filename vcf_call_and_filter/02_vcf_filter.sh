name=“”
out_name=“”
gunzip Chr_1.$name.vcf.gz
cat Chr_1.$name.vcf > ${out_name}.vcf

for f in 2 3 4 5 6 7 8 9 10 11 12 13 14
do 
	gunzip Chr_$f.$name.vcf.gz
	grep -v "\#" Chr_$f.$name.vcf >> ${out_name}.vcf
done

bgzip ${out_name}.vcf
tabix  ${out_name}.vcf.gz



bcftools filter ER.combo.bcf.vcf -e 'QUAL<30 || INFO/RPB<0.01 || INFO/MQB<0.01 || INFO/BQB<0.01 ||  INFO/MQSB<0.01  || INFO/MQ<20' > ER.1.vcf

#RPB   1      Float   Mann-Whitney U test of Read Position Bias (bigger is better)                                        
#MQB   1      Float   Mann-Whitney U test of Mapping Quality Bias (bigger is better)                                      
#BQB   1      Float   Mann-Whitney U test of Base Quality Bias (bigger is better)                                         
#MQSB  1      Float   Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)         

bcftools view -h ER.1.vcf > ER.2.vcf
bedtools intersect -v -a ER.1.vcf -b /datacommons/willislab/repeats_tol.bed >> ER.2.vcf

bedtools intersect -a ER.2.vcf -b /datacommons/willislab/Mguttatus_genome_5/coding_sequences_tol.bed > ER.CDS.vcf
