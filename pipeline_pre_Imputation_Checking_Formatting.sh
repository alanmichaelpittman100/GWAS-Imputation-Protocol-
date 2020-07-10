#!/bin/bash

#V1
#Alan Pittman 2018

###check the input plink file prior to imputation
###amend and/or discard any inconsisten alleles/ref-alleles/strand issues
###convert each plink file to a .vcf file ready for Impuation on the Michigan Imputation Server
###autosomes 
###X-chromosome separate

file=$1 #specify your input file here  

#--------------------------------------------------------------------------------------------
#resources

plink=/resourses/plink_linux_x86_64/plink
SHAPEIT=/resourses/shapeit.v2.904/bin/shapeit
checkVCF=/imputation/checkVCF.py

#Dirs

inPLINK=/imputation/starting_plink_files
freqOUT=/imputation/frequency_outputs
updatedOUT=/imputation/updated_plink_files
PhasedOUT=/imputation/shapeit_phase_out
VCF_out=/imputation/pre-phased_vcfs

#--------------------------------------------------------------------------------------------

chrNo="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

#--------------------------------------------------------------------------------------------

echo " "


echo " "
echo \
'''---------------------------------------------------------------------------------------------

AUTOMATED IMPUTATION CHECKING AND FORMATING

---------------------------------------------------------------------------------------------'''
date
echo " "
sleep 1

sleep 1
echo " "
echo \
'''---------------------------------------------------------------------------------------------

CALCULATING ALLELE FREQUENCIES

---------------------------------------------------------------------------------------------'''
echo " "
sleep 1

$plink --bfile $inPLINK/${file} --freq --out $freqOUT/$file

cp $inPLINK/${file}.* ./

sleep 1
echo " "
echo \
'''---------------------------------------------------------------------------------------------

EXECUTING HRC IMPUTATION CHECKING TOOL

---------------------------------------------------------------------------------------------'''
echo " "
sleep 1

###Execute checking script:
perl /imputation/HRC-1000G-check-bim-v4.2_AP_MODIFIED.pl \
 -b $inPLINK/${file}.bim \
 -f $freqOUT/${file}.frq \
 -r /resourses/HRC_reference/1000GP_Phase3_combined.legend \
 -p EUR -g
 
sleep 1
echo " "
echo \
'''---------------------------------------------------------------------------------------------

EXECUTING PLINK COMMANDS TO REPAIR ALLELE CONFIGURATION

---------------------------------------------------------------------------------------------'''
echo " "
sleep 1

##run this self -generating script running plink commands to fix any snps/alleles/positions/names etc.
sh Run-plink.sh

#------------------------------------------------------------------------------------------------
#housekeeping

mv *.txt log_files/
mv *.log log_files/
mv /homes/athosnew/Alan/Alan_Bioinformatics_Projects/imputation/${file}-updated* updated_plink_files/

rm $file.fam $file.bed $file.bim
 
#------------------------------------------------------------------------------------------------

sleep 1
echo " "
echo \
'''---------------------------------------------------------------------------------------------

PRE PHASING HAPLOTYPES FOR IMPUTATION

---------------------------------------------------------------------------------------------'''
echo " "
sleep 1

for chr in $chrNo; do

	$SHAPEIT --input-bed $updatedOUT/${file}-updated-chr$chr.bed $updatedOUT/${file}-updated-chr$chr.bim $updatedOUT/${file}-updated-chr$chr.fam \
	-M /homes/athosnew/Genetics_Centre_Bioinformatics/resourses/shapeit.v2.904/bin/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt \
	--force \
	-O $PhasedOUT/${file}_${chr}.phased \
	-T 10

done

sleep 3
echo " "
echo \
'''---------------------------------------------------------------------------------------------

GENERATE VCF FILES

---------------------------------------------------------------------------------------------'''
echo " "
sleep 1

for chr in $chrNo; do

	$SHAPEIT -convert \
	--input-haps $PhasedOUT/${file}_${chr}.phased \
	--output-vcf $VCF_out/${file}_${chr}.phased.vcf
		
done

sleep 1
echo " "
echo \
'''---------------------------------------------------------------------------------------------

CHECK OUR FINAL VCF FOR INTEGRITY

---------------------------------------------------------------------------------------------'''
echo " "
sleep 1


for chr in $chrNo; do

	python checkVCF.py \
	-r /homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/human_g1k_v37.fasta \
	-o $VCF_out/${file}_${chr}.phased.vcf.ref.check $VCF_out/${file}_${chr}.phased.vcf
		
done

sleep 1
echo " "
echo \
'''---------------------------------------------------------------------------------------------

ZIPPING AND INDEXING READY FOR IMPUTATION SERVER UPLOAD

---------------------------------------------------------------------------------------------'''
echo " "
sleep 1

for chr in $chrNo; do

/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/tabix-0.2.6/bgzip $VCF_out/${file}_${chr}.phased.vcf

/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/tabix-0.2.6/tabix $VCF_out/${file}_${chr}.phased.vcf.gz

done


#final folder housekeeping

mv *.log log_files/
mv *.mm log_files/



sleep 1
echo " "
echo \
'''---------------------------------------------------------------------------------------------

YOUR FILES ARE NOW READY! HAPPY GWASSING !!!!

---------------------------------------------------------------------------------------------'''
echo " "
sleep 1

exit
exit
exit



#now for the X chromosome


















