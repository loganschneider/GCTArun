#!/bin/sh

#NOTE: run this script from within the directory containing the QC'ed files
#NOTE: this file must be run after QC.sh as it uses the trimmed_pruned files
echo "Have you run QC.sh? Respond Y or N, followed by [ENTER]: "
read QCdone
#NOTE: there must be an individual list with FID and IID columns (no column labels needed)
echo "Do you have an individual list with two columns (FID and IID) in this folder? Respond Y or N, followed by [ENTER]: "
read indlist

if [ $QCdone == "Y" ] && [ $indlist == "Y" ]
then
	echo "good to go"
else
	echo "You need to prepare the files for this script to run properly"
fi

module load plink
module load r

echo "Enter name of study population (e.g. WSC, MrOS, APOE), followed by [ENTER]: "
read study

echo "Enter number of genotype/chip/array pseudocohorts, followed by [ENTER]: "
read cohortnum

# gather genotype/chip/array pseudocohort(s) names into an array called "list"
for i in {1..$cohortnum}
do
	echo "Enter name(s) of genotype/chip/array pseudocohort(s), separated by spaces, followed by [ENTER]: "
	read cohortnames
	list=($cohortnames)
done

# determine phenotype name (as designated in the ${cohortname}_pheno.txt file)
#	can generate the pheno file, using the fam2pheno.R script: https://www.dropbox.com/s/g8e5zzwkvdc8ny0/fam2pheno.R?dl=0
echo "Enter BINARY phenotype name as it appears in the {cohortname}_pheno.txt file, followed by [ENTER]: "
read pheno

workDir=$PWD

mkdir -p $workDir/${study}_GCTA

# loop over cohort names in "list" to submit jobs
for k in "${list[@]}"
do

	# make directory to receive files for pseudocohort output
	mkdir -p $workDir/${k}_GCTA
	# copy in cleaned files (not trimmed or pruned)
	cp ${k}_clean2.* ${k}_GCTA
	# copy in the GCTA program
	cp -n /srv/gsfs0/projects/mignot/PLMGWAS/gcta64 .
	# copy pheno files into the respective pseudocohort GCTA directory
	cp ${k}_subset/${k}_pheno.txt ${k}_GCTA
	#generate study-of-interest only subset
	#-includes usage of phenotype from the ${k}_pheno.txt file with case=1,ctrl=0,else=-9/NA (but output is case=2,ctrl=1,else=0/-9/NA)
	plink --file $workDir/${k}_GCTA/${k}_clean2 --keep ${study}.indlist --allow-no-sex --pheno $workDir/${k}_GCTA/${k}_pheno.txt --pheno-name ${pheno} --1 --make-bed --out $workDir/${k}_GCTA/${k}_${study}only
	# run GCTA
	./gcta64 --bfile $workDir/${k}_GCTA/${k}_${study}only --autosome --maf 0.01 --make-grm --out $workDir/${k}_GCTA/${k}_${study}only --thread-num 10
	# this determines which column the input binary ${pheno} is in
	col=$(head -1 $workDir/${k}_GCTA/${k}_pheno.txt | tr -s ' ' '\n' | nl -nln |  grep "${pheno}" | cut -f1)
	phenocol=$((col - 2))
	# makes GCTA-compliant phenotype file - by removing header row
	sed '1d' $workDir/${k}_GCTA/${k}_pheno.txt > $workDir/${k}_GCTA/${k}.phen
	./gcta64 --grm $workDir/${k}_GCTA/${k}_${study}only --pheno $workDir/${k}_GCTA/${k}.phen --mpheno ${phenocol} --reml --out $workDir/${k}_GCTA/${k}_${study}only --thread-num 10
	
	# this will combine all pseudocohorts to run GCTA as well...larger N allows for better heritability estimate
	if [ -e $workDir/${study}_GCTA/${study}.fam ]
	then
		plink --bfile $workDir/${study}_GCTA/${study} --bmerge $workDir/${k}_GCTA/${k}_${study}only.bed $workDir/${k}_GCTA/${k}_${study}only.bim $workDir/${k}_GCTA/${k}_${study}only.fam --allow-no-sex --make-bed --out $workDir/${study}_GCTA/${study}
	else
		plink --bfile $workDir/${k}_GCTA/${k}_${study}only --allow-no-sex --make-bed --out $workDir/${study}_GCTA/${study}
	fi
	# generates GCTA-compliant phenotype file from the ${study}.fam file
	# also ensures case=1,ctrl=0,else=-9/NA
	sed -e 's/\<0\>/\<-9\>/g;s/NA/NA/g' $workDir/${study}_GCTA/${study}.fam | awk -v s=1 '{print $1, $2, $6-s}' > ${study}_GCTA/${study}_01.phen
	./gcta64 --bfile $workDir/${study}_GCTA/${study} --autosome --maf 0.01 --make-grm --out $workDir/${study}_GCTA/${study} --thread-num 10
	./gcta64 --grm $workDir/${study}_GCTA/${study} --pheno ${study}_GCTA/${study}_01.phen --reml --out $workDir/${study}_GCTA/${study} --thread-num 10
done

