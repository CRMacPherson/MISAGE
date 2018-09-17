#!/bin/bash
SNP=$1
OUTPUT_PATH='/Users/cmacpherson/Work/Milieu_Interieur/R_Projects/Shiny_Projects/Age_Sex_Genetics/newData'
#PLINK='/pasteur/projets/policy01/LabExMI/GWAS/LMM/FACS_GWAS/bin/plink1.9'
#BFILE='/pasteur/projets/policy01/LabExMI/GWAS/LMM/LabExMI_imputation_1000x5699237'
PLINK='plink'
BFILE='/Volumes/LabExMI/GWAS/LMM/LabExMI_imputation_1000x5699237'
PLINK_OUTPUT=$OUTPUT_PATH"/"$SNP

$PLINK --bfile $BFILE --snp $SNP --recodeAD --out $PLINK_OUTPUT #this gives the .raw files that you used before
cat $PLINK_OUTPUT.raw | cut -f  2,7 -d " " > $PLINK_OUTPUT.sel #this gives only donorID and SNP value
rm $PLINK_OUTPUT.raw
rm $PLINK_OUTPUT.log
