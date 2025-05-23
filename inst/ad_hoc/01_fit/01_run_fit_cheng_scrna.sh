#!/bin/bash

START=$(date +%s)

source ${rattaca_args}

echo ""
echo "------------------------------------------"
echo "----------- RATTACA MODEL FITS -----------"
echo "------------------------------------------"
echo ""

((line=SLURM_ARRAY_TASK_ID))
trait=$(head -n ${line} ${trait_list} | tail -n 1)

echo "$(date +"%Y%m%d-%H:%M:%S")"
echo ""
echo "RATTACA generation: ${generation}"
echo "Fitting model for trait: ${trait}"
echo "Using script: ${fit_Rscript}"
echo ""
echo "Phenotype data: ${train_pheno_file}"
echo "Genotype data: ${train_geno_input}"
echo "Test genotypes: ${test_geno_input}" 


# run predictions on the trait
source activate rattaca

# cd ${fit_dir}
# fit_Rscript=$(basename ${fit_Rscript})

if [ ${sample_type} == "ldprune" ]; then 

    Rscript ${fit_Rscript} \
        --pred_type ${pred_type} \
        --sample_type ${sample_type} \
        --generation ${generation} \
        --trait ${trait} \
        --trait_dict ${trait_dict} \
        --rfid_list ${rfid_list} \
        --train_pheno_file ${train_pheno_file} \
        --train_geno_input ${train_geno_input} \
        --test_geno_input ${test_geno_input} \
        --proj_dir ${proj_dir} \
        --bpar_dir ${bpar_dir} \
        --plink1 ${plink1} \
        --plink2 ${plink2} \
        --maf_cutoff ${maf_cutoff} \
        --missing_cutoff ${missing_cutoff} \
        --hwe_cutoff ${hwe_cutoff} \
        --ld_r2 ${ld_r2} \
        --ld_window_size ${ld_window_size} \
        --ld_step_size ${ld_step_size} \
        rattaca

elif [ ${sample_type} == "random" ]; then 

    Rscript ${fit_Rscript} \
        --pred_type ${pred_type} \
        --sample_type ${sample_type} \
        --generation ${generation} \
        --trait ${trait} \
        --trait_dict ${trait_dict} \
        --rfid_list ${rfid_list} \
        --train_pheno_file ${train_pheno_file} \
        --train_geno_input ${train_geno_input} \
        --test_geno_input ${test_geno_input} \
        --proj_dir ${proj_dir} \
        --bpar_dir ${bpar_dir} \
        --plink1 ${plink1} \
        --plink2 ${plink2} \
        --maf_cutoff ${maf_cutoff} \
        --missing_cutoff ${missing_cutoff} \
        --hwe_cutoff ${hwe_cutoff} \
        --n_samples ${n_samples} \
        --n_snps ${n_snps} \
        rattaca


fi
conda deactivate


END=$(date +%s)
echo "Time elapsed for ${trait} predictions: $(( $END - $START )) seconds"
