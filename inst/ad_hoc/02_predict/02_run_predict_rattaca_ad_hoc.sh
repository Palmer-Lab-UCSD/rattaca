#!/bin/bash

START=$(date +%s)
source ${rattaca_args}
source activate ${conda_env}

echo ""
echo "-------------------------------------------"
echo "-------- RATTACA TRAIT PREDICTIONS --------"
echo "-------------------------------------------"
echo ""
echo ""
echo "$(date +%Y%m%d-%H:%M:%S)"
echo "user: ${user}"
echo ""
echo ""

((line=SLURM_ARRAY_TASK_ID))
trait=$(head -n ${line} ${trait_list} | tail -n 1)
bpar_file=${proj_dir}/results/${trait}/${trait}.bpar

echo "RATTACA generation: ${generation}"
echo "Predicting trait: ${trait}"
echo "Using script: ${pred_Rscript}"
echo "Using bpar file: ${bpar_file}"
echo ""


# run predictions on the trait
Rscript ${pred_Rscript} \
    --pred_type ${pred_type} \
    --trait ${trait} \
    --rfid_list ${rfid_list} \
    --train_pheno_file ${train_pheno_file} \
    --train_geno_input ${train_geno_input} \
    --test_geno_input ${test_geno_input} \
    --proj_dir ${proj_dir} \
    --bpar_file ${bpar_file} \
    --generation ${generation} \
    --plink1 ${plink1} \
    --plink2 ${plink2} \
    rattaca


conda deactivate


END=$(date +%s)
echo ""
echo "Time elapsed for ${trait} predictions: $(( $END - $START )) seconds"
