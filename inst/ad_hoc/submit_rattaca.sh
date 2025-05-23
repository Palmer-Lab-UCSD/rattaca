#!/bin/bash

# use this script to run the RATTACA genomic predictions pipeline for a given project
# usage: bash submit_rattaca.sh 

#### MANUAL INPUT ####
rattaca_args=/path/to/rattaca_args.cfg
####

# read in arguments
source ${rattaca_args}


#### step 1: fit ####

job_name=rattaca_fit_${generation}

# create log output directory
logs_dir=${log_dir}/01_fit
mkdir -p ${logs_dir}

# run job array based on the number of traits to predict
n_traits=$(cat ${trait_list} | wc -l)

# submit the job
sbatch -J ${job_name} -o ${logs_dir}/${job_name}-%A-%a.o -e ${logs_dir}/${job_name}-%A-%a.o \
	-p condo -q condo -N 1 -c ${cpu_fit} --mem-per-cpu ${mem_fit} \
	--array 1-${n_traits} --time 4:00:00 --account ${allocation} --mail-type END,FAIL \
	--export rattaca_args="${rattaca_args}",generation="${generation}",test_geno_input="${test_geno_input}" \
	${fit_script}


#### step 2: predict ####

job_name=rattaca_pred_${generation}

# create log output directory
logs_dir=${log_dir}/02_predict
mkdir -p ${logs_dir}

# run job array based on the number of traits to predict
n_traits=$(cat ${trait_list} | wc -l)

# submit the job
sbatch -J ${job_name} -o ${logs_dir}/${job_name}-%A-%a.o -e ${logs_dir}/${job_name}-%A-%a.o \
	-p condo -q condo -N 1 -c ${cpu_pred} --mem-per-cpu ${mem_pred} \
	--array 1-${n_traits} --time 1:00:00 --account ${allocation} --mail-type END,FAIL \
	--export rattaca_args="${rattaca_args}",generation="${generation}",test_geno_input="${test_geno_input}" \
	${pred_script}


#### step 3: summarize ####

source activate ${conda_env}

# create log output directory
logs_dir=${log_dir}/03_summary
mkdir -p ${logs_dir}
log=${logs_dir}/rattaca_${generation}_summary
touch ${log}

results_dir=${proj_dir}/results
preds_file=${results_dir}/rattaca_${generation}_merged_predictions.csv
	
# merge predictions and summarize results
Rscript ${summary_Rscript} \
	--trait_list ${trait_list} \
	--trait_dict ${trait_dict} \
	--results_dir ${results_dir} \
	--generation ${generation} \
	rattaca > ${log}

if [ ! -z ${composite_trait_list+x} ]; then

	Rscript ${composite_trait_Rscript} \
		--trait_list ${trait_list} \
		--trait_dict ${trait_dict} \
		--predictions ${preds_file} \
		--composite_trait_dir ${composite_trait_dir} \
		--composite_trait_list ${composite_trait_list} \
		--results_dir ${results_dir} \
		--generation ${generation} \
		rattaca >> ${log}
fi

done

conda deactivate