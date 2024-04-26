
# conda create --name tf python=3.9
# conda activate tf
# conda install tensorflow

# input_peptides  = "/Users/wsun/research/data/TCR/TCGA/_results/step8a_TCGA_peptides_mut.tsv"
# input_alist     = "/Users/wsun/research/data/TCR/TCGA/_results/step8a_TCGA_hla_list.tsv"
# input_mhc_seq   = "/Users/wsun/research/GitHub/PEPPRMINT/data/NetMHCpan4_1_train/MHC_pseudo.dat"
# data_label      = "TCGA"
# model           = "/Users/wsun/research/GitHub/PEPPRMINT/results/PEPPRMINT/MA_800_split4.h5"
# model_tag       = "800_split4"
# results_dir     = "/Users/wsun/research/data/TCR/TCGA/_results/PEPPRMINT"
# save_all_pred   = True
# input_pep_hla   = False
# encode9AA       = False
# neoantigen      = True


#!/bin/bash

# Directory containing model files
model_dir="/Users/wsun/research/GitHub/PEPPRMINT/results/PEPPRMINT"
output_dir="/Users/wsun/research/data/TCR/TCGA/_results"

# Loop over each file in the model directory
for model_file in ${model_dir}/*.h5; do
    echo "${model_file}"

    # Only process files (skip directories)
    if [ -f "$model_file" ]; then
        # Extract model_tag
        model_file=$(basename "$model_file")

        model_tag=$(echo "$model_file" | cut -d'.' -f1 | cut -d'_' -f2,3)
        
        echo "${model_tag}"

        # Run the Python script with model file and output file as arguments
        python ~/research/GitHub/PEPPRMINT/python/_prediction_PEPPRMINT.py \
        --input_peptides "${output_dir}/step8a_TCGA_peptides_mut.tsv" \
        --input_alist "${output_dir}/step8a_TCGA_hla_list.tsv" \
        --input_mhc_seq /Users/wsun/research/GitHub/PEPPRMINT/data/NetMHCpan4_1_train/MHC_pseudo.dat \
        --data_label TCGA \
        --model "${model_dir}/${model_file}" \
        --model_tag "${model_tag}" \
        --results_dir "${output_dir}/PEPPRMINT" \
        --neoantigen --save_all_pred \
        > "${output_dir}/PEPPRMINT/TCGA_MA_${model_tag}.log"
    fi
done

