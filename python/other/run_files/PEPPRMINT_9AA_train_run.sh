
python3 mixPep_v2.py \
    --input_train train_v4_el_multi_HLA_9AA_0.txt.gz \
    --input_valid train_v4_el_multi_HLA_9AA_1.txt.gz \
    --input_test ../pMHCpan_data/test_v4_el_single_HLA_9AA.txt.gz \
    --init_model pMHCpan_800_bs_32_lr_0.001_e_5_layer_1_split0_9AA_w_pep_len_8to11_Aug28_wsun_best_valid.h5 \
    --hidden_size1 800 --hidden_size2 400 -L 1 -e 3 --n_iter 10 \
    --olabel 9AA_w_pep_len_8to11_wsun_Aug28 \
    --save_all_iterations --save_all_pred --save_model \
    > logfiles/mixpep_800_layer_1_split0_9AA_w_pep_len_8to11_Aug28_wsun.log &
    



python3 mixPep_v2.py \
    --input_train train_v4_el_multi_HLA_9AA_0.txt.gz \
    --input_valid train_v4_el_multi_HLA_9AA_1.txt.gz \
    --input_test ../pMHCpan_data/test_v4_el_single_HLA_9AA.txt.gz \
    --init_model pMHCpan_800_bs_32_lr_0.001_e_5_layer_1_split0_9AA_w_pep_len_Aug24_wsun_best_valid.h5 \
    --hidden_size1 800 --hidden_size2 400 -L 1 -e 3 --n_iter 10 \
    --olabel 9AA_w_pep_len_Aug24_wsun \
    --save_all_iterations --save_all_pred --save_model \
    > logfiles/mixpep_800_layer_1_split0_9AA_w_pep_len_Aug24_wsun.log 


python3 mixPep_v2.py \
    --input_train train_v4_el_multi_HLA_balanced_0.txt.gz \
    --input_valid train_v4_el_multi_HLA_balanced_1.txt.gz \
    --input_test ../pMHCpan_data/test_v4_el_single_HLA.txt.gz \
    --init_model pMHCpan_800_bs_32_lr_0.001_e_50_layer_1_split0_one_input_Jun29_wsun_best_valid.h5 \
    --hidden_size1 800 --hidden_size2 400 -L 1 -e 2 \
    --olabel pMHC_balanced_MA_balanced_0801_wsun \
    --save_all_iterations --save_all_pred --save_model \
    > logfiles/mixpep_800_layer_1_split0_pMHC_balanced_MA_balanced_0801_wsun.log 
