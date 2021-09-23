

python3 pMHCpan_v2.py \
    --input-train train_v4_el_single_HLA_9AA_0.txt.gz \
    --input-validate train_v4_el_single_HLA_9AA_1.txt.gz \
    --hidden-size1 800 --hidden-size2 400 -L 1 \
    --olabel split0_9AA_w_pep_len_Aug24_wsun \
    -e 5 --n_iter 5 --save_validate_pred \
    > logfiles/pMHCpan_v2_800_split0_9AA_w_pep_len_Aug24_wsun.log &


python3 pMHCpan_v2.py \
    --input-train train_v4_el_single_HLA_9AA_0.txt.gz \
    --input-validate train_v4_el_single_HLA_9AA_1.txt.gz \
    --hidden-size1 800 --hidden-size2 400 -L 1 \
    --olabel split0_9AA_w_pep_len_8to11_Aug28_wsun \
    -e 5 --n_iter 5 --save_validate_pred \
    > logfiles/pMHCpan_v2_800_split0_9AA_w_pep_len__8to11_Aug28_wsun.log &

