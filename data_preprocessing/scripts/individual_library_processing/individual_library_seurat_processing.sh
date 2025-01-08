#!/usr/bin/bash

working_dir=$(pwd)
scripts_dir=...../scripts_for_sharing/data_preprocessing/scripts/individual_library_processing

macs2_path=/diskmnt/Projects/Users/rliu/Software/miniconda/envs/Seurat_RNA_ATAC_ST/bin/macs2
chr_size=/diskmnt/Projects/Users/rliu/Projects/normal_mouse_kidney_development/Resources/genome_reference/mm10.chrom.sizes.txt
out_path=${working_dir}/seurat_out/
LOG_DIR=${working_dir}/seurat_logs/


#for sample_id in 20210319-NMKE16_5F 20220130-E16_5F 20210218-NMKE16_5M 20220130-E16_5M 20210125-NMK0F 20220129-NMK0F 20210125-NMK0M 20220129-NMK0M 20210228-NMK3F 20210419-NMK3F GPAD-20210419-NMK3F2 20201201-NMK3M 20210419-NMK3M GPAD-20210419-NMK3M2 20210129-AKICT3F 20210423-NMK12F1 GPAD-20210423-NMK12F2 20210129-AKICT3M 20210423-NMK12M1 GPAD-20210423-NMK12M2 20201213-NMK52F GPAD-20210422-NMK52F1 20210422-NMK52F2 20201213-NMK52M GPAD-20210422-NMK52M1 20210422-NMK52M2 20210908-NMK92F1 20210908-NMK92F2 GPAD-20220124-NMK92F2 20210907-NMK92M1 20210907-NMK92M2 GPAD-20220124-NMK92M1;do
sample_id = 20210319-NMKE16_5F
scrublet_path=${working_dir}/scrublet_out/combined/${sample_id}/${sample_id}_combo_scrublet_output_table.csv
cellranger_path=/diskmnt/Projects/Mouse_kidney_development/combo_snRNA_snATAC/cellranger_arc_2.0/${sample_id}
Rscript ${scripts_dir}/src/seurat_pipeline_with_emptydroplet.R \
--sample=$sample_id --data=$cellranger_path \
--macs2_path=$macs2_path --output_folder=$out_path --chrom_size=$chr_size \
--prf_min=1000 --pct_min=15 --ns_max=5 --pc_first=2 --pc_num=50 --tss=3 \
--scrublet=$scrublet_path >> ${LOG_DIR}/${sample_id}.log 2>&1
#done
