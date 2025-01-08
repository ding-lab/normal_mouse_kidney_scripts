#!/usr/bin/bash

# using scrublet for doublet detection, this includes the doublet detection for RNA and ATAC portion. The cells that are defined as doublets from either RNA or ATAC modality would be defined as doublets
working_dir=$(pwd)
scripts_dir=...../scripts_for_sharing/data_preprocessing/scripts/individual_library_processing

sample_id=20210319-NMKE16_5F
input_dir=/diskmnt/Projects/Mouse_kidney_development/combo_snRNA_snATAC/cellranger_arc_2.0/${sample_id} ## input would be a directory to a Cellranger-ARC output. I am using one library as an example
input_library=$(basename "$input_dir")

cd ${working_dir}
output_dir=${working_dir}/scrublet_out/${input_library}
mkdir -p ${output_dir}
log_dir=${working_dir}/scrublet_logs/${input_library}
mkdir -p ${log_dir}

cd ${output_dir}
mkdir RNA ATAC combined
cd RNA;${scripts_dir}/src/scrublet-RNA-auto.sh ${input_dir} 1>${log_dir}/${input_library}_RNA.log 2>${log_dir}/${input_library}_RNA.err;cd ..
cd ATAC;${scripts_dir}/src/scrublet-ATAC-auto.sh ${input_dir} 1>${log_dir}/${input_library}_ATAC.log 2>${log_dir}/${input_library}_ATAC.err;cd ..
cd combined;bash ${working_dir}/src/scrublet-auto-combining.sh ${input_dir} ${output_dir}/RNA ${output_dir}/ATAC 1>${log_dir}/${input_library}_combined.log 2>${log_dir}/${input_library}_combined.err;cd ..

