#!/bin/bash/python
#python3 combining-sn-and-atac-scrublet.py --sample insert_sample_name_here --rna_input /diskmnt/Projects/PDX_scRNA_analysis/matching_snRNAseq/PDAC/scrubletv3-RNA --atac_input /diskmnt/Projects/PDX_scRNA_analysis/matching_snRNAseq/PDAC/scrubletv3-ATAC

import sys
import pandas as pd

sample = sys.argv[2]
rna_input_dir = sys.argv[4]
atac_input_dir = sys.argv[6]
input_rna_path = rna_input_dir +"/"+ sample +"/"+sample+"_scrublet_output_table.csv"
input_atac_path = atac_input_dir +"/"+ sample +"/"+sample+"_scrublet_output_table.csv"

rna_in = pd.read_csv(input_rna_path, sep=",", header = 0)
atac_in = pd.read_csv(input_atac_path, sep=",", header = 0)

combo_out=rna_in.merge(atac_in, how="inner", on='Barcodes')
print(combo_out[0:5])
combo_out.columns = ['Barcodes','doublet_score_rna','predicted_doublet_rna','doublet_score_atac','predicted_doublet_atac']

combo_out_path = sample +"_combo_scrublet_output_table.csv"
combo_out.to_csv(combo_out_path, sep = ",", index=False)

