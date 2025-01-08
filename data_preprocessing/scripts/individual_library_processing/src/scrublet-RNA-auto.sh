
#scrublet-RNA-auto.sh
#must be run in conda environment containg python3.9, pandas, numpy, scipy.io, matplotlib, random, sklearn.cluster
#Example command call: bash scrublet-RNA-auto.sh /diskmnt/Projects/PDX_scRNA_analysis/matching_snRNAseq/PDAC/cellranger-arc2
#file structure for when command is run:
#incoming-path/
            # /cellranger #provide the fullpath here
            #            /Sample-1
            #                     /outs
            #            /sample-2
            #                     /outs
            #            /3rd-sample
            #                       /outs
            # /scrublet #this directory can actually be anywhere you want it to be so long as it exists. It also can be named something else. provide the full path here second
            #            /Sample-1
            #                     /Sample-1_scrublet_output_table.csv
            #            /sample-2
            #                     /sample-2_scrublet_output_table.csv
            #            /3rd-sample
            #                       /3rd-sample_scrublet_output_table.csv
#
folder=$(pwd)
cellrangerpath=$1
#scrublet=$2
for dir in $cellrangerpath/*; do
    dir=${dir%*/}
    dir="${dir##*/}"
    echo $dir
    if [ -d "$cellrangerpath"/"$dir" ]
    then
        if [ -d "$cellrangerpath"/"$dir"/outs ]
        then
            echo "Cellranger for $dir exists"
            mkdir $dir
            cd $dir
            python3 /diskmnt/Projects/Users/austins2/tools/scrublet-comboRNA-auto.py -s $dir -c $cellrangerpath/$dir -u 0.2
            cd ..
        else
            echo "Cellranger outs for $dir does not exist at the provided location"
        fi
    fi
done

