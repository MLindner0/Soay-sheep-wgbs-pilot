#!/bin/sh

## prepare sequencing data files for the bioinformatics pipelines

# activate conda env
module load Anaconda3
source activate rename

# first, create symbolic links (THE RAW DATA IS NEVER CHNAGED)
for f in *.gz; do ln -s /shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed/$f /shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_renamed/; done
# Second, remove the tissue abbreviation and extra numbers (for repeat extractions that are not consistent across samples)
cd /shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_renamed/

rename 's/_1_BC//' *.gz #remove '_1_BC'
rename 's/_1_NB//' *.gz #remove '_1_NB'
rename 's/_BC//' *.gz #remove '_BC'
rename 's/_ET//' *.gz #remove '_ET'
rename 's/_NB//' *.gz #remove '_NB'

# deactivate conda env
conda deactivate

# third, exchange the lane identifies (3 per sample) to 1, 2 and 3 to be consistent across samples
# this time, rename symbolic links based on .txt file
while read line ; do 
    cp -R $line 
done < src/FileRename.txt # see pipeline-data_from_filenames.R for how file was created
