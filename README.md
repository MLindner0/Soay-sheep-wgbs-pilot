# Soay sheep wgbs pilot data

This is a repository for the wgbs pilot data from the St. Kilda Soay Sheep Project. inlcuding the bioinformatic processing of the wgbs data, the preprocessing of the methylation count data, and the preliminary analysis of the methylation count data.

## Brief description of the wgbs dataset

The data set inlcudes sequencing data for **224 samples** (mostly from lambs with some repeat samples from sheep of older age). Wgbs sequencing was performed at the NEOF facility in Liverpool. The wgbs raw data (already adapter trimmed) is stored at `/shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed`. There are nine files per sample; samples were run on 3 lanes and data is split into R1, R2, and R0. (R0 files, reads without a mate, are not considered for analysis.)

## Bioinformatic analysis

### Preparations

To prepare the data for bioinformatic processing, I created symbolic links of the raw data files and renamed the links (see `pipeline-prepare_data.sh`) to have consistent file names between samples - **note**: never rename the raw data files, instead create symbolic links and rename the links. The renamed symbolic links are stored at `/shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_renamed_2`.

