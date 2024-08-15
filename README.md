# Soay sheep wgbs pilot data

This is a repository for the wgbs pilot data from the St. Kilda Soay Sheep Project. inlcuding the bioinformatic processing of the wgbs data, the preprocessing of the methylation count data, and the preliminary analysis of the methylation count data.

## Brief description of the wgbs dataset

The data set inlcudes sequencing data for **224 samples** (mostly from lambs with some repeat samples from sheep of older age). Wgbs sequencing was performed at the NEOF facility in Liverpool. The wgbs raw data (already adapter trimmed) is stored at `/shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed`. There are nine files per sample; samples were run on 3 lanes and data is split into R1, R2, and R0. (R0 files, reads without a mate, are not considered for analysis.)

## Bioinformatic analysis

### Preparations

#### raw data
To prepare the data for bioinformatic processing, I created symbolic links of the raw data files and renamed the links (see `src/pipeline-prepare_data.sh` & `pipeline-data_from_filenames.R`) to have consistent file names between samples - **note**: never rename the raw data files, instead create symbolic links and rename the links. The renamed symbolic links are stored at `/shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_renamed_2`.

Additionally, I extracted information from the file names regarding read groups (lane, library, etc., see `pipeline-data_from_filenames.R` & `make_.pkl.py`) to add read group information to the alignments of all aliquots (3x per sample) before all aliquots per sample are merged into one alignment file. 

#### reference genome
To perform the alignment to wgbs data, a in silico bisulfite conversion of the reference genome must be perfromed. First, a copy of the reference genome with the .fa extension is created:
```
cp /shared/slate_group1/Shared/genomes/ARS-UI_Ramb_v2.0/GCF_016772045.1/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna /shared/slate_group1/Shared/genomes/ARS-UI_Ramb_v2.0/GCF_016772045.1/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fa
```
Second, the in silico bisulfite conversion of the reference genome is performed with Bismark (v0.24.0):
```
bismark_genome_preparation /shared/slate_group1/Shared/genomes/ARS-UI_Ramb_v2.0/GCF_016772045.1/
```

### Run snakemake pipeline

The snakemake pipeline (`snakefile`) includes the trimming, alignment, deduplication, and addition of read group information for all aliquots (i.e. 3x per sample). Then, alignments are merged for methylation calling. Additionally, the pipeline extracts alignment stats, average base coverage (depth), genome coverage (breadth) at multiple base coverage thresholds, and conservative estimates of the bisulfite conversion efficiency from the alignment files.

To run the pipeline, I executed a batch job on Stanage's HPC using the follwoing submission script:

```
#!/bin/bash

# load conda and activate environment
module load Anaconda3
source activate methylation_calling

# load snakemake
module use /usr/local/modulefiles/staging/eb/all
module load snakemake/6.1.0-foss-2020b

# run the pipeline
snakemake -j 20
```

**Note**: Before running the pipeline perform a dry-run to make sure the pipeline is executable and executes the correct rules: `snakemake -n`

To submit the script to the HPC I used the following specifications in the submission file:
```
#!/bin/bash

#SBATCH --comment=run_pipeline
#SBATCH --output=logs/run_pipeline.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100GB
#SBATCH --time=96:00:00
#SBATCH --mail-user=my.name@sheffield.ac.uk
#SBATCH --mail-type=all

submit/script.sh
```

**Note**: If the batch job is terminated before the pipeline is finished (e.g. out of memory, time out etc.), make sure that there are no unfinished output files and unlock the pipeline before you submit the batch job (with adjusted specifications) again: `snakemake --unlock`


