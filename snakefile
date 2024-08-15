## snakemake pipeline ##

# bioinformatic processing of wgbs pilot data from the St. Kilda Soay Sheep Project
# The dataset inlcudes 224 samples sequenced in 3 aliquots each
# QC is performed outside the pipeline
# In silico conversion of reference genome (Rambouillet v2.0) is performed outside the pipeline using Bismark
# Trimming, alignment, and deduplication are performed for each aliquot (after which read groups are added). Then, alignments are merged for methylation calling.
# Alignment stats, genome coverage at multiple base coverage thresholds, and estimates of bisulfite conversion efficiency are extracted from the alignment files

# EXEPTION: the pipeline can be used to call methylationfrom all aliquots (i.e. 3x per sample) for 10 samples (i.e. 30 aliquots in total). This cam be useful to assess the repeatability of aliquots. To run the pipeline this way, the wildcards and the rule all must be adjusted. (Alternative definition of wildcards and the rule all are included in the script but "out-commented".)

## TOOLS ##
# Trimming: trim_galore v0.6.10
# Alignment: Bismark v0.24.0
# Deduplication: Bismark v0.24.0
# Addition of read groups: Picard v2.20.4
# Merging of alignments: Picard v2.20.4
# Methylation calling: Bismark v0.24.0
# Formatting alignments and extracting information from alignments: Samtools v1.9 & custom code

## prepare pipeline
# load python library
import pickle

# 1. get read group information

# load .pkl file with read group information
# .pkl file was created using make_.pkl.py
with open('/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/src/dict_ReadGroups.pkl', 'rb') as fp:
    dict_ReadGroups = pickle.load(fp)

print('dict_ReadGroups length: ' + str(len(dict_ReadGroups)))

# 2. define wildcards:
# Wild cards "run" (aliquot) and "name" (sample)
RUN, = glob_wildcards("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_2/{run}_R1_val_1.fq.gz")
NAME, = glob_wildcards("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_2/{name}_1_R1_val_1.fq.gz")

# Exeption: here, you find the wildcard definition for calling methylation from all aliquots (i.e. 3x per sample) for 10 samples (i.e. 30x aliquots in total).
# RUN = ["1-1_1","1-1_2","1-1_3","10-11_1","10-11_2","10-11_3","100-87_1","100-87_2","100-87_3","101-88_1","101-88_2","101-88_3","102-89_1","102-89_2","102-89_3","3-3_1","3-3_2","3-3_3","30-16_1","30-16_2","30-16_3","31-17_1","31-17_2","31-17_3","9-10_1","9-10_2","9-10_3","90-77_1","90-77_2","90-77_3"]
# NAME = ["1-1","10-11","100-87","101-88","102-89","3-3","30-16","31-17","9-10","90-77"]

# 3. constrain wildcards
wildcard_constraints:
    run = '|'.join([re.escape(x) for x in RUN]),
    name = '|'.join([re.escape(x) for x in NAME])

# define rule all
rule all:
    input:
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/meth/{name}.deduplicated_merged.bismark.cov.gz", name=NAME),
        # expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/meth/{run}.deduplicated.bismark.cov.gz", run=RUN),
        # expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}_samtools_stats_out.txt", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.AverageCoverage", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_10x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_05x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_15x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_20x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_25x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_30x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_02x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_04x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_06x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_08x", name=NAME),
        expand("/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.BS-conversion", run=RUN)

# the following rules are performed on each aliquot per sample (i.e. 3x per sample)

# run trimming of raw sequencing reads. Note that adapter trimming was already performed at Liverpool.
rule trimming:
    input:
        R1="/shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_renamed_2/{run}_R1.fastq.gz",
        R2="/shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_renamed_2/{run}_R2.fastq.gz"
    output:
        R1Out="/fastdata/bi1ml/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_2/{run}_R1_val_1.fq.gz",
        R2Out="/fastdata/bi1ml/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_2/{run}_R2_val_2.fq.gz"
    log:
        "/shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/logs/{run}_trimming.log"
    threads: 4
    shell:
        "(trim_galore --2colour 20 -j 4 --paired --gzip -o /fastdata/bi1ml/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_2/ {input.R1} {input.R2}) 2> {log}"

# run alignments
rule alignment:
    input:
        R1Tr="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_2/{run}_R1_val_1.fq.gz",
        R2Tr="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/Trimmed_2/{run}_R2_val_2.fq.gz"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}_R1_val_1_bismark_bt2_pe.bam",
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}_R1_val_1_bismark_bt2_PE_report.txt"
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{run}_alignment.log"
    threads: 4
    shell:
        "(bismark -X 1000 --parallel 4 --genome /mnt/parscratch/users/bi1ml/public/genomes/ARS-UI_Ramb_v2.0/GCF_016772045.1/ -1 {input.R1Tr} -2 {input.R2Tr} --temp_dir /mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/temp/ -o /mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/) 2> {log}"

# prior to deduplication, alignments must be converted to the .sam format and sorted to endure optimal performance of bismark during deduplication.

# convert to sam and sort
rule samsort:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}_R1_val_1_bismark_bt2_pe.bam"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.sam"
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{run}_sam.log"
    threads: 1
    shell:
        "(samtools view -h {input} | samtools sort -n -O sam -T /mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/temp/ > {output}) 2> {log}"

# run deduplication
rule dedup:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.sam"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.deduplicated.bam"
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{run}_deduplication.log"
    threads: 1
    shell:
        "(deduplicate_bismark -p --output_dir /mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/ -o {wildcards.run} {input}) 2> {log}"

# convert alignments back to .bam format
rule samsort2:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.deduplicated.bam"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.deduplicated.coordinate.bam"
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{run}_samsort2.log"
    threads: 1
    shell:
        "(samtools view -h {input} | samtools sort -O bam -T /mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/temp/ > {output}) 2> {log}"

# add read groups to alignments using the dictionary created in lines 28-31
rule readgr:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.deduplicated.coordinate.bam"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.deduplicated.withRG.bam"
    params:
        Lane=lambda wildcards: dict_ReadGroups[wildcards.run]['Lane'],
        Flowcell=lambda wildcards: dict_ReadGroups[wildcards.run]['Flowcell']
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{run}_RG.log"
    threads: 1
    shell:
        "(picard -Xmx4096m AddOrReplaceReadGroups I={input} O={output} TMP_DIR=/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/temp/ ID={params.Flowcell}.{params.Lane} LB=NEB_EM-Seq PL=illumina PU={params.Flowcell}.{params.Lane} SM={wildcards.run} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=queryname) 2> {log}"

# the alignments of all aliquots per sample (i.e. 3x per sample) are now merged, which means that the following rules are performed on each sample (i.e. 1x per sample)

## EXEPTION: if you call methylation from all aliquots (i.e. 3x per sample) for 10 samples (i.e. 30x aliquots in total), rule merge is skipped and methylation is called for each of the aliquots separately (i.e. 30x in total). 

# merge alignments of all aliquots per sample
rule merge:
    input:
       R1="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}_1.deduplicated.withRG.bam",
       R2="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}_2.deduplicated.withRG.bam",
       R3="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}_3.deduplicated.withRG.bam"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_merged.bam"
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_merge.log"
    threads: 1
    shell:
        "(picard -Xmx4096m MergeSamFiles I={input.R1} I={input.R2} I={input.R3} O={output} TMP_DIR=/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/temp/ SORT_ORDER=queryname) 2> {log}"

# call methylation from merged alignments
rule meth:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_merged.bam"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/meth/{name}.deduplicated_merged.bismark.cov.gz"
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_methcalling_test.log"
    threads: 4
    shell:
        "(bismark_methylation_extractor -p --parallel 4 --no_overlap --report --bedGraph --scaffolds --cytosine_report --ignore 3 --ignore_r2 4 --ignore_3prime 3 --ignore_3prime_r2 2 -o /mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/meth/ --genome_folder /mnt/parscratch/users/bi1ml/public/genomes/ARS-UI_Ramb_v2.0/GCF_016772045.1/ {input}) 2> {log}"

# call methylation from alignments of all aliquots (i.e. 3x per sample) for 10 samples
# note, while this is the "EXEPTION" rule, there is no need to "out-comment" the rule as it won't be executed unless the output files are defined in the rule all
rule meth_repeatability:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.deduplicated.bam"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/meth/{run}.deduplicated.bismark.cov.gz"
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{run}_methcalling_test.log"
    threads: 4
    shell:
        "(bismark_methylation_extractor -p --parallel 4 --no_overlap --report --bedGraph --scaffolds --cytosine_report --ignore 3 --ignore_r2 4 --ignore_3prime 3 --ignore_3prime_r2 2 -o /mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/meth/ --genome_folder /mnt/parscratch/users/bi1ml/public/genomes/ARS-UI_Ramb_v2.0/GCF_016772045.1/ {input}) 2> {log}"

# Extract (1) alignment stats, (1) genome coverage at multiple base coverage thresholds, and (3) estimates of bisulfite conversion efficiency from the alignments

# (1) alignment stats

# run samtools (rule sam_stats) and format the output (rule get_sam_stats)
rule sam_stats:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_merged.bam"
    output:
         "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.samtools_stats.txt"
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_samstats1.log"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"

rule get_sam_stats:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.samtools_stats.txt"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}_samtools_stats_out.txt"
    shell:
        """
        cat {input} | grep ^SN | awk -F '\t' -v OFS='\t' '{{print $2,$3}}' 1> {output}
        """
# extract average coverage and (2) genome coverage at multiple base coverage thresholds (i.e. the number of bases covered at a threshold divide by genome length)
# sort alignments (rule cov_sort), extract average coverage (samtools depth) and genome coverage (samtools mpileup) at coverage threshold of 10x (rule coverage10)
# repeat extraction of genome coverage at various coverage thresholds (5x, 15x, 20x, 25x, 30x & 2x, 4x, 6x, 8x)

rule cov_sort:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_merged.bam"
    output:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    log:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_picard_sort.log"
    shell:
        "(picard -Xmx4096m SortSam I={input} O={output} SORT_ORDER=coordinate TMP_DIR=/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/temp/) 2> {log}"

rule coverage10:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        depth="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.AverageCoverage",
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_10x"
    log:
        depth="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_depth.log",
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_10x.log"
    shell:
        """
        samtools depth {input} | awk '{{sum+=$3}} END {{ print "Average = ",sum/NR}}' 1> {output.depth} 2> {log.depth}
        MIN_COVERAGE_DEPTH=10
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

rule coverage05:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_05x"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_05x.log"
    shell:
        """
        MIN_COVERAGE_DEPTH=5
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

rule coverage15:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_15x"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_15x.log"
    shell:
        """
        MIN_COVERAGE_DEPTH=15
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

rule coverage20:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_20x"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_20x.log"
    shell:
        """
        MIN_COVERAGE_DEPTH=20
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

rule coverage25:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_25x"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_25x.log"
    shell:
        """
        MIN_COVERAGE_DEPTH=25
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

rule coverage30:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_30x"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_30x.log"
    shell:
        """
        MIN_COVERAGE_DEPTH=30
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

rule coverage02:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_02x"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_02x.log"
    shell:
        """
        MIN_COVERAGE_DEPTH=2
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

rule coverage04:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_04x"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_04x.log"
    shell:
        """
        MIN_COVERAGE_DEPTH=4
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

rule coverage06:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_06x"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_06x.log"
    shell:
        """
        MIN_COVERAGE_DEPTH=6
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

rule coverage08:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.deduplicated_sorted_coordinate.bam"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{name}.BasesCovered_08x"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{name}_sam_mpileup_08x.log"
    shell:
        """
        MIN_COVERAGE_DEPTH=8
        samtools mpileup {input} | awk -v X="${{MIN_COVERAGE_DEPTH}}" '$4>=X' | wc -l 1> {output.mpileup} 2> {log.mpileup}
        """

# Extract conservative estimates of bisulfite conversion efficiency from the alignments based on non-CpG methylation in alignment report

rule BS_Conversion:
    input:
        "/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}_R1_val_1_bismark_bt2_PE_report.txt"
    output:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/alignment/{run}.BS-conversion"
    log:
        mpileup="/mnt/parscratch/users/bi1ml/public/methylated_soay/soay_wgbs_pilot_mar2023/logs/{run}_BS.conv.log"
    shell:
        """
        grep 'C methylated in CHG context:' {input} | sed -e 's/%//g' | awk -v OFS='\t' '{{print $6,(100-$6),"{wildcards.run}"}}' > {output}
        """
