"""
Pipeline for Nanopore.
"""
import os

# a wrapper rule
rule all:
    input:
        nanoplot_raw_output_dir = directory(os.path.join(os.getcwd(), 'nanoplot_raw')), # output of rule 'quality_control_raw_reads'
        clean_reads = 'clean_reads_nanopore.fastq.gz',  # output of rule 'filter_raw_reads'
        nanoplot_clean_output_dir = directory(os.path.join(os.getcwd(), 'nanoplot_clean')),  # output from rule 'quality_control_filtered_reads'
        sam_file = 'minimap2-nanopore.sam',     # output from rule 'mapping'
        bam_file = 'minimap2-nanopore.bam',    # output from the rule 'sam_to_bam'
        sorted_bam_file = 'minimap2-nanopore.sorted.bam',   # output from rule 'process_bam_file'
        indexed_bam_file = 'minimap2-nanopore.sorted.bam.bai',   # output from rule 'process_bam_file'
        primerclipped_sorted_bam_file = 'minimap2-nanopore.sorted.primerclipped.bam', # output from rule 'primer_clipping'
        vcf = 'please.vcf'
        #hdf = 'medaka-nanopore.consensus.hdf',    # output from rule 'variant_calling_preparation'
        #vcf = 'medaka-nanopore.vcf'    # output from rule 'variant_calling'
        #vcf = 'medaka-nanopore.vcf',
        #vcf_annotated = 'medaka-nanopore.annotate.vcf'

# quality control of raw reads
rule quality_control_raw_reads:
    input:
        nanopore_fastq_raw = os.path.join(os.getcwd(),'nanopore.fastq.gz')
    output:
        nanoplot_raw_output_dir = directory(os.path.join(os.getcwd(), 'nanoplot_raw'))
    shell:
        """
        mkdir -p {output.nanoplot_raw_output_dir}
        NanoPlot -t 4 --fastq {input.nanopore_fastq_raw} --title "Raw reads" --color darkslategrey --N50 --loglength -o {output.nanoplot_raw_output_dir}
        """

# filtering of raw reads
# THINK ABOUT PASSING THE filtlong PARAMETERS AS PARAMS INSIDE THE SNAKEMAKE RULE
rule filter_raw_reads:
    input:
        nanopore_fastq_raw = os.path.join(os.getcwd(),'nanopore.fastq.gz')
    output:
        #filtered_reads = directory("/home/qbaliu/envs/DAY1/PROJECT_2/SNAKEMAKE_POLIGON/FILTERING")
        clean_reads = 'clean_reads_nanopore.fastq.gz'
    shell:
        """
        filtlong --min_length 800 --max_length 1400 {input.nanopore_fastq_raw} | gzip - > {output.clean_reads}
        """ 


# quality control of filtered reads
rule quality_control_filtered_reads:
    input:
        clean_reads = "clean_reads_nanopore.fastq.gz"
    output:
        nanoplot_clean_output_dir = directory(os.path.join(os.getcwd(), 'nanoplot_clean'))
    shell:
        """
        mkdir -p {output.nanoplot_clean_output_dir}
        NanoPlot -t 4 --fastq {input.clean_reads} --title "Clean reads" --color darkslategrey --N50 --loglength -o {output.nanoplot_clean_output_dir}
        """


# mapping the filtered reads to the reference genome
rule mapping:
    input:
        clean_reads = 'clean_reads_nanopore.fastq.gz',
        reference = 'SARS_Wuhan_reference_genome.fasta'
    output:
        sam_file = 'minimap2-nanopore.sam'
    shell:
        """
        minimap2 -x map-ont -t 4 -a -o {output.sam_file} {input.reference} {input.clean_reads}
        """


# convert bam file to sam file
rule sam_to_bam:
    input:
        sam_file = 'minimap2-nanopore.sam'
    output:
        bam_file = 'minimap2-nanopore.bam'
    shell:
        """
        samtools view -bS {input.sam_file} > {output.bam_file}
        """

# sort and index the bam file
rule process_bam_file:
    input:
        bam_file = 'minimap2-nanopore.bam' 
    output:
        sorted_bam_file = 'minimap2-nanopore.sorted.bam',
        indexed_bam_file = 'minimap2-nanopore.sorted.bam.bai'
    shell:
        """
        samtools sort {input.bam_file} > {output.sorted_bam_file}
        samtools index {output.sorted_bam_file} > {output.indexed_bam_file}
        """


# primer clipping
rule primer_clipping:
    input:
        bed_file = 'nCoV-2019.bed',
        sorted_bam_file = 'minimap2-nanopore.sorted.bam'
    output:
        primerclipped_sorted_bam_file = 'minimap2-nanopore.sorted.primerclipped.bam'
    shell:
        """
        old_name=$(head -n 4 {input.bed_file}| tail -n 1 | awk '{{print $1}}')
        sed "s/$old_name/NC_045512.2/g" {input.bed_file} > nCov-2019.corrected.bed
        python primerbed2bedpe.py nCov-2019.corrected.bed --forward_identifier _LEFT --reverse_identifier _RIGHT -o nCoV-2019.bedpe
        bamclipper.sh -b {input.sorted_bam_file} -p nCoV-2019.bedpe -n 4
        """


rule variant_calling_please:
    input:
        sorted_primerclipped_bam_file = 'minimap2-nanopore.sorted.primerclipped.bam',
        #sorted_primerclipped_bam_file_index = 'minimap2-nanopore.sorted.primerclipped.bam.bai',
        reference = 'SARS_Wuhan_reference_genome.fasta'
    output:
        vcf = 'please.vcf'
    shell:
        """
        java -jar pilon-1.24.jar --genome {input.reference} --frags {input.sorted_primerclipped_bam_file} --output {output.vcf} --threads 4 --vcf
        """

# variant calling
'''
# this works well
rule variant_calling_preparation:
    input:
        primerclipped_sorted_bam_file = 'minimap2-nanopore.sorted.primerclipped.bam',
        reference = 'SARS_Wuhan_reference_genome.fasta'
    output:
        #vcf_annotated = 'medaka-nanopore.annotate.vcf',
        hdf = 'medaka-nanopore.consensus.hdf'
        #vcf = 'medaka-nanopore.vcf'
    shell:
        
        medaka consensus --model r103_fast_snp_g507 --threads 4 --chunk_len 800 --chunk_ovlp 400 {input.primerclipped_sorted_bam_file} {output.hdf}
        

# this doesnt work well
rule variant_calling:
    input:
        reference = 'SARS_Wuhan_reference_genome.fasta',
        hdf = 'medaka-nanopore.consensus.hdf'
    output:
        vcf = 'medaka-nanopore.vcf'
    shell:
        
        medaka variant -i {input.reference} -d {input.hdf} -o {output.vcf}
        
'''