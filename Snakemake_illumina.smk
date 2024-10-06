"""
Illumina pipeline
"""
import os

# wrapper rule
rule all:
    input:
        report_R1 = 'illumina.R1_fastqc.html',  # output from rule 'quality_control_raw_reads'
        report_R2 = 'illumina.R2_fastqc.html',   # output from rule 'quality_control_raw_reads'
        clean_reads_R1 = 'clean_reads.R1.fastq.gz',  # output from the rule 'trimming_raw_reads'
        clean_reads_R2 = 'clean_reads.R1.fastq.gz',   # output from the rule 'trimming_raw_reads'
        report_clean_R1 = 'clean_reads.R1_fastqc.html',  # output from rule 'quality_control_clean_reads'
        report_clean_R2 = 'clean_reads.R2_fastqc.html',  # output from rule 'quality_control_clean_reads'
        zip_folder_clean_R1 = 'clean_reads.R1_fastqc.zip',  # output from rule 'quality_control_clean_reads'
        zip_folder_clean_R2 = 'clean_reads.R2_fastqc.zip',  # output from rule 'quality_control_clean_reads'
        sam_file = 'minimap2-illumina.sam',   # output from rule 'mapping'
        bam_file = 'illumina.bam',  # output from the rule 'convert_sam_to_bam'
        sorted_bam_file = 'illumina.sorted.bam',  # output from the rule 'convert_sam_to_bam'
        indexed_sorted_bam_file = 'illumina.sorted.bam.bai',  # output from the rule 'convert_sam_to_bam'
        amplicons_corrected = 'cleanplex-corrected.amplicons.bedpe',   # output from rule 'primer_clipping'
        vcf = 'freebayes-illumina.vcf',   # output from rule 'variant_calling'
        vcf_gzipped = 'freebayes-illumina.vcf.gz',  # output from rule 'consensus_generation'
        consensus = 'consensus-illumina.fasta'    #  # output from rule 'consensus_generation'


# quality control of raw reads
rule quality_control_raw_reads:
    input:
        R1 = 'illumina.R1.fastq.gz',
        R2 = 'illumina.R2.fastq.gz'
    output:
        report_R1 = 'illumina.R1_fastqc.html',
        report_R2 = 'illumina.R2_fastqc.html'
    shell:
        """
        fastqc -t 4 {input.R1} {input.R2}
        """

# trimming of raw reads
rule trimming_raw_reads:
    input:
        R1 = 'illumina.R1.fastq.gz',
        R2 = 'illumina.R2.fastq.gz'
    output:
        clean_reads_R1 = 'clean_reads.R1.fastq.gz',
        clean_reads_R2 = 'clean_reads.R2.fastq.gz'

    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o clean_reads.R1.fastq.gz -O clean_reads.R2.fastq.gz --thread 4 --qualified_quality_phred 20 --length_required 50
        """


# quality control of trimmed reads
rule quality_control_clean_reads:
    input:
        clean_reads_R1 = 'clean_reads.R1.fastq.gz',
        clean_reads_R2 = 'clean_reads.R2.fastq.gz'
    output:
        report_clean_R1 = 'clean_reads.R1_fastqc.html',
        report_clean_R2 = 'clean_reads.R2_fastqc.html',
        zip_folder_clean_R1 = 'clean_reads.R1_fastqc.zip',
        zip_folder_clean_R2 = 'clean_reads.R2_fastqc.zip'
    shell:
        """
        fastqc -t 2 {input.clean_reads_R1}
        fastqc -t 2 {input.clean_reads_R2}
        """


# mapping
rule mapping:
    input:
        reference = 'SARS_Wuhan_reference_genome.fasta',
        clean_reads_R1 = 'clean_reads.R1.fastq.gz',
        clean_reads_R2 = 'clean_reads.R2.fastq.gz',
    output:
        sam_file = 'minimap2-illumina.sam'
    shell:
        """
        minimap2 -x sr -t 4 -a -o {output.sam_file} {input.reference} {input.clean_reads_R1} {input.clean_reads_R2}
        """


# convert sam file to bam file and process the bam file
rule convert_sam_to_bam:
    input:
        sam_file = 'minimap2-illumina.sam'
    output:
        bam_file = 'illumina.bam',
        sorted_bam_file = 'illumina.sorted.bam',
        indexed_sorted_bam_file = 'illumina.sorted.bam.bai'
    shell:
        """
        samtools view -bS {input.sam_file} > {output.bam_file}
        samtools sort {output.bam_file} > {output.sorted_bam_file}
        samtools index {output.sorted_bam_file} > {output.indexed_sorted_bam_file}
        """


# primer clipping
rule primer_clipping:
    input:
        amplicons = 'cleanplex.amplicons.bedpe',
        sorted_bam_file = 'illumina.sorted.bam'
    output:
        amplicons_corrected = 'cleanplex-corrected.amplicons.bedpe'
    shell:
        """
        bedpe_name=$(head -n 4 nCoV-2019.bedpe| tail -n 1 | awk '{{print $1}}')
        sed "s/$bedpe_name/NC_045512.2/g" {input.amplicons} > {output.amplicons_corrected}
        bamclipper.sh -b {input.sorted_bam_file} -p {output.amplicons_corrected} -n 4
        """


# variant calling
rule variant_calling:
    input:
        reference = 'SARS_Wuhan_reference_genome.fasta',
        primerclipped_bam_file = 'illumina.sorted.primerclipped.bam'
    output:
        vcf = 'freebayes-illumina.vcf'
    shell:
        """
        samtools faidx {input.reference}

        freebayes -f {input.reference} --min-alternate-count 10 \
        --min-alternate-fraction 0.1 --min-coverage 20 --pooled-continuous \
        --haplotype-length -1 {input.primerclipped_bam_file} > {output.vcf}
        """


rule consensus_generation:
    input:
        vcf = 'freebayes-illumina.vcf',
        reference = 'SARS_Wuhan_reference_genome.fasta'
    output:
        vcf_gzipped = 'freebayes-illumina.vcf.gz',
        consensus = 'consensus-illumina.fasta'
    shell:
        """
        bgzip -f {input.vcf} -o {output.vcf_gzipped}
        tabix -f -p vcf {output.vcf_gzipped}
        bcftools consensus -f {input.reference} {output.vcf_gzipped} -o {output.consensus}
        old_header=$(grep '>' {output.consensus} | awk '{{print $1}}')
        sed -i "s/$old_header/>Consensus-Nanopore/g" {output.consensus}
        """