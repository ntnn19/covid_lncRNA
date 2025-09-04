# Snakefile: BWA exact matching workflow with uppercase rules and rule-specific output folders

import os

# Config
OUTPUT_DIR   = "output/align_to_hg38"
GENOME_FASTA = "input/fasta/GRCh38.p14.genome.fa"
SOURCES,  = glob_wildcards("output/annotate_probe_seqs/{source}.fasta")

rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "rule_SAM_TO_TSV", "{source}_alignments.tsv"),source=SOURCES)

# -----------------------------------
# Build BWA index
# -----------------------------------
rule BWA_INDEX:
    input:
        fasta = GENOME_FASTA
    container:
        "~/singularity_containers/bwa/bwa.sif"
    output:
        expand(os.path.join(OUTPUT_DIR, "rule_BWA_INDEX", "GRCh38.p14.genome.fa.{ext}"),
               ext=["amb","ann","bwt","pac","sa"])
    shell:
        """
        bwa index {input.fasta}
        mv input/fasta/GRCh38.p14.genome.fa.* {OUTPUT_DIR}/rule_BWA_INDEX
        mv input/fasta/GRCh38.p14.genome.fa {OUTPUT_DIR}/rule_BWA_INDEX
        """

# -----------------------------------
# Align queries with exact matches
# -----------------------------------
rule BWA_ALIGN:
    input:
        fasta = "output/annotate_probe_seqs/{source}.fasta",
        index = rules.BWA_INDEX.output
    output:
        sam = os.path.join(OUTPUT_DIR, "rule_BWA_ALIGN", "{source}_alignments.sam")
    container:
        "~/singularity_containers/bwa/bwa.sif"
    params:
        sai    = os.path.join(OUTPUT_DIR, "rule_BWA_ALIGN", "{source}_temp.sai"),
        genome    = os.path.join(OUTPUT_DIR, "rule_BWA_INDEX", "GRCh38.p14.genome.fa")
    shell:
        """
        bwa aln -n 0 {params.genome} {input.fasta} > {params.sai}
        bwa samse {params.genome} {params.sai} {input.fasta} > {output.sam}
        rm {params.sai}
        """

# -----------------------------------
# Convert SAM to TSV
# -----------------------------------
rule SAM_TO_TSV:
    input:
        sam = rules.BWA_ALIGN.output.sam
    container:
        "~/singularity_containers/bioawk/bioawk.sif"
    output:
        tsv = os.path.join(OUTPUT_DIR, "rule_SAM_TO_TSV", "{source}_alignments.tsv")
    shell:
        """
        bioawk -c sam '{{print $qname, $rname, $pos, $seq}}' {input.sam} > {output.tsv}
        """
