SAMPLES, = glob_wildcards("../output/convert_fasta_to_fastq/{sample}.fastq.gz")
def get_fq(wildcards):
    return {
        "fq1": "../output/convert_fasta_to_fastq/{sample}.fastq.gz".format(**wildcards)
    }

rule all:
    input:
        expand("results/star/{sample}/Aligned.sortedByCoord.out.bam",sample=SAMPLES),
        expand("results/star/{sample}/Aligned.sortedByCoord.out.sam",sample=SAMPLES),
        #expand("results/star/{sample}/ReadsPerGene.out.tab",sample=SAMPLES)

rule align:
    input:
        unpack(get_fq),
        index="../output/align_to_hg38_star/index",
        gtf="../input/gtf/GCF_000001405.40_GRCh38.p14_genomic.gtf",
    output:
        aln="results/star/{sample}/Aligned.sortedByCoord.out.bam",
    log:
        "logs/star/{sample}.log",
    params:
        idx=lambda wc, input: input.index,
        extra=lambda wc, input: f'--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {input.gtf} {config["params"]["star"]}',
    threads: 24
    container:
        "/home/nadmin/singularity_containers/star/star.sif"
    shell:
        """
        STAR  --runThreadN {threads} --genomeDir /index --readFilesIn  /annotate_probe_seq/{wildcards.sample}.fastq.gz  --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile /gtf/Homo_sapiens.GRCh38.115.gtf  --outTmpDir /tmp/tmp_cvhehzg/STARtmp --outFileNamePrefix /tmp/tmp_cvhehzg/ --outStd BAM_SortedByCoordinate > results/star/{wildcards.sample}/Aligned.sortedByCoord.out.bam  2> logs/star/vds_orphan_lncrna.log
        """

rule sam:
    input:
        rules.align.output[0]
    output:
        sam="results/star/{sample}/Aligned.sortedByCoord.out.sam",
    log:
        "logs/star/{sample}_sam.log",
    threads: 24
    shell:
        """
        samtools view -h {input} > {output}
        """
