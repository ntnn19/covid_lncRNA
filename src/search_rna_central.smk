from Bio import SeqIO
import os
output_dir = "output/search_rna_central"
fasta_files = [os.path.join("output", "annotate_probe_seqs/vds_orphan_lncrna.fasta"),  os.path.join("output","annotate_probe_seqs/vds_orphan_mrna.fasta")]

for f, biotype in zip(fasta_files,["lncrna","mrna"]):
    os.makedirs(output_dir,exist_ok=True)
    fasta = SeqIO.parse(f,"fasta")
    for rec in fasta:
        with open(os.path.join(output_dir,f"{rec.id}_biotype-{biotype}.fasta"),"w") as record_file:
            rec.id = f"{rec.id}_biotype-{biotype}"
            rec.description = f""
            SeqIO.write(rec,record_file,"fasta")

FASTA_NAMES, = glob_wildcards(os.path.join(output_dir,"{fasta_name_}.fasta"))
print(FASTA_NAMES[:10])

rule all:
    input:
        expand(os.path.join(output_dir,"rule_SEARCH_RNA_CENTRAL","{fasta_name_}.json"), fasta_name_=FASTA_NAMES),
        expand(os.path.join(output_dir,"rule_EXTRACT_TOP_HIT","{fasta_name_}.tsv"), fasta_name_=FASTA_NAMES)

rule SEARCH_RNA_CENTRAL:
    input:
        os.path.join(output_dir,"{fasta_name_}.fasta")
    output:
        os.path.join(output_dir,"rule_SEARCH_RNA_CENTRAL","{fasta_name_}.json")
    shell:
        """
        python src/search_rna_central.py {input}
        mv results/{wildcards.fasta_name_}.json {output}
        """

rule EXTRACT_TOP_HIT:
    input:
        rules.SEARCH_RNA_CENTRAL.output[0]
    output:
        os.path.join(output_dir,"rule_EXTRACT_TOP_HIT","{fasta_name_}.tsv")
    shell:
        """
        python src/extract_top_hit.py extract-top-hit {input} --format tsv -f rnacentral_id --fields rna_type  --fields identity --fields query_coverage --species human > {output}
        """
