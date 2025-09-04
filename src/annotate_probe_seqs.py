

import click
# Save original echo
_original_echo = click.echo

def extract_gene_id(t):
    if not t:  # empty tuple
        return None
    results = []
    for item in t:
        if isinstance(item, str):
            parts = item.split("|")
            if len(parts) > 1:
                results.append(parts[1])
    return tuple(results) if results else ()

def contains_any_queries(x, queries):
    if isinstance(x, str):
        return any(q in x for q in queries)
    elif isinstance(x, (tuple, list)):
        return any(any(q in str(item) for q in queries) for item in x)
    else:
        return False

def timestamped_echo(message=None, **kwargs):
    ts = datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    if message is None:
        _original_echo(ts, **kwargs)
    else:
        _original_echo(f"{ts} {message}", **kwargs)

def find_matches(sequence_to_find, seq_dict,source="GENCODE"):
    """
    Search sequences in FASTA_FILE for any of the query sequences provided with --seq.
    Prints the IDs of matching sequences.
    """
    # Load sequences into a dict (id -> SeqRecord)
    if source == "GENCODE":
        matches = [key for key in seq_dict if sequence_to_find in str(seq_dict[key].seq)]
    if source == "REFSEQ":
        matches = [key for key in seq_dict if sequence_to_find in str(seq_dict[key].seq)]
    #matches = list(map(lambda x: x if sequence_to_find in str(seq_dict[x].seq) else None, list(seq_dict.keys())))
    #matches = [m for m in matches if m]
    return matches

def majority_vote_with_tie(lst):
    if not lst:
        return None
    counts = Counter(lst)
    max_count = max(counts.values())
    winners = [item for item, count in counts.items() if count == max_count]
    if len(winners) == 1:
        return [winners[0]]
    else:
        return sorted(winners)


def find_annotations(list_of_hits,source="GENCODE",seq_dict=None):
    """
    Search sequences in FASTA_FILE for any of the query sequences provided with --seq.
    Prints the IDs of matching sequences.
    """
    # Load sequences into a dict (id -> SeqRecord)
    if source == "GENCODE":
        matches = list(map(lambda x: x.split("|")[-2] if x else [], list_of_hits))
    if source == "REFSEQ":
        matches = list(map(lambda x: seq_dict[x].description.split(",")[-1].strip() if x else [], list_of_hits))
    return matches


click.echo = timestamped_echo

@click.command()
@click.argument('source', type=str)  #GENECODE: wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz, wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.basic.annotation.gtf.gz
@click.argument('lncrna', type=click.Path(exists=True), required=True)
@click.argument('mrna', type=click.Path(exists=True), required=True)
def main(source,lncrna,mrna):
    pandarallel.initialize(progress_bar=True)
    click.echo("Loading lncRNA and mRNA tables")
    lncrna_df = pd.read_excel(
        lncrna,
        skiprows=66)
    mrna_df = pd.read_excel(
        mrna,
        skiprows=37)
    click.echo("lncRNA and mRNA tables loaded")
    vds_lncrna_annotation_df = pd.DataFrame(lncrna_df["Sequence"].drop_duplicates().copy())
    vds_mrnarna_annotation_df = pd.DataFrame(mrna_df["Sequence"].drop_duplicates().copy())
    gencode_fasta = "input/fasta/gencode.v48.transcripts.fa"
    refseq_fasta = "input/fasta/GCF_000001405.40_GRCh38.p14_rna.fna"
    gencode_seq_d = SeqIO.to_dict(SeqIO.parse(gencode_fasta, "fasta"))
    refseq_seq_d = SeqIO.to_dict(SeqIO.parse(refseq_fasta, "fasta"))
    gencode_df = pd.DataFrame.from_dict({key:str(gencode_seq_d[key].seq) for key in gencode_seq_d.keys()},orient="index").reset_index()
    refseq_df = pd.DataFrame.from_dict({key:str(refseq_seq_d[key].seq) for key in refseq_seq_d.keys()},orient="index").reset_index()
    gencode_df.columns =["id","Sequence"]
    refseq_df.columns =["id","Sequence"]
    if source == "GENCODE" or source == "gencode":
        if not os.path.exists("output/annotate_probe_seqs/vds_lncrna_annotation.csv") or not os.path.exists("output/annotate_probe_seqs/vds_mrna_annotation.csv"):
            for df,out_path in zip([vds_lncrna_annotation_df,vds_mrnarna_annotation_df],["output/annotate_probe_seqs/vds_lncrna_annotation.csv","output/annotate_probe_seqs/vds_mrna_annotation.csv"]):
                click.echo(f"Writing VDS annotations to output/annotate_probe_seqs")
                click.echo(f"Searching for matches")
                df["GENCODE_hit"] = df["Sequence"].parallel_apply(lambda x: find_matches(x, gencode_seq_d))
                click.echo(f"Search completed")
                df["GENCODE_hit"] = df["GENCODE_hit"].apply(lambda x: tuple(x) if x else ())
                df["n_GENCODE_hits"] = df["GENCODE_hit"].apply(len)
                df["GENCODE_annotations"] = df["GENCODE_hit"].apply(lambda y: find_annotations(y))
                df["GENCODE_majority_vote_annotation"] = df["GENCODE_annotations"].apply(majority_vote_with_tie)
                df["GENCODE_majority_vote_annotation"] = df["GENCODE_majority_vote_annotation"].apply(lambda x: tuple(x) if x else ())
                annotations_with_lncrna = [a for a in df["GENCODE_majority_vote_annotation"].value_counts().index if "lncRNA" in a and "protein_coding" not in a and "protein_coding_CDS_not_defined" not in a]

                click.echo(f"Number of sequences annotated by GENCODE as lncRNA and not as protein_coding annotations:\n{df.GENCODE_majority_vote_annotation.value_counts().loc[annotations_with_lncrna]}")
                df["GENCODE_hit"] = df["GENCODE_hit"].apply(lambda x: tuple(x) if x else ())
                df["GENCODE_annotations"] = df["GENCODE_annotations"].apply(
                    lambda x: tuple(x) if x else ())

                df['GENCODE_gene_id'] = df['GENCODE_hit'].apply(extract_gene_id)
                df.to_csv(out_path, index=False)
        else:
            lncrna_vds_df = pd.read_csv("output/annotate_probe_seqs/vds_lncrna_annotation.csv")
            mrnarna_vds_df = pd.read_csv("output/annotate_probe_seqs/vds_mrna_annotation.csv")
            lncrna_vds_not_in_gencode_df = lncrna_vds_df[lncrna_vds_df.n_GENCODE_hits==0].copy()
            mrnarna_vds_not_in_gencode_df = mrnarna_vds_df[mrnarna_vds_df.n_GENCODE_hits==0].copy()
            if not os.path.exists("output/annotate_probe_seqs/vds_lncrna_refseq_annotation.csv") or not os.path.exists(
                    "output/annotate_probe_seqs/vds_mrna_refseq_annotation.csv"):
                for df, out_path in zip([lncrna_vds_not_in_gencode_df, mrnarna_vds_not_in_gencode_df],
                                        ["output/annotate_probe_seqs/vds_lncrna_refseq_annotation.csv",
                                         "output/annotate_probe_seqs/vds_mrna_refseq_annotation.csv"]):
                    click.echo(f"Writing VDS annotations for REFSEQ to output/annotate_probe_seqs")
                    click.echo(f"Searching for matches")
                    df["refseq_hit"] = df["Sequence"].parallel_apply(lambda x: find_matches(x, refseq_seq_d))
                    click.echo(f"Search completed")
                    df["refseq_hit"] = df["refseq_hit"].apply(lambda x: tuple(x) if x else ())
                    df["n_refseq_hits"] = df["refseq_hit"].apply(len)
                    df["refseq_annotations"] = df["refseq_hit"].apply(lambda y: find_annotations(y,source="REFSEQ",seq_dict=refseq_seq_d))
                    df["refseq_majority_vote_annotation"] = df["refseq_annotations"].apply(majority_vote_with_tie)
                    df["refseq_majority_vote_annotation"] = df["refseq_majority_vote_annotation"].apply(
                        lambda x: tuple(x) if x else ())
                    annotations_with_lncrna = [a for a in df["refseq_majority_vote_annotation"].value_counts().index if
                                               "lncRNA" in a and "protein_coding" not in a and "protein_coding_CDS_not_defined" not in a]

                    click.echo(
                        f"Number of sequences annotated by refseq as lncRNA and not as protein_coding annotations:\n{df.refseq_majority_vote_annotation.value_counts().loc[annotations_with_lncrna]}")
                    df["refseq_hit"] = df["refseq_hit"].apply(lambda x: tuple(x) if x else ())
                    df["refseq_annotations"] = df["refseq_annotations"].apply(
                        lambda x: tuple(x) if x else ())

                    #df['GENCODE_gene_id'] = df['GENCODE_hit'].apply(extract_gene_id)
                    df.to_csv(out_path, index=False)
            else:
                print("Loading vds refseq annotations")
                lncrna_vds_not_in_gencode_df = pd.read_csv("output/annotate_probe_seqs/vds_lncrna_refseq_annotation.csv")
                mrnarna_vds_not_in_gencode_df = pd.read_csv("output/annotate_probe_seqs/vds_mrna_refseq_annotation.csv")

                seq_2probe_d_lncrna = dict(zip(lncrna_df["Sequence"], lncrna_df["ProbeName"]))
                seq_2probe_d_mrna = dict(zip(mrna_df["Sequence"], mrna_df["ProbeName"]))

                lncrna_vds_not_in_gencode_df["ProbeName"] = lncrna_vds_not_in_gencode_df.Sequence.map(seq_2probe_d_lncrna)
                mrnarna_vds_not_in_gencode_df["ProbeName"] = mrnarna_vds_not_in_gencode_df.Sequence.map(seq_2probe_d_mrna)



                lncrna_vds_not_in_gencode_df_orphan_probes = lncrna_vds_not_in_gencode_df[lncrna_vds_not_in_gencode_df.refseq_majority_vote_annotation=="()"]
                mrna_vds_not_in_gencode_df_orphan_probes = mrnarna_vds_not_in_gencode_df[mrnarna_vds_not_in_gencode_df.refseq_majority_vote_annotation=="()"]

                # Assuming df has columns "id" and "sequence"
                lncrna_vds_not_in_gencode_df_orphan_probes_seq_l = [SeqRecord(Seq(seq), id=seq_id, description="")
                           for seq_id, seq in zip(lncrna_vds_not_in_gencode_df_orphan_probes["ProbeName"], lncrna_vds_not_in_gencode_df_orphan_probes["Sequence"])]

                mrna_vds_not_in_gencode_df_orphan_probes_seq_l = [SeqRecord(Seq(seq), id=seq_id, description="")
                           for seq_id, seq in zip(mrna_vds_not_in_gencode_df_orphan_probes["ProbeName"], mrna_vds_not_in_gencode_df_orphan_probes["Sequence"])]

                SeqIO.write(lncrna_vds_not_in_gencode_df_orphan_probes_seq_l, "output/annotate_probe_seqs/vds_orphan_lncrna.fasta", "fasta")
                SeqIO.write(mrna_vds_not_in_gencode_df_orphan_probes_seq_l, "output/annotate_probe_seqs/vds_orphan_mrna.fasta", "fasta")
                exp_cols = [c for c in lncrna_df.columns if "normalized" in c]
                lncrna_exp_levels_plot = lncrna_df[["ProbeName"]+exp_cols]
                mrna_exp_levels_plot = mrna_df[["ProbeName"]+exp_cols]
                metadata = pd.read_csv("snakemake-workflows-rna-seq-star-deseq2-e103c1c/config/samples.tsv",sep="\t")
                mrna_exp_levels_plot_long = mrna_exp_levels_plot.reset_index(drop=True).melt(
                    id_vars='ProbeName',
                    var_name='sample_col',
                    value_name='expression'
                )
                lncrna_exp_levels_plot_long = lncrna_exp_levels_plot.reset_index(drop=True).melt(
                    id_vars='ProbeName',
                    var_name='sample_col',
                    value_name='expression'
                )
                mrna_exp_levels_plot_long['sample_name'] = mrna_exp_levels_plot_long['sample_col'].str.extract(r'\[(.*?)\]')
                lncrna_exp_levels_plot_long['sample_name'] = lncrna_exp_levels_plot_long['sample_col'].str.extract(r'\[(.*?)\]')

                # -----------------------------
                meta_long = metadata.melt(
                    id_vars='sample_name',
                    var_name='Condition',
                    value_name='treatment'
                )

                lncrna_exp_levels_plot_long_meta = lncrna_exp_levels_plot_long.merge(meta_long, on='sample_name')
                mrna_exp_levels_plot_long_meta = mrna_exp_levels_plot_long.merge(meta_long, on='sample_name')

                lncrna_exp_levels_plot_long_meta.loc[lncrna_exp_levels_plot_long_meta.ProbeName.isin(lncrna_vds_not_in_gencode_df_orphan_probes.ProbeName),"is_orphan"] = "yes"
                mrna_exp_levels_plot_long_meta.loc[mrna_exp_levels_plot_long_meta.ProbeName.isin(mrna_vds_not_in_gencode_df_orphan_probes.ProbeName),"is_orphan"] = "yes"
                lncrna_exp_levels_plot_long_meta.loc[~lncrna_exp_levels_plot_long_meta.ProbeName.isin(lncrna_vds_not_in_gencode_df_orphan_probes.ProbeName),"is_orphan"] = "no"
                mrna_exp_levels_plot_long_meta.loc[~mrna_exp_levels_plot_long_meta.ProbeName.isin(mrna_vds_not_in_gencode_df_orphan_probes.ProbeName),"is_orphan"] = "no"

                lncrna_exp_levels_plot_long_meta_avg_expr = lncrna_exp_levels_plot_long_meta.groupby(['Condition', 'is_orphan', 'sample_name'])['expression'].mean().reset_index()
                mrna_exp_levels_plot_long_meta_avg_expr = mrna_exp_levels_plot_long_meta.groupby(['Condition', 'is_orphan', 'sample_name'])['expression'].mean().reset_index()

                def plot_histograms(df, title):
                    conditions = df['Condition'].unique()
                    gene_groups = df['is_orphan'].unique()

                    plt.figure(figsize=(14, 8))

                    for i, condition in enumerate(conditions):
                        plt.subplot(2, (len(conditions) + 1) // 2, i + 1)
                        for group in gene_groups:
                            subset = df[(df['Condition'] == condition) & (df['is_orphan'] == group)]
                            plt.hist(subset['expression'], bins=20, alpha=0.6, label=group)
                        plt.title(condition)
                        plt.xlabel('Expression')
                        plt.ylabel('Count')
                        plt.legend(title='Gene Group')

                    plt.suptitle(title, fontsize=16)
                    plt.tight_layout(rect=[0, 0, 1, 0.96])
                    plt.show()

                # Plot for lncRNA
                plot_histograms(lncrna_exp_levels_plot_long_meta_avg_expr,
                                'Distribution of Average Expression per Gene Group Across Conditions (lncRNA)')
if __name__ == "__main__":
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord, Seq
    import pandas as pd
    from pandarallel import pandarallel
    import os
    import datetime
    import ast
    from collections import Counter
    import numpy as np
    import re
    import seaborn as sns
    import matplotlib.pyplot as plt
    main()



