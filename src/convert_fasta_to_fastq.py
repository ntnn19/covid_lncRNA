import click

#records = []
#for rec in SeqIO.parse("output/annotate_probe_seqs/vds_orphan_lncrna.fasta", "fasta"):
#    rec.letter_annotations["phred_quality"] = [40] * len(rec.seq)  # Q40 = high quality
#    records.append(rec)

#count = SeqIO.write(records, "output/annotate_probe_seqs/vds_orphan_lncrna.fastq", "fastq")
#print(f"Converted {count} records")

#for rec in SeqIO.parse("output/annotate_probe_seqs/vds_orphan_mrna.fasta", "fasta"):
#    rec.letter_annotations["phred_quality"] = [40] * len(rec.seq)  # Q40 = high quality
#    records.append(rec)

#count = SeqIO.write(records, "output/annotate_probe_seqs/vds_orphan_mrna.fastq", "fastq")
#print(f"Converted {count} records")


@click.command()
@click.option("--input-xlsx", required=True, type=click.Path(exists=True), help="Input Excel file with columns ProbeName and Sequence")
@click.option("--output-dir", required=True, type=click.Path(), help="Output directory for FASTA/FASTQ files")
@click.option("--suffix", type=str, default=None, help="Optional suffix for output file")
@click.option("--skip-rows", type=int, default=None, help="last row of metadata")
def create_probe_files(input_xlsx, output_dir, suffix,skip_rows):
    """
    Convert probe sequences from an Excel file into deduplicated FASTA and FASTQ files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load Excel
    df = pd.read_excel(input_xlsx,skiprows=skip_rows)
    if "ProbeName" not in df.columns or "Sequence" not in df.columns:
        raise click.BadParameter("Excel must contain 'ProbeName' and 'Sequence' columns")

    # Deduplicate by Sequence
    df_unique = df.drop_duplicates(subset="Sequence")
    click.echo(f"Loaded {len(df)} probes, {len(df_unique)} unique sequences after deduplication")

    # If type_column is provided, split by type

    write_probe_files(df_unique, output_dir, suffix=suffix)

def write_probe_files(df, output_dir, suffix="all"):
    records = []
    for _, row in df.iterrows():
        rec = SeqRecord(
            Seq(row["Sequence"]),
            id=str(row["ProbeName"]),
            description=""
        )
        rec.letter_annotations["phred_quality"] = [40] * len(rec.seq)
        records.append(rec)

    fasta_file = output_dir / f"probes_{suffix}.fasta"
    fastq_file = output_dir / f"probes_{suffix}.fastq"

    count_fasta = SeqIO.write(records, fasta_file, "fasta")
    count_fastq = SeqIO.write(records, fastq_file, "fastq")

    click.echo(f"Written {count_fasta} records to {fasta_file}")
    click.echo(f"Written {count_fastq} records to {fastq_file}")

if __name__ == "__main__":
    import pandas as pd
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    from pathlib import Path
    create_probe_files()

