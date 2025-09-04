STAR_CONTAINER=$1
STAR_INDEX_DIR=$2
GENOME_FASTA_FILES=$3
GTF_FILE=$4
CPUS=$5
MEM=$6
mkdir -p $STAR_INDEX_DIR

singularity exec $STAR_CONTAINER STAR \
 --runMode genomeGenerate \
 --genomeDir $STAR_INDEX_DIR \
 --genomeFastaFiles $GENOME_FASTA_FILES \
 --sjdbGTFfile $GTF_FILE \
 --runThreadN $CPUS \
 $MEM \
 --genomeSAindexNbases 14 \
 --sjdbOverhang 59


# example cmd:
# bash src/genomeGenerate.sh /home/nadmin/singularity_containers/star/star.sif  output/align_to_hg38_star/index /home/nadmin/projects/groups/pfänder/jule/covid_lncRNA/input/fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /home/nadmin/projects/groups/pfänder/jule/covid_lncRNA/input/gtf/Homo_sapiens.GRCh38.115.gtf 24 30
#
