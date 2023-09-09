#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition normal
#SBATCH --time=2:00:00
#SBATCH --mem=140G
#SBATCH --job-name=metagenomics_pipeline
#SBATCH --error ./log/%x_%j.err.txt
#SBATCH --output ./log/%x_%j.out.txt

## Parameters ############################################################################
R1_base=`basename $1 .fastq.gz`
R2_base=`basename $2 .fastq.gz`

READ_LENGTH=75

TMP_DIR=./tmp/${R1_base}
OUT_DIR=./output/${R1_base}
mkdir tmp
mkdir output
mkdir "$TMP_DIR"
mkdir "$OUT_DIR"
mkdir "$OUT_DIR"/alignment
mkdir "$OUT_DIR"/fastqc
mkdir "$OUT_DIR"/arriba
mkdir "$OUT_DIR"/kraken
mkdir "$OUT_DIR"/bracken
mkdir "$OUT_DIR"/bracken/report
mkdir "$OUT_DIR"/centrifuge
mkdir "$OUT_DIR"/scripts

# Copy markdown scripts in output folder
cp ./scripts/metagenome_report.Rmd "$OUT_DIR"/scripts/metagenome_report.Rmd
cp ./scripts/metagenome_report.R "$OUT_DIR"/scripts/metagenome_report.R

# Bowtie index (human reference)
BOWTIE_INDEX=/path/to/bowtie_index/GRCh38_noalt_as

# Kraken databases
DB_VIRAL=/path/to/k2_viral_20230314
DB_STANDARD=/path/to/k2_standard_20230314
DB_EUPATH=/path/to/k2_eupathdb48_20230407

# Centrifuge index
CENTRIFUGE_INDEX=/path/to/centrifuge/indices/nt

# Arriba & STAR
BASE_DIR=/path/to/arriba_v2.4.0
STAR_INDEX_DIR=$BASE_DIR/STAR_index_GRCh37viral_GENCODE19
ANNOTATION_GTF=$BASE_DIR/GENCODE19.gtf
ASSEMBLY_FA=$BASE_DIR/GRCh37viral.fa
BLACKLIST_TSV=$BASE_DIR/database/blacklist_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz
KNOWN_FUSIONS_TSV=$BASE_DIR/database/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz
TAGS_TSV="$KNOWN_FUSIONS_TSV"
PROTEIN_DOMAINS_GFF3=$BASE_DIR/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3

# Threads
THREADS_KRAKEN=8
THREADS=24


##########################################################################################
## TRIM ##################################################################################
##########################################################################################
if [ -f ""$TMP_DIR"/${R1_base}.trimmed.fq.gz" ]; then

    echo "Trimming for ${R1_base} already performed."

else

fqtrim -l 30 -D -p $THREADS -q 20 -P33 \
	-r "$OUT_DIR"/alignment/fqtrim_report.txt \
	-A -o trimmed.fq.gz --outdir "$TMP_DIR" "$1","$2"

fi


##########################################################################################
## FASTQC ################################################################################
##########################################################################################
if [ -f "${OUT_DIR}/fastqc/${R1_base}_fastqc.html" ]; then

    echo "Fastqc for ${R1_base} already performed."

else

fastqc --extract -t $THREADS -o "$OUT_DIR"/fastqc "$1" "$2"
fastqc --extract -t $THREADS -o "$OUT_DIR"/fastqc "$TMP_DIR"/*.trimmed.fq.gz

fi


##########################################################################################
## Host removal ##########################################################################
##########################################################################################
if [ -f "${TMP_DIR}/${R1_base}_nonhuman_reads.fastq.gz" ]; then

    echo "Host removal for ${R1_base} already performed."

else

bowtie2 -x $BOWTIE_INDEX -p $THREADS -1 "$TMP_DIR"/${R1_base}.trimmed.fq.gz \
	-2 "$TMP_DIR"/${R2_base}.trimmed.fq.gz \
        --un-conc-gz "$TMP_DIR"/${R1_base}_nonhuman_reads \
        -S "$TMP_DIR"/${R1_base}_host_alignment.sam

mv "$TMP_DIR"/${R1_base}_nonhuman_reads.1 "$TMP_DIR"/${R1_base}_nonhuman_reads.fastq.gz
mv "$TMP_DIR"/${R1_base}_nonhuman_reads.2 "$TMP_DIR"/${R2_base}_nonhuman_reads.fastq.gz

fi

echo $(zcat "$1"|wc -l)/4|bc > "$OUT_DIR"/alignment/total_reads.txt
echo $(zcat "$TMP_DIR"/${R1_base}_nonhuman_reads.fastq.gz|wc -l)/4|bc > "$OUT_DIR"/alignment/non_human_reads.txt


##########################################################################################
## Kraken2 ###############################################################################
##########################################################################################
if [ -f "${OUT_DIR}/kraken/${R1_base}.kraken.standard.txt" ]; then

   echo "Kraken output for ${R1_base} already exists."

else

kraken2 \
	--db $DB_STANDARD \
	--threads $THREADS_KRAKEN \
	--minimum-hit-groups 3 \
	--use-names \
	--report-minimizer-data \
	--output "$TMP_DIR"/${R1_base}.kraken \
	--report "$OUT_DIR"/kraken/${R1_base}.kraken.standard.txt \
	--paired "$TMP_DIR"/${R1_base}_nonhuman_reads.fastq.gz \
             "$TMP_DIR"/${R2_base}_nonhuman_reads.fastq.gz

kraken2 \
	--db $DB_VIRAL \
	--threads $THREADS_KRAKEN \
	--minimum-hit-groups 3 \
	--use-names \
	--report-minimizer-data \
	--output "$TMP_DIR"/${R1_base}.kraken \
	--report "$OUT_DIR"/kraken/${R1_base}.kraken.viral.txt \
	--paired "$TMP_DIR"/${R1_base}_nonhuman_reads.fastq.gz \
             "$TMP_DIR"/${R2_base}_nonhuman_reads.fastq.gz

kraken2 \
	--db $DB_EUPATH \
	--threads $THREADS_KRAKEN \
	--minimum-hit-groups 3 \
	--use-names \
	--report-minimizer-data \
	--output "$TMP_DIR"/${R1_base}.kraken \
	--report "$OUT_DIR"/kraken/${R1_base}.kraken.eupath.txt \
	--paired "$TMP_DIR"/${R1_base}_nonhuman_reads.fastq.gz \
             "$TMP_DIR"/${R2_base}_nonhuman_reads.fastq.gz

fi

if [ -f "${OUT_DIR}/bracken/${R1_base}.standard_genus.txt" ]; then

   echo "Bracken output for ${R1_base} already exists."

else

# Standard Database
bracken \
	-d $DB_STANDARD \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.standard.txt \
	-r $READ_LENGTH \
	-l G \
	-t 10 \
	-o "$OUT_DIR"/bracken/${R1_base}.standard_genus.txt \
	-w "$OUT_DIR"/bracken/report/${R1_base}.standard_genus.report.txt

bracken \
	-d $DB_STANDARD \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.standard.txt \
	-r $READ_LENGTH \
	-l S \
	-t 10 \
	-o "$OUT_DIR"/bracken/${R1_base}.standard_species.txt \
	-w "$OUT_DIR"/bracken/report/${R1_base}.standard_species.report.txt

# Viral Database
bracken \
	-d $DB_VIRAL \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.viral.txt \
	-r $READ_LENGTH \
	-l G \
	-t 10 \
	-o "$OUT_DIR"/bracken/${R1_base}.viral_genus.txt \
	-w "$OUT_DIR"/bracken/report/${R1_base}.viral_genus.report.txt

bracken \
	-d $DB_VIRAL \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.viral.txt \
	-r $READ_LENGTH \
	-l S \
	-t 10 \
	-o "$OUT_DIR"/bracken/${R1_base}.viral_species.txt \
	-w "$OUT_DIR"/bracken/report/${R1_base}.viral_species.report.txt

# EUPATHDB48 Database
bracken \
	-d $DB_EUPATH \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.eupath.txt \
	-r $READ_LENGTH \
	-l G \
	-t 10 \
	-o "$OUT_DIR"/bracken/${R1_base}.eupath_genus.txt \
	-w "$OUT_DIR"/bracken/report/${R1_base}.eupath_genus.report.txt

bracken \
	-d $DB_EUPATH \
	-i "$OUT_DIR"/kraken/${R1_base}.kraken.eupath.txt \
	-r $READ_LENGTH \
	-l S \
	-t 10 \
	-o "$OUT_DIR"/bracken/${R1_base}.eupath_species.txt \
	-w "$OUT_DIR"/bracken/report/${R1_base}.eupath_species.report.txt

fi


##########################################################################################
## Centrifuge ############################################################################
##########################################################################################
if [ -f ""$OUT_DIR"/centrifuge/${R1_base}.centrifuge.txt" ]; then

   echo "Centrifuge output for ${R1_base} already exists."

else

centrifuge -q \
-1 "$TMP_DIR"/${R1_base}_nonhuman_reads.fastq.gz \
-2 "$TMP_DIR"/${R2_base}_nonhuman_reads.fastq.gz \
--report-file "$OUT_DIR"/centrifuge/${R1_base}.centrifuge.txt \
--threads $THREADS \
-x $CENTRIFUGE_INDEX

fi


##########################################################################################
## Arriba ################################################################################
##########################################################################################
if [ -f ""$OUT_DIR"/arriba/${R1_base}virus_expression1.tsv" ]; then

   echo "Arriba analysis for ${R1_base} already exists."

else


STAR \
	--runThreadN "$THREADS" \
	--outFileNamePrefix "$TMP_DIR"/ \
	--genomeDir "$STAR_INDEX_DIR" \
	--genomeLoad NoSharedMemory \
	--readFilesIn "$TMP_DIR"/${R1_base}_nonhuman_reads.fastq.gz "$TMP_DIR"/${R2_base}_nonhuman_reads.fastq.gz \
	--readFilesCommand zcat \
	--outStd BAM_Unsorted \
	--outSAMtype BAM Unsorted \
	--outSAMunmapped Within \
	--outBAMcompression 0 \
	--outFilterMultimapNmax 50 \
	--peOverlapNbasesMin 10 \
	--alignSplicedMateMapLminOverLmate 0.5 \
	--alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 \
	--chimOutType WithinBAM HardClip \
	--chimJunctionOverhangMin 10 \
	--chimScoreDropMax 30 \
	--chimScoreJunctionNonGTAG 0 \
	--chimScoreSeparation 1 \
	--chimSegmentReadGapMax 3 \
	--chimMultimapNmax 50 |

tee "$TMP_DIR"/${R1_base}_Aligned.out.bam |

"$BASE_DIR/arriba" \
	-x /dev/stdin -I \
	-o "$OUT_DIR"/arriba/${R1_base}fusions.tsv \
	-a "$ASSEMBLY_FA" -g "$ANNOTATION_GTF" -b "$BLACKLIST_TSV" \
	-k "$KNOWN_FUSIONS_TSV" -t "$TAGS_TSV" -p "$PROTEIN_DOMAINS_GFF3"

samtools sort -@ "$THREADS" -m $((40000/THREADS))M -T tmp \
-O bam "$TMP_DIR"/${R1_base}_Aligned.out.bam > "$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam

rm -f "$TMP_DIR"/${R1_base}_Aligned.out.bam
samtools index "$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam

$BASE_DIR/scripts/quantify_virus_expression.sh "$TMP_DIR"/${R1_base}_Aligned.sortedByCoord.out.bam \
"$OUT_DIR"/arriba/${R1_base}virus_expression.tsv

fi


##########################################################################################
## De novo assembly ######################################################################
##########################################################################################
if [ -f ""$OUT_DIR"/assembly/${R1_base}.contigs.fa" ]; then

   echo "De novo assambly for ${R1_base} already exists."

else

megahit -1 "$TMP_DIR"/${R1_base}_nonhuman_reads.fastq.gz \
        -2 "$TMP_DIR"/${R2_base}_nonhuman_reads.fastq.gz \
        -o "$OUT_DIR"/assembly --min-contig-len 1000 \
        --out-prefix ${R1_base} -t $THREADS

fi

rm -R "$OUT_DIR"/assembly/intermediate_contigs/*


##########################################################################################
## Kraken of de novo assembled contigs ###################################################
##########################################################################################
if [ -f "$OUT_DIR"/assembly/${R1_base}.kraken.report.txt ]; then

   echo "Kraken output for ${R1_base} (de novo assembly) already exists."

else

kraken2 \
	--db $DB_STANDARD \
	--threads $THREADS_KRAKEN \
	--minimum-hit-groups 3 \
	--use-names \
	--report-minimizer-data \
	--output "$OUT_DIR"/assembly/${R1_base}.kraken.output.txt \
	--report "$OUT_DIR"/assembly/${R1_base}.kraken.report.txt \
	"$OUT_DIR"/assembly/${R1_base}.contigs.fa

fi


##########################################################################################
## Centrifuge of de novo assembled contigs ###############################################
##########################################################################################
if [ -f "$OUT_DIR"/assembly/${R1_base}.centrifuge.txt ]; then

   echo "Centrifuge output for ${R1_base} (de novo assembly) already exists."

else

centrifuge -f \
	--report-file "$OUT_DIR"/assembly/${R1_base}.centrifuge.report.txt \
	--threads $THREADS \
	-S "$OUT_DIR"/assembly/${R1_base}.centrifuge.txt \
	-x $CENTRIFUGE_INDEX \
	"$OUT_DIR"/assembly/${R1_base}.contigs.fa

fi


##########################################################################################
## Report ################################################################################
##########################################################################################
Rscript ./output/${R1_base}/scripts/metagenome_report.R ./output/${R1_base}/scripts

