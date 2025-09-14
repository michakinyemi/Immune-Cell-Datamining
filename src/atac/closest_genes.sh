source "$(dirname "$0")/env.sh"

TYPE="top_down"

REGIONS_BED=$curDir"results/atac/$TYPE.bed"
GTF_FILE=$curDir"datasets/genomes/gencode.vM36.chr_patch_hapl_scaff.annotation.gtf"
WINDOW_DISTANCE=10000

GENE_BED=$(mktemp)
GENE_NAMES=$(mktemp)
GENE_COORDS=$(mktemp)

# Extract gene features from the GTF.
# We assume that gene features have $3=="gene" and that the attributes field ($9) contains "gene_name".
# If gene_name is not available, the script will try using "gene_id" instead.
awk '$3=="gene"' "$GTF_FILE" | \
  sed -n 's/.*gene_name "\([^"]\+\)".*/\1/p' > "$GENE_NAMES"

awk '$3=="gene"' "$GTF_FILE" | \
  awk '{ print $1 "\t" ($4 - 1) "\t" $5 }' > "$GENE_COORDS"

paste "$GENE_COORDS" "$GENE_NAMES" > "$GENE_BED"
rm "$GENE_COORDS" "$GENE_NAMES"

# Use bedtools window to find gene features within the specified distance (window) of the regions.
# This command will attach gene information from the gene BED file to each region in the regions BED.
bedtools window \
  -a <(sed 's/^/chr/' "$REGIONS_BED") \
  -b "$GENE_BED" \
  -w "$WINDOW_DISTANCE" > $curDir"results/atac/window_$TYPE.txt"

# From the output, extract the gene names.
# BEDTools window outputs first the columns of file A (your regions)
# then the columns of file B (the gene BED), where the 4th column of file B is the gene name.
# Since file A likely has 3 columns, the gene name will be in column 7.
cut -f7 temp_window.txt | sort | uniq > $curDir"results/atac/$TYPE.txt"

echo "Genes within ${WINDOW_DISTANCE} bp of regions in ${REGIONS_BED}:"
cat $curDir"results/atac/$TYPE.txt"

# Clean up temporary files
rm "$GENE_BED" temp_window.txt
