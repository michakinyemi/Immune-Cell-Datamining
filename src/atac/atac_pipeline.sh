
trimDir="$HOME/bin/Trimmomatic-0.39"
curDir="/home/michael/Research/Nguyen-Lab/ATAC-Seq-Analysis/"
resultsDir="$(pwd)/results/atac/"

sampleName=$1

if [ "$sampleName" == "Sample1" ]; then
  sampleDir="$(pwd)/datasets/Foxk1/FT-SA330767/"
  file1="QD1456-1.SE9236_SA330767_S32_L001_R1.fastq.gz"
  file2="QD1456-1.SE9236_SA330767_S32_L001_R2.fastq.gz"

elif [ "$sampleName" == "Sample2" ]; then
  sampleDir="$(pwd)/datasets/Foxk1/FT-SA330768/"
  file1="QD1456-2.SE9236_SA330768_S33_L001_R1.fastq.gz"
  file2="QD1456-2.SE9236_SA330768_S33_L001_R2.fastq.gz"

elif [ "$sampleName" == "Sample3" ]; then
  sampleDir="$(pwd)/datasets/Foxk1/FT-SA330769/"
  file1="QD1456-3.SE9236_SA330769_S34_L001_R1.fastq.gz"
  file2="QD1456-3.SE9236_SA330769_S34_L001_R2.fastq.gz"

elif [ "$sampleName" == "Sample4" ]; then
  sampleDir="$(pwd)/datasets/Foxk1/FT-SA330770/"
  file1="QD1456-4.SE9236_SA330770_S35_L001_R1.fastq.gz"
  file2="QD1456-4.SE9236_SA330770_S35_L001_R2.fastq.gz"

else
  echo "invalid sample, exiting..."
  exit
fi

genomeDir="$(pwd)/datasets/genomes/"
genomeFile="Mus_musculus.GRCm39.dna.primary_assembly.fa"

nextera_r1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
nextera_r2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"


# --== Genome Configuration ==--

# Isolate Chromosomes of Interest
# samtools faidx $genomeDir$genomeFile 5 8 10 11 17 19 X > $genomeDir"GRCm39_subset.fa"

# Index Genome
# bowtie2-build $genomeDir$genomeFile $genomeDir"GRCm39"
# OR
# bowtie2-build $genomeDir"GRCm39_subset.fa" $genomeDir"GRCm39_subset"

# --== Read Alignment & Coverage Prep ==--

# Adapter Trimming (via Trimmomatic)
java -jar $trimDir"/trimmomatic-0.39.jar" PE -threads $(nproc) -phred33 \
  $sampleDir$file1 $sampleDir$file2 \
  $sampleDir"trimmed_1.fastq" $sampleDir"unpaired_1.fastq" \
  $sampleDir"trimmed_2.fastq" $sampleDir"unpaired_2.fastq" \
  "ILLUMINACLIP:$(pwd)/src/nextera.fa:2:30:10" \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Align Reads to Genome
bowtie2 -x $genomeDir"GRCm39" --no-unal -1 $sampleDir"trimmed_1.fastq" -2 $sampleDir"trimmed_2.fastq" \
        -S $sampleDir$sampleName".sam" -p $(nproc)
rm $sampleDir"trimmed_1.fastq" $sampleDir"trimmed_2.fastq"


# --== Generating BAM File ==--

# Convert SAM to BAM
samtools view -bS -o $sampleDir$sampleName".bam" $sampleDir$sampleName".sam"
rm $sampleDir$sampleName".sam"

# Sort BAM File (queryname)
samtools sort -n $sampleDir$sampleName".bam" -o $sampleDir$sampleName".sorted.bam"
rm $sampleDir$sampleName".bam"

# Verify MC Tag
samtools fixmate -m $sampleDir$sampleName".sorted.bam" $sampleDir$sampleName".fixmate.bam"
rm $sampleDir$sampleName".sorted.bam"

# Sort BAM File (coord)
samtools sort $sampleDir$sampleName".fixmate.bam" -o $sampleDir$sampleName".sorted.bam"
rm $sampleDir$sampleName".fixmate.bam"

# Remove Duplicates
samtools markdup -r -s $sampleDir$sampleName".sorted.bam" $sampleDir$sampleName".bam"
rm $sampleDir$sampleName".sorted.bam"

# Filter Regions of Interest (Optional)
# samtools view -b $sampleDir"aligned_reads.sorted.bam" -L $genomeDir"target_regions.bed" > $sampleDir"aligned_reads.bam"
# rm $sampleDir"aligned_reads.dedup.bam"

# # Index BAM File
samtools index $sampleDir$sampleName".bam"


# --== Peak Calling ==--

# Call Peaks 
macs3 callpeak -t $sampleDir$sampleName".bam" -f BAMPE -g mm --nomodel --shift -100 --extsize 200 -B --SPMR -n $sampleDir$sampleName
cut -f1-3 $sampleDir$sampleName"_peaks.narrowPeak" > $sampleDir$sampleName".bed"
sed -i 's/^/chr/' $sampleDir$sampleName".bed"

# Convert to BigWig
bamCoverage -b $sampleDir$sampleName".bam" -o $sampleDir$sampleName".bw" --binSize 10 \
  --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreDuplicates --extendReads --centerReads -p max

