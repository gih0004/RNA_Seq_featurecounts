#!/bin/bash
#SBATCH --nodes=2
#SBATCH --mem=64G


# Check if the number of arguments is less than three
if [ "$#" -lt 3 ]; then
  echo "Error: Insufficient arguments. Please provide three arguments."
  exit 1
fi

# Assigning the command-line arguments to variables
species=$1  #can be latin species name or common name
reference_fa=$2 #reference genome for alignment 
reference_gtf=$3 #reference gtf file 


# Using the variables
echo "Species: ${species}"
echo "reference genome file: ${reference_fa}"
echo "referencegtf file: ${reference_gtf}"


### Module Dependcy####
module load FastQC/0.11.9-Java-11
module load HISAT2/2.2.1-foss-2020a
module load StringTie/2.2.1-GCC-11.2.0
module load SAMtools/1.16.1-GCC-11.3.0 
module load fastp/0.23.2-GCC-11.2.0
module load Subread/2.0.3-GCC-11.2.0 


mkdir ./FASTQC
mkdir ./FASTQC_filtered
mkdir FASTQC.html
mkdir FASTQC_filtered.html
mkdir featurecounts
mkdir sample_alignmnets_summary

#Creates File that is continuosly updated on pipeline progress
touch Progress_pipeline


### Step 1: RNA-seq read pre-proccesing ###
# FASTQC is a tool used for quality checking of raw sequencing data (in this case, RNA reads). It provides a report on the quality of the data which can help identify if further processing (like trimming) is needed.

fastqc *.fq.gz -o ./FASTQC/
# Run FastQC on all .fq.gz files (compressed FASTQ files) in the current directory and output the results to a directory named FASTQC.
# The '*.fq.gz' is a wildcard for selecting all files ending with .fq.gz.



# Fastp is a tool for fast preprocessing of sequencing data. It can trim adapters and filter reads based on quality among other functionalities.

# Loop through each file that matches the pattern *_1.fq.gz (typically one of paired-end reads from sequencing data), perform adapter trimming and quality filtering, and output filtered data.
#Adapter trimming is automatically detected by fastp with ouption 'detect_adapter_for_pe'   (paired end)
#Quality is set for phred value of 20, adjust if desired 

for fq1 in ./*_1.fq.gz
do 
    base="${fq1%_1.fq.gz}"
   
    fastp -i $fq1 -I ${base}_2.fq.gz -o ${base}.filtered.1.fq.gz -O ${base}.filtered.2.fq.gz --detect_adapter_for_pe --qualified_quality_phred 20 -h ${base}_fastp.html -j ${base}_fastp.json


done


# Organize output: Copy of the filtered data is moved to this directory.
cp *filtered* ./FASTQC_filtered
cd ./FASTQC_filtered



# Run FastQC on the filtered data to assess quality after preprocessing. YOu shpould notice no adapters detected from the fastqc after fastp step
fastqc *.filtered.1*  
fastqc *.filtered.2* 
cd ..


# Copies HTML reports from FastQC and moves the reports to these directories for better organization.
cp ./FASTQC/*.html ./FASTQC.html
cp ./FASTQC_filtered/*.html ./FASTQC_filtered.html


###STEP 2: RUN HISAT2 hisat2-build ####
#This step is an alignment of your reads with genome indexes. It first creates the index based on a GFF file and then aligns reads to the created index

#Usage:
#hisat2-build [options]* <reference_in> <ht2_base>
#<reference_in> A comma-separated list of FASTA files containing the reference sequences to be aligned to, or, if -c is specified, the sequences themselves
#<ht2-base> The basename of the index files to write. By default, hisat2-build writes files named NAME.1.ht2, NAME.2.ht2 where NAME is <ht2_base>
#--large-index  Force hisat2-build to build a large index, even if the reference is less than ~ 4 billion nucleotides long.
#-f The reference input files (specified as <reference_in>) are FASTA files (usually having extension .fa, .mfa, .fna or similar).


hisat2-build ${reference_fa} ${species}



#STEP 3: Run hisat2 alignment step 


# Run HISAT2 alignment
# -p 8: Use 8 processors for the alignment process
# -q: Input files are in FASTQ format
# -x potato: Specify the base name of the index for the reference genome
# -1 and -2: Specify the files for read 1 and read 2 of paired-end data
# -S: Specify the output file name for the SAM format
# --summary-file: Specify the file to write alignment summary statistics

for fq1 in ./*.filtered.1.fq.gz #<>This should be changed into whichever string you have last in sample names common between all samples 
do

    base="${fq1%.filtered.1.fq.gz}"

    hisat2 -p 8 -q -x potato -1 "$fq1" -2 "${base}.filtered.2.fq.gz" -S "${base}_aligned.sam" --summary-file "${base}_summary.txt"
    echo "Finished hisat2 alignment ${base}" >> Progress_pipeline
done



# This for loop block takes the alignment summaries created by hisat2 with your reference samples and creates a summary file that lets you know the sample, total reads, total aligned reads, aligned score (%) and reads aligned 0 times 
touch "summary_file"

echo -e "File\tTotal Reads\tAlignment Score (%)\tPairs Aligned Concordantly 0 times" > summary_file  
for file in ./*summary.txt; do
  # Extract file name
  base="${file%_summary.txt}"

  total_reads=$(head -1 "$file" | grep -o '^[0-9]\+')
  # Extract final alignment score (last line)
  alignment_score=$(tail -n 1 "$file" | grep -oP '\d+\.\d+(?=%)')
  # Extract pairs_aligned_0 value and remove leading space
  pairs_aligned_0=$(grep -n '^' "$file" | grep -E '^3:' | cut -d':' -f2- | sed 's/ aligned concordantly 0 times//g' | sed 's/^[[:space:]]*//')

  echo -e "${base#./}\t$total_reads\t$alignment_score\t$pairs_aligned_0" >> summary_file
done



# STEP 4: RUN SAMTOOLS

#To do anything meaningful with alignment data you must swith from SAM to its binary counterpart BAM file. This binary format is much easier for computer programs such as StringTie to work with.

#Basic usage: 
#$ samtools <command> [options] Samtools has a vast amount of commands, we will use the sort command to sort our alignment files 
#-o gives the output file name

 
for i in ./*_aligned.sam  
do
    base="${i%_aligned.sam}"
    samtools sort  ${i}  -o ${base}.bam -O bam
    
    echo "Finished SAM to BAM conversion of ${base}" >> Progress_pipeline
done



### Step 5: RUN FEATURECOUNTS
#options used for featurecounts in this base code:
#-F: specifies annotation type file (GTF)
#-p: is paired end
#-g: Specify the attribute type used to group features (eg. exons) into meta-features (eg. genes) when GTF annotation is provided. ‘gene id’ by default. This attribute type is usually the gene identifier.
#-O: allow multioverlap of reads 
#-t: specifies feature types to count, can be more than one
#-T: specifies number of threads
#-M: count multimapping reaads
#-a: Provide name of annotation file
#-o: Provide name for outputfile 


# Loop through all BAM files in the current directory
for file in *.bam
do 
    # Extract the base name of the file by removing the '.bam' suffix
    tag="${file%.bam}"
    
    # Run featureCounts to count reads that align to features in the provided GTF file
    featureCounts -p -O -T 2 -g ID -f -F GTF -t mRNA -M --countReadPairs  -a ${reference_gtf} -o ${tag}.txt $file 
    
    # Remove temporary files generated by featureCounts
    rm *temp
done


# Create a directory to store feature counts results
mkdir featurecounts

# Create a directory to store sample alignment summaries
mkdir sample_alignments_summary

# Move all files ending with "_summary.txt" to the sample alignments summary directory
mv *_summary.txt ./sample_alignments_summary

# Copy all files ending with ".txt" or ".txt*" to the feature counts directory
cp *.txt* ./featurecounts








