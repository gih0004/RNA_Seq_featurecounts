# RNA-Seq data processing using featurecounts # 
This repository describes a RNA seq pipeline designed for paired end reads using featurecounts module for producing gene count matrix neccesary for visualizing transcriptomic data with DESEQ2. 

To use this script, put all fastq files within a single directory containing additionally your reference genome fna file and your reference gtf file.  

Fastq file names should be formated as follows: sample_1.fq.gz sample_2.fq.gz

This script requires three arguments and in specific order within the command line when submititng the script to HPC and shouuld look like: 
```ruby
$ featurecounts_pipeline <species name> <reference genome fna> <reference gtf> 
```
where species name can be either a common name for the species or scientific, should be one string in total.
where refernce genome fna is a fasta file that is used as the reference genome for sample alignment and indexing 
where reference gtf file is a gene transfer format file from your species of interest, neccesary for generarting featurecounts output file 

```ruby 

#These  are the configurations I usewhen running this script on my HPC
#SBATCH --nodes=2
#SBATCH --mem=64G
#SBATCH --time=96:00:00 # change this with desired wall time limit

```


Step one is preprocessing of RNA raw reads. This block takes the fastq RNA raw reads and does three things:  
1) runs fastqc on the raw reads   
2) trims of adapters using fastp and an automatic adapter recognition option  
3) runs fastqc on the filtered and adapterterless RNA reads created by fastp 



## Step 1: RNA-seq read pre-proccesing ###
```ruby
fastqc *.fq.gz -o ./FASTQC/
# Run FastQC on all .fq.gz files (compressed FASTQ files) in the current directory and output the results to a directory named FASTQC.
# The '*.fq.gz' is a wildcard for selecting all files ending with .fq.gz.



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
```
Results of Step 1: You should have three new directories:  
a. FASTQC = has original fastqc from rna raw reads  
b. FASTQC_filtered  = has fastqc results from adapterless and quality filtered reads   
c. FASTQC_filtered.html = contains all html files for viewing fastqc results in browser from the filtered and adapterless reads 




## STEP 2: RUN HISAT2 hisat2-build

Hisat first creates the index based on a fna file and then aligns reads to the created index.   
This step is an alignment of your reads with genome indexes. It first creates the index based on a GTF/GFF file which is required to align reads to the created index of the reference genome.

Usage:  
hisat2-build [options]* <reference_in> <ht2_base>  
<reference_in> A comma-separated list of FASTA files containing the reference sequences to be aligned to, or, if -c is specified, the sequences themselves = {reference_fa}   
<ht2-base> The basename of the index files to write. By default, hisat2-build writes files named NAME.1.ht2, NAME.2.ht2 where NAME is <ht2_base> = {species}  

```ruby 
hisat2-build ${reference_fa} ${species}
echo "Finished creating index of ${species}" >> Pogress_pipeline
```
Note the Progress_ file, its purpose is to document where within the pipeline is the HPC at that time. To view you can do `less Progress_` on the comand line. This is the same for the execution of the rest of the script. 


## Step 3: HISAT2 alignment step 
After creating indices from genome, we can then run alignment of the samples against the indices recently created:  
```ruby
for fq1 in ./*.filtered.1.fq.gz #<>This should be changed into whichever string you have last in sample names common between all samples 
do

    base="${fq1%.filtered.1.fq.gz}"

    hisat2 -p 8 -q -x potato -1 "$fq1" -2 "${base}.filtered.2.fq.gz" -S "${base}_aligned.sam" --summary-file "${base}_summary.txt"
    echo "Finished hisat2 alignment ${base}" >> Progress_pipeline
done
```
Since alignment is recursive in hisat and time consuming, once a library sample has been correctly algined, it will stated on the Progress_pipeline file

  
This for loop block takes the alignment summaries created by hisat2 with your reference samples and creates a summary file that lets you know the sample, total reads, total aligned reads, aligned score (%) and reads aligned 0 times 
```ruby
touch "summary_file"
# This for loop block takes the alignment summaries created by hisat2 with your reference samples and creates a summary file that lets you know the sample, total reads, total aligned reads, aligned score (%) and reads aligned 0 times 
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

```
To view alignment summary, use `less summary_file` , this file can also be downloaded as it is a tab delimited table.  

## Step 4: Converting SAM files to BAM files
To do anything meaningful with alignment data you must swith from SAM to its binary counterpart BAM file. This binary format is much easier for computer programs to work with.

Basic usage: 
$ samtools <command> [options] Samtools has a vast amount of commands, we will use the sort command to sort our alignment files 
-o gives the output file name
```ruby
for i in ./*_aligned.sam  
do
    base="${i%_aligned.sam}"
    samtools sort  ${i}  -o ${base}.bam -O bam
    echo "Finished sam-> bam of ${base}" >> Progress_
done

```

To convert  gff to gtf file, which wont always be neccesary, use following command:
```ruby
module load gffread/0.9.8
gffread genomic.gff -T -o genomic.gtf
#this code block is not within the executable script 
```


## STEP 5: Run featurecounts to generate gene counts 

```ruby 
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
```

This is to organize the output from step 5
```ruby

# Create a directory to store feature counts results
mkdir featurecounts

# Create a directory to store sample alignment summaries
mkdir sample_alignments_summary

# Move all files ending with "_summary.txt" to the sample alignments summary directory
mv *_summary.txt ./sample_alignments_summary

# Copy all files ending with ".txt" or ".txt*" to the feature counts directory
cp *.txt* ./featurecounts


```
  
With these 5 steps you should have all files for a trasncriptomic analysis of your RNA-seq data. At the end of this pipeline you should have:  
1) Directories containing important files : FASTQC , FASTQC_filtered , FASTQC_filtered.html ; featurecounts
2) Progress_file that indicates what was done and what was not
3) Alignment_summary file which is tab delimited and can be used to evaluate alignment of your dataset
4) Genecount for samples will be in featurecounts directory

