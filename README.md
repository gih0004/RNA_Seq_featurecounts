# RNA_Seq_featurecounts# HTSEQ_counts
RNA seq pipeline designed for paired end reads using HTSEQ for producing files needed for visualizing transcriptomic data with DESEQ2. To use this script, put all fastq files within a single directory containing additionally your reference genome fna file and your reference gtf file.  
Fastq file names should be formated as follows: sample_1.fq.gz sample_2.fq.gz

This script requires three arguments and in specific order within the command line when submititng the script to HPC and shouuld look like: 
```ruby
$ htseq_counts <species name> <reference genome fna> <reference gtf> 
```
where species name can be either a common name for the species or scientific, should be one string in total.
where refernce genome fna is a fasta file that is used as the reference genome for sample alignment and indexing 
where reference gtf is a gene transfer format  file neccesary for generarting HTSEQ count file 

```ruby 
#!/bin/bash
#Make sure to change all <pwd> with the current working directory where you have all fastq raw reads and your gft and reference genome 
# Format for Fastq raw read file names : <sample>_1.fq.gz <sample>_2.fq.gz

#SBATCH --job-name=HTSEQ_counts
#SBATCH --ntasks=10
#SBATCH --partition=bigmem2
#SBATCH --export=ALL
#SBATCH --array=1-49
#SBATCH --time=96:00:00
#SBATCH --error=/<pwd>/error.err
#SBATCH --output=/<pwd>/output.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<easley associated email> 

These SBATCH commands are condtions specific for Easley HPC, the only important and universal command from this block is the shebang line. if using easley, to submit the script:
$ sbatch HTSEQ_counts <species> <reference genome> <reference gtf> 
```


Step one is preprocessing of rna raw reads. This block takes the fastq rna raw reads and does three things:  
1) runs fastqc on the raw reads   
2) trims of adapters using fastp and an automatic adapter recognition option  
3) runs fastqc on the filtered AND adapterterless rna reads created by fastp 
Adapters can also be specified within the fastp command

### STEP 1: Run fastqc and Filter Raw reads
```ruby
module load fastqc
mkdir ./FASTQC
fastqc *.fq.gz  -o ./FASTQC/  


for fq1 in ./*_1.fq.gz 
do 
    base="${fq1%_1.fq.gz}"
   
    fastp -i $fq1 -I ${base}_2.fq.gz -o ${base}.filtered.1.fq.gz -O ${base}.filtered.2.fq.gz --detect_adapter_for_pe --qualified_quality_phred 20 -h ${base}_fastp.html -j ${base}_fastp.json


done
mkdir ./FASTQC_filtered

cp *filtered* ./FASTQC_filtered
cd ./FASTQC_filtered

fastqc *.filtered.1*  
fastqc *.filtered.2* 
cd ..


mkdir FASTQC.html
mkdir FASTQC_filtered.html

cp ./FASTQC/*.html ./FASTQC.html
cp ./FASTQC_filtered/*.html ./FASTQC_filtered.html

```
Results of Step 1: You should have three new directories:  
a. FASTQC = has original fastqc from rna raw reads  
b. FASTQC_filtered  = has fastqc results from adapterless and quality filtered reads   
c. FASTQC_filtered.html = contains all html files for viewing fastqc results in browser from the filtered and adapterless reads 




### STEP 2: Run HISAT2

Hisat firs creates the index based on a fna file and then aligns reads to the created index.   
Usage:  
hisat2-build [options]* <reference_in> <ht2_base>  
<reference_in> A comma-separated list of FASTA files containing the reference sequences to be aligned to, or, if -c is specified, the sequences themselves = {reference_fa}   
<ht2-base> The basename of the index files to write. By default, hisat2-build writes files named NAME.1.ht2, NAME.2.ht2 where NAME is <ht2_base> = {species}  
```ruby 

module load hisat2
hisat2-build ${reference_fa} ${species}
echo "Finished index creation" >> Prrogress_ 
```
Note the Progress_ file, its purpose is to document where within the pipeline is the HPC at that time. To view you can do `less Progress_` on the comand line. This is the same for the execution of the rest of the script. 


### Step 3: HISAT2 alignment step 
After creating indices from genome, we can then run alignment of the samples against the indices recently created:  
```ruby
touch "hisat_alginment_process"
for fq1 in ./*.filtered.1.fq.gz 
do

    base="${fq1%.filtered.1.fq.gz}"

   #  gzip -t R2.fq.gz && echo ok || echo bad # checks tgat the file is good, would add as sanity check 

    hisat2 -p 8 -q -x ${species} -1 "$fq1" -2 "${base}.filtered.2.fq.gz" -S "${base}_aligned.sam" --summary-file "${base}_summary.txt"
    echo "Finished hisat2 alignment ${base}" >> hisat_alignment_process
done
echo "finished hisat alignments" >> Progress_

#-dta is for downstream aplications such as Stringtie
#-p is for processors being used
#-S is for the output sam file
#-x is for the indices built and being used for alignment
```
Note that a hisat_alignment_process fike is made, alignment is recursive in hisat and time consuming, hence a new file is made just to monitor the progress of the alignment step. 
Again, once the step is completed, the Progress_ file is updated and it states the progress of the pipeline  
  
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

### Step 4: converting SAM files to BAM files
#To do anything meaningful with alignment data you must swith from SAM to its binary counterpart BAM file. This binary format is much easier for computer programs such as StringTie to work with.

#Basic usage: 
#$ samtools <command> [options] Samtools has a vast amount of commands, we will use the sort command to sort our alignment files 
#-o gives the output file name
```ruby
module load samtools
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


### STEP 5: Generate featurecounts 

```ruby 
module load python
for file in *.bam;
  do base=${file%.bam};
#echo $tag
htseq-count -i gene_id -f bam -s no -r pos ${base}.bam ${reference_gtf} > ${base}_HTSEQ_counts
    echo "Progress on ${base}" >> Progress_ 
done
```

Tis is to fine tmp files that have nothing inside and elimate them. When I ran this code some temporary files from some of the samples werent elimated after completion. This command eliminates them only if empty  `find . -name '*bam.tmp*' -size 0 -exec rm {} + `
```ruby
mdkir HTSEQ_Counts
cp *HTSEQ_counts* ./HTSEQ_counts
```
  
With these 5 steps you should have all files for a holsitic understanding of your rna seq data. The following step would be visualitation, which should be done with DESEQ2 and which script is not provided at the time of creation of this repository. At the end of this pipeline you should have:  
1) Directories containing important files : FASTQC , FASTQC_filtered , FASTQC_filtered.html ; HTSEQ_counts
2) Progress_ file that indicates what was done and what was not
3) and error file, if specified within the sbatch options - this is specific for easley
4) alignment_summary file which is tab delimited and can be used to evaluate alignment of your dataset


