# Tutorial for Bowtie2, SAMtools, and Visualization with IGV 
# Date: Feb 13th

## Required software: 
1. Bowtie2 and SAMtools should already be downloaded
Open by:

``` bash
module load bowtie2 
module load bcftools
```

2. IGV needs to be downloaded onto your computer. Here are the links:
   - Desktop: [Download IGV](https://software.broadinstitute.org/software/igv/download)
   - Web-based: [Use IGV Web App](https://igv.org/app/)

---

### --- three dashes splits sections in markdown documents (md). For this to work, the file name should have ".md" at the end

### the following synthax allow students to copy text directly from the tutorial (must remove # to work)
### ``` bash
### code
### ```
---
# First steps:
Fork Student-Led-Tutorial-1 repository, this is where you can send your data after it is analyzed
#In order to get all files in your workspace use:
``` bash
  git clone repositorylinkexample.com
```
Log into your bridges account in gitbash and open your personal folder in ocean
/ocean/projects/agr250001p/your-psc-username

Download the refrence genome (FASTA)
``` bash
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz
```
then gunzip the file  # are we working on human data??
   - Following this link (https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/) will show other refrence genomes that can be used. We used the first one.Also download two sequencing reads. These can be searched for through https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186573 It shows publically avalable wgs-sew datasets from GEO.

If you were to chose your own:
- Scroll down from the link page and hit 'SRA run selector'
   The table shown includes GSM samples, giving information such as sample accessions (begin with gsm) and metadata like age, sex, and health condition.
   One sample for this tutorial should be a healthy donor (NC)
  Find the SRR acession according to your chosen sample and click on it
   Navigate to the FASTA/FASTQ download tab

Use the SRA toolkit to download paired-end-FASTQ files
``` bash
   module load sra-toolkit
```
``` bash
fastq-dump --split-files SRR11412215
```
Replace SRR11412215 with your SRR acession of your sample and do for both samples
Ours:
Healthy donor: SRR16574675
Stage 2A CRC donor: SRR16574651

You can chose your own, but we provided examples here for your conveinence.

SLURM:
Edit the slurm script to include your email, username, and SRR acession number in the appropriate spots like we have done in class before. The code to edit:
``` bash
vi sra-download.slurm
```
You will also need to change the output name for each dataset
Navigate to the directory where the script is located (if not already there):
``` bash
cd /path/to/your/script
```
Submit the job to the SLURM scheduler:
``` bash
sbatch sra-download.slurm
```

Once the job completes, check the output log file:
``` bash
cat/ocean/projects/agr250001p/your-username/download_sra.log
``` 
THE FOLLOWING IS COPIED FROM THE ASSIGNMENT REQUIREMENTS AND WILL NEED TO BE EDITIED. WE LIKELY NEED TO ADD MIRE DETAIL!

Index the reference genome using Bowtie2:
``` bash
bowtie2-build reference.fasta reference_index
``` 
replace refrence.fasta with the file name for the refrence data

I HAVE NOT TESTED ANYTHING BEYOND THIS POINT

# Part 2: Mapping Reads

 Align the sequencing reads to the reference genome using Bowtie2:
``` bash
bowtie2 -x reference_index -1 reads_1.fastq -2 reads_2.fastq -S aligned.sam
```
# Process the SAM file using SAMtools

1. Convert to BAM:
``` bash
samtools view -Sb aligned.sam > aligned.bam
```
2. Sort and index the BAM file:
``` bash
samtools sort aligned.bam -o sorted.bam
samtools index sorted.bam
``` 
3.Compute basic statistics:
``` bash
samtools flagstat sorted.bam
``` 
# Part 3: Visualization:

1. Load the reference genome and sorted BAM file in IGV (Desktop or Web-based).
2. Highlight:
 - Coverage across interesting regions (e.g., BRCA1).
 - Mismatches or SNPs in the alignments

From workling through the tutorial I found that the web version of IGV was easier for me to navigate, so maybe we encourage that during our tutorial 

NOTES:
ADD SLURM CODE TO THIS REPOSITORY
DO WE NEED TO ADD PARTS ON SUBMITTING A JOB TO THE REPOSITORY? i think thats needed in order to open the data in IGV??? It does not been to be submitted, just need to add the reference genome and the data file

