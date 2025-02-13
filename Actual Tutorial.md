# Tutorial for Bowtie2, SAMtools, and Visualization with IGV 
# Date: Feb 13th

## Required software: 
1. Bowtie2 and SAMtools should already be downloaded
Open by:

``` bash
module load bowtie2 
module load samtools
```

2. IGV needs to be downloaded onto your computer. Here are the links:
   - Desktop: [Download IGV](https://software.broadinstitute.org/software/igv/download)
   - Web-based: [Use IGV Web App](https://igv.org/app/)

3. Windows users may need WinSCP to upload their files to IGV.
   - [Download WinSCP](https://winscp.net/eng/index.php) 

---

### --- three dashes splits sections in markdown documents (md). For this to work, the file name should have ".md" at the end

### the following synthax allow students to copy text directly from the tutorial (must remove # to work)
### ``` bash
### code
### ```
---
# Part 1:
Fork Student-Led-Tutorial-1 repository, this is where you can send your data after it is analyzed
#In order to get all files in your workspace use:
``` bash
  git clone repositorylinkexample.com
```
Log into your bridges account in gitbash and open your personal folder in ocean
``` bash
ssh your-username@bridges2.psc.edu
```
To avoid a slurm script we are going to use interact:
``` bash
interact -t 3:00:00 –-ntasks-per-node=16  --mem=31G
```
The reference genome has been provided for you guys for effeciency, but if you wer to download your own reference genome you would follow the steps below:

**Here are the steps if you were to do this yourself. We have already done this, so you can skip ahead**
``` bash
cd /ocean/projects/agr250001p/your-psc-username
```
Download the refrence genome (FASTA)
``` bash
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz
```
then gunzip the file
``` bash
gunzip Homo_sapiens.GRCh38.dna.alt.fa.gz
```
   - Following this link (https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/) will show other refrence genomes that can be used.
   - Also need to download two sequencing reads. These can be searched for through https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186573
   - It shows publically avalable wgs-sew datasets from GEO.
You can choose your own, but we provided examples here for your conveinence.

If you were to choose your own:
- Scroll down from the link page and hit 'SRA run selector'
   The table shown includes GSM samples, giving information such as sample accessions (begin with gsm) and metadata like age, sex, and health condition.
   One sample for this tutorial should be a healthy donor (NC)
  Find the SRR acession according to your chosen sample and click on it
   Navigate to the FASTA/FASTQ download tab

Use the SRA toolkit to download paired-end-FASTQ files
``` bash
   module load sra-toolkit
```
Healthy Example:
``` bash
fastq-dump --split-files SRR16574675
```
Stage 2A CRC Donor Example:
``` bash
fastq-dump --split-files SRR16574651
```

**HERE IS WHERE WE WILL PICK BACK UP**

# Copying the necessary files to your computer 
1. You should still be in the interact verion of the HPC, but if not then:
   ``` bash
   exit
   ``` 
2. Now enter this to copy the necessary files:
   ``` bash
    cp -r /ocean/projects/agr250001p/shared/tutorial-data/tutorial_1_data/ .
    ```
3. Ensure the files were properly copied over:
   ``` bash
   ls
   ```
   

# Mapping Reads

  Process a subset of each dataset with Bowtie2. Repeat for both healthy (NC) and cancer (CRC) donors:
``` bash
seqtk sample -s100 SRR16574675_1.fastq 0.5 > NC_1.fastq
seqtk sample -s100 SRR16574675_2.fastq 0.5 > NC_2.fastq
```
``` bash
seqtk sample -s100 SRR16574651_1.fastq 0.5 > CRC_1.fastq
seqtk sample -s100 SRR16574651_2.fastq 0.5 > CRC_2.fastq
```

 Align the sequencing reads to the reference genome using Bowtie2:
``` bash
module load bowtie
``` 
 
- For our the helathy genome (NC)
``` bash
bowtie2 --very-fast-local -p 16 -x human_bowtie_reference -1 NC_1.fastq -2 NC_2.fastq -S NC.sam
```
- For the cancer genome (CRC)
``` bash
bowtie2 --very-fast-local -p 16 -x human_bowtie_reference -1 CRC_1.fastq -2 CRC_2.fastq -S CRC.sam
```
  
# Process the SAM file using SAMtools
We now need to convert our .sam files to .bam files so that they can be used by the IGV software.

1. Convert to BAM:
``` bash
samtools view -Sb NC.sam > NC.bam
```
``` bash
samtools view -Sb CRC.sam > CRC.bam
```
2. Sort and index the BAM file for both healthy and cancer donors:
``` bash
samtools sort CRC.bam -o sorted-CRC.bam
samtools index sorted-CRC.bam
```
``` bash
samtools sort NC.bam -o sorted-NC.bam
samtools index sorted-NC.bam
``` 
3.Compute basic statistics:
``` bash
samtools flagstat sorted-CRC.bam
samtools flagstat sorted-NC.bam
``` 
# Part 3: Visualization:

- Files must be transferred to your computer to be visualized:
     1. If you are not in the interact computer then use:
      ``` bash
        exit
        ```
2. To ensure that you are able to find your downloads use:
   ``` bash
   pwd
   ```
   3. Please note where you will be locally dowloading these files.

 We now need to downlaod the files to our computers for them to be used in IGV.
 1. ``` bash
scp your-username@bridges2.psc.edu:/ocean/projects/agr250001p/your-username/sorted-sample.bam .
``` 
``` bash
scp your-username@bridges2.psc.edu:/ocean/projects/agr250001p/your-username/sorted-sample.bam.bai .
```
Here you will have to cahnge the "sorted-sample" to the sorted-NC and sorted-CRC.

- Windows: Use your Bridges2 credentials to log in to WinSCP.
File protocol: SCP
Host name: bridges2.psc.edu
- If not in our shared directory, click on "Find Files" and in the "Search in" box, type /ocean/projects/agr250001p/your-user
- Find your sorted sample files ending in .bam and .bam.bai and download them.

1. Load the reference genome and sorted BAM file in IGV (Desktop or Web-based).
2. Highlight:
 - Coverage across interesting regions (e.g., BRCA1).
 - Mismatches or SNPs in the alignments



NOTES:
ADD SLURM CODE TO THIS REPOSITORY
DO WE NEED TO ADD PARTS ON SUBMITTING A JOB TO THE REPOSITORY? i think thats needed in order to open the data in IGV??? It does not been to be submitted, just need to add the reference genome and the data file

