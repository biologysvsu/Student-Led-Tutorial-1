# Student-Led-Tutorial-1
# Task: Tutorial for Bowtie2, SAMtools, and Visualization with IGV (Web/Desktop) or Tablet
# Date: Feb 13th
## More information at: https://github.com/biologysvsu/biol443_syllabus/blob/main/tutorial.md

## **Objective**
Prepare a hands-on tutorial for your peers on:
1. Aligning short reads to a reference genome using Bowtie2.
2. Processing alignment files with SAMtools.
3. Visualizing the results with IGV (web-based or desktop) or Tablet.

Your tutorial should:
- Provide clear instructions with examples.
- Include explanations of each step and the purpose of the tools used.
- Engage the audience by making it relevant

---

## **Data to Use**

### **Input Files**
1. **Reference Genome**:
   - **GRCh38 (Human Genome Assembly)**:
     - Download the reference genome (FASTA) and annotation file (GTF) from Ensembl:
     - wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz
       - [FASTA file](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/)
   - Ensure the reference genome is indexed before alignment.

2. **Sequencing Reads**:
   - Use a publicly available wgs-seq dataset from GEO: [GSE186573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186573).
   - Explore the dataset, which includes samples for:
     - **CRC**: Colorectal cancer.
     - **STAD**: Stomach adenocarcinoma.
     - **NC**: Healthy donors (Normal Control).

---

### **Steps to Choose Samples**
1. Scroll down on the GEO page and run the **SRA Run Selector**.
2. Review the table of **GSM samples**, which provides information such as:
   - **Sample accessions**: Begin with `GSM`.
   - **Subject metadata**: Includes age, sex, and health condition.
3. Identify and select two samples for analysis:
   - At least one sample should be from a **healthy donor (NC)**.
4. Locate the **SRR accessions** corresponding to your selected samples.
   - Note: The SRR accession is required to download sequencing data.

---

### **Downloading Sequencing Reads**
1. Click on the SRR accession for your chosen sample.
   - The page may take a few minutes to load. Refresh if necessary.
2. Navigate to the **FASTA/FASTQ download tab**.
   - Most SRR files require the **SRA Toolkit** for download.
3. Use the SRA Toolkit to download paired-end FASTQ files:
   ```bash
   # Example: Download paired-end reads for your chosen SRR accession
   fastq-dump --split-files SRR11412215
- Replace SRR11412215 with the actual SRR accession for your sample.
4. Repeat this process for all selected samples.
---

## **Software Requirements**
1. **Bowtie2**: Align reads to the human reference genome.
   - [Bowtie2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
2. **SAMtools**: Convert, sort, and index SAM/BAM files.
   - [SAMtools Documentation](http://www.htslib.org/doc/)
3. **IGV**:
   - **Desktop**: [Download IGV](https://software.broadinstitute.org/software/igv/download)
   - **Web-based**: [Use IGV Web App](https://igv.org/app/)

---

## **Tasks and Deliverables**
### **Part 1: Data Preparation**
1. Access the provided links and download:
   - Reference genome.
   - Sample FASTQ files
2. Index the reference genome using Bowtie2:
   ```bash
   bowtie2-build reference.fasta reference_index
### **Part 2: Mapping Reads**

1. **Align the sequencing reads to the reference genome using Bowtie2**:
   ```bash
   bowtie2 -x reference_index -1 reads_1.fastq -2 reads_2.fastq -S aligned.sam

### **Process the SAM file using SAMtools**

1. **Convert to BAM**:
   ```bash
   samtools view -Sb aligned.sam > aligned.bam

2. **Sort and index the BAM file**:  

   ```bash
   samtools sort aligned.bam -o sorted.bam
   samtools index sorted.bam

3. **Compute basic statistics**:
   ```bash
   samtools flagstat sorted.bam

### **Part 4: Visualization**:
1. Load the reference genome and sorted BAM file in IGV (Desktop or Web-based).
2. Highlight:
- Coverage across interesting regions (e.g., BRCA1).
- Mismatches or SNPs in the alignments.
