# Student-Led-Tutorial-1
# Task: Tutorial for Bowtie2, SAMtools, and Visualization with IGV (Web/Desktop) or Tablet
# Date: Feb 13th
## More informatation at: https://github.com/biologysvsu/biol443_syllabus/blob/main/tutorial.md

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
- **Input File**: Human genome data (subset to keep it manageable).
  - Reference genome: **GRCh38 (Human Genome Assembly)**.
  - Sequencing reads: A publicly available dataset from the 1000 Genomes Project (e.g., chr22).
    - Download link for reads: [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/)
    - Use FASTQ files for paired-end reads.
  - Download link for reference: [Ensembl Reference Genomes](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/).

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
   - Sample FASTQ files (subset of chr22 or another interesting region).
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
