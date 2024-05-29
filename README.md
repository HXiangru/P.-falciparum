# Align Oxford Nanopore reads to the *P. falciparum* ASM276v2 assembly 

## Preparation
1. Download the fastq file containing a few thousand Oxford Nanopore reads
2. Download the *P. falciparum* ASM276v2 assembly from [NCBI](https://www.ncbi.nlm.nih.gov/assembly/?term=txid36329[Organism:noexp)

## Alignment
1. Install minimap2
```linux
git clone https://github.com/lh3/minimap2
cd minimap2 && make
```
2. Alignment
```linux
minimap2 -ax map-ont GCF_000002765.6_GCA_000002765_genomic.fna reads.fastq > alignment.sam
```

## Subsequent processing
1. Install Samtools
```linux
wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
tar -jxvf samtools-1.20.tar.bz2
```
2. Convert SAM file into sorted BAM file
```linux
samtools view -Sb alignment.sam > aln.bam
samtools sort -@ 2 -m 4G -O bam -o aln.sorted.bam aln.bam
```
3. Create index
```linux
samtools faidx GCF_000002765.6_GCA_000002765_genomic.fna
samtools index aln.sorted.bam
```

## Visualisation in R
1. Library packages and import the file
```r
library(Rsamtools)
library(ggplot2)
bam <- scanBam("aln.sorted.bam")
```
2. Convert data to df and filter unaligned reads
```r
tmp=as.data.frame(do.call(cbind,lapply(bam[[1]], as.character)))
tmp=tmp[tmp$flag!=4,] 
```
3. Extract and count the mapping quality
```r
names(bam[[1]])
mapq <- tmp$mapq
mapq_counts <- table(tmp$mapq)
mapq_df <- data.frame(mapq = as.numeric(names(mapq_counts)), count = as.numeric(mapq_counts))
```
4. Creat a histogram of the mapping qualities
```r
ggplot(mapq_df, aes(x = mapq, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Mapping Quality Distribution", x = "Mapping Quality", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave('mapping_quality.png', width = 16, height = 14, units = 'cm', dpi = 300)
```






