#Creating bam files
setwd('D:\\Work\\RNA seq\\Final RNAseq')

#Install Rsubread from Bioconductor open the package.
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsubread", "limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))

install.packages("fastqcr")

library(Rsubread)
library(fastqcr)
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

#Making list of paired end read files
fastq1 <- list.files(path= './raw_data', pattern = "*.1.fq.gz$", full.names = TRUE)
fastq2 <- list.files(path= './raw_data', pattern = "*.2.fq.gz$", full.names = TRUE)

fastqc(fq.dir = base_dir, qc.dir = "./data", threads = 4)

#building index file using Rsubread. fasta.fa is genomic fasta from NZ_JAHDTB000000000.1
buildindex(basename="./index/cv017", reference="./genomic/fasta.fa")


#Aligning reads with the index
Rsubread::align(index = "./index/cv017",
      readfile1 = fastq1,
      readfile2 = fastq2,
      type = "rna",
      input_format = "gzFASTQ",
      output_format = "BAM",
      PE_orientation = "rf")


fc <- Rsubread::featureCounts(c("./raw_data/C601_1.fq.gz.subread.BAM",
                                "./raw_data/C602_1.fq.gz.subread.BAM",
                                "./raw_data/C603_1.fq.gz.subread.BAM",
                                "./raw_data/C801_1.fq.gz.subread.BAM",
                                "./raw_data/C802_1.fq.gz.subread.BAM",
                                "./raw_data/C803_1.fq.gz.subread.BAM",
                                "./raw_data/NO01_1.fq.gz.subread.BAM",
                                "./raw_data/NO02_1.fq.gz.subread.BAM",
                                "./raw_data/NO03_1.fq.gz.subread.BAM"
                                ),
            annot.ext = "./genomic/genomic.gtf",
            isPairedEnd = TRUE,
            isGTFAnnotationFile = TRUE,
            genome = "./genomic/fasta.fa",
            GTF.featureType = "gene")


#creating a file with counts of each gene instance in all BAM files
write.table(x = data.frame(fc$annotation[, c("GeneID", "Length")],
        fc$counts,
        stringsAsFactors = FALSE),
        file = "counts.txt",
        quote = FALSE,
        sep = "\t",
        row.names = FALSE)

#Importing counts data and sample information
seqdata <- read.delim("counts.txt", stringsAsFactors = FALSE)
sampleinfo <- read.delim("SampleInfo.txt", stringsAsFactors = TRUE)

#Removing the first two columns from counts table
countdata <- seqdata[,-(1:2)]

#Using gene_id as row name
rownames(countdata) <- seqdata[,1]

#Trimming column names to make it readable
colnames(countdata) <- substr(colnames(countdata),start=1,stop=4)

#Creating a DGEList object to store counts data
y <- DGEList(countdata)

#Creating groups of replicated
group <- paste(sampleinfo$Signal,sampleinfo$OD600,sep=".")
group <- factor(group)
y$samples$group <- group

#Filtering lowly expressed genes. Looking at counts per million and removing 
#all genes that have CPM of less than 1. Then we filtered out genes which has
#CPM of less than 3 in 3 or less sample.
CPM_sample <- cpm(countdata)

threshold <- CPM_sample > 1
keep <- rowSums(threshold) >= 3
logcounts <- cpm(y,log=TRUE)

#Normalining for composition bias
y <- calcNormFactors(y)

#Now creating desing matrix for limma-voom transformation
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#Performing voom transformation
v <- voom(y,design,plot = TRUE)

#Testing for differential equation. Fitting linear model using lmfit
fit <- lmFit(v)

#Making contrast matrix to compare between groups
contrast_matrix1 <- makeContrasts(C6VSNO=C6.2 - None.2,levels=design)
contrast_matrix2 <- makeContrasts(C8VSNO=C8.2 - None.2,levels=design)

#apply the contrasts matrix to the fit object from limma package
fitted_contrast1 <- contrasts.fit(fit, contrast_matrix1)
fitted_contrast2 <- contrasts.fit(fit, contrast_matrix2)

# empirical Bayes shrinkage on the variances, and estimates 
# moderated t-statistics and the associated p-values
fitted_contrast1 <- eBayes(fitted_contrast1)
fitted_contrast2 <- eBayes(fitted_contrast2)

#Calculating P-values based on logFC cutoff of 2
filtered_treated1 <- treat(fitted_contrast1,lfc=2)
filtered_treated2 <- treat(fitted_contrast2,lfc=2)

#Generating the table 
C6vsNone <- topTable(filtered_treated1,coef=1,number = 4530, sort.by="p")
C8vsNone <- topTable(filtered_treated2,coef=1,number = 4530, sort.by="p")

#Exporting the table
write.csv(C6vsNone, 'C6up.csv')
write.csv(C8vsNone, 'C8up.csv')

#We manually removed all genes that had Log2 fold change of 
# less than 2 and P>0.05 
