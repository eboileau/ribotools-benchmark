#!/usr/bin/env Rscript --vanilla

source(file.path(Sys.getenv("SCRIPT_DIR"), "utils.R"))
opts <- parse_args()

library(DESeq2)

ribo <- as.matrix(read.csv(opts$ribo, row.names = 1, check.names = FALSE))
rna  <- as.matrix(read.csv(opts$rna, row.names = 1, check.names = FALSE))
samples <- read.csv(opts$samples)

alpha <- as.numeric(opts$alpha)

# Wrangle to match published protocol
merge <- cbind(ribo, rna)
coldata <- samples[,c("sampleName", "condition", "assay")]
colnames(coldata) <- c("SampleID", "Condition", "SeqType")
coldata <- as.data.frame(apply(coldata,2,as.factor))
coldata$SeqType <- toupper(coldata$SeqType)

# no batch
batch <- 0

# Taken from https://github.com/SGDDNB/translational_regulation/commit/60b61ca2061546999de1598f4d280c86383031b9
# "patch" script to handle any conditions - default lexicographic ordering
conditions <- as.character(levels(as.factor((coldata$Condition))))
# -------------------------------------------------------------------------------------------------------------

if(batch == 1){
  ddsMat <- DESeqDataSetFromMatrix(countData = merge,
                                   colData = coldata, design =~ Batch + Condition + SeqType + Condition:SeqType)
}else if(batch == 0){
  ddsMat <- DESeqDataSetFromMatrix(countData = merge,
                                   colData = coldata, design =~ Condition + SeqType + Condition:SeqType)
}else{
  stop("Batch presence should be indicated by 0 or 1 only", call.=FALSE)
}

ddsMat$SeqType = relevel(ddsMat$SeqType,"RNA")
ddsMat <- DESeq(ddsMat)

# Condition2.SeqTypeRibo.seq means Changes in Ribo-seq levels in Condition2 vs
# Condition1 accounting for changes in RNA-seq levels in Condition2 vs Condition1
res <- results(ddsMat, contrast=list(paste0("Condition", conditions[2], ".SeqTypeRIBO")))

# DESeq2 object with batch for Ribo-seq
ind = which(coldata$SeqType == "RIBO")
coldata_ribo = coldata[ind,]

if(batch == 1){
  ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                        colData = coldata_ribo, design =~ Condition + Batch)
  }else if(batch ==0){
  ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                        colData = coldata_ribo, design =~ Condition)
}

ddsMat_ribo <- DESeq(ddsMat_ribo)
res_ribo <- results(ddsMat_ribo, contrast=c("Condition", conditions[2], conditions[1]))
res_ribo <- lfcShrink(ddsMat_ribo, coef=2,res=res_ribo,type="apeglm")

# DESeq2 object with batch for RNA-seq
ind = which(coldata$SeqType == "RNA")
coldata_rna = coldata[ind,]

if(batch == 1){
  ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                       colData = coldata_rna, design =~ Condition + Batch)
}else if(batch ==0){
  ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                       colData = coldata_rna, design =~ Condition)
}

ddsMat_rna <- DESeq(ddsMat_rna)
res_rna <- results(ddsMat_rna, contrast=c("Condition", conditions[2], conditions[1]))
res_rna <- lfcShrink(ddsMat_rna, coef=2,type="apeglm",res=res_rna)

# Classes of genes
forwarded = rownames(res)[which(res$padj > alpha & res_ribo$padj < alpha & res_rna$padj < alpha)]
exclusive = rownames(res)[which(res$padj < alpha & res_ribo$padj < alpha & res_rna$padj > alpha)]
both = which(res$padj < alpha & res_ribo$padj < alpha & res_rna$padj < alpha)
intensified = rownames(res)[both[which(res[both,2]*res_rna[both,2] > 0)]]
buffered = rownames(res)[both[which(res[both,2]*res_rna[both,2] < 0)]]
buffered = c(rownames(res)[which(res$padj < alpha & res_ribo$padj > alpha & res_rna$padj < alpha)],buffered)

# -------------------------------------------------------------------------------------------------------------

# unify output
res$effect <- res$log2FoldChange
res$score <- -log10(res$padj + 1e-300)

res$effect[rownames(res) %in% buffered] <- res_rna$log2FoldChange[rownames(res) %in% buffered] 
res$score[rownames(res) %in% buffered] <- -log10(res_rna$padj[rownames(res) %in% buffered] + 1e-300)
res$effect[rownames(res) %in% forwarded] <- res_rna$log2FoldChange[rownames(res) %in% forwarded] 
res$score[rownames(res) %in% forwarded] <- -log10(res_rna$padj[rownames(res) %in% forwarded] + 1e-300)

res$class <- "background"
res$class[rownames(res) %in% forwarded] <- "abundance"
res$class[rownames(res) %in% intensified] <- "abundance"
res$class[rownames(res) %in% exclusive] <- "translation"
res$class[rownames(res) %in% buffered] <- "buffering"

res$identifier <- rownames(res)
res <- res[,c("identifier", "class", "effect", "score")]

write.csv(res, opts$out, row.names = FALSE, quote = FALSE)
