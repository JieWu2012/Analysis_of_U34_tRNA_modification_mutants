## DE analysis for 5p-seq data

# Read raw counts. 
p5_d=read.table("/Users/jwu/Downloads/all_count_rmdup.txt",header = T,row.names = 1)
head(p5_d)
colData <- data.frame (row.names = colnames(p5), condition = c("wt","wt","wt","h2","h2","h2","n2e6","n2e6","n2e6","n2e6h2","n2e6h2","n2e6h2"))
dds <- DESeqDataSetFromMatrix (countData = p5_d,colData = colData, design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 5, ]
dds$condition <- relevel (dds$condition, ref = "wt")
dds <- DESeq(dds)

# Generate PCA plot
vsd=vst(dds,blind = F)
pcaData=plotPCA(vsd, intgroup=c("condition"),returnData=T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group))+geom_point(size=0.1)+geom_jitter()+xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+scale_color_manual(values = color2)+theme_bw()+theme(text=element_text(family="Helvetica", size=12),panel.grid.minor = element_blank())

# For comparison between hel2 and wt
res=results(dds,contrast=c("condition","n2e6","wt"),alpha=0.05)
resultsNames(dds)
res <- lfcShrink(dds, coef=3,type = "apeglm",res=res)
head(res)
res <- res[order(res$padj), ]
res=data.frame(gene_name=rownames(res),res)
write.table(res,"/Users/jwu/Downloads/all_gene_5p_ncs2elp6_vs_wt_p0.05_shrink_rmdup.txt",row.names=F,sep="\t",col.names=T,quote=F)

# For comparison between hel2 and wt
res=results(dds,contrast=c("condition","h2","wt"),alpha=0.05)
res <- lfcShrink(dds, coef=2,type = "apeglm",res=res)
res <- res[order(res$padj), ]
res=data.frame(gene_name=rownames(res),res)
write.table(res,"/Users/jwu/Downloads/all_gene_5p_hel2_vs_wt_p0.05_shrink_rmdup.txt",row.names=F,sep="\t",col.names=T,quote=F)


# For comparison between n2e6h2 and n2e6
res=results(dds,contrast=c("condition","n2e6h2","n2e6"),alpha=0.05)
resultsNames(dds)
res <- lfcShrink(dds, coef=3,type = "apeglm",res=res)
res <- res[order(res$padj), ]
res=data.frame(gene_name=rownames(res),res)
write.table(res,"/Users/jwu/Downloads/all_gene_5p_ncs2elp6hel2_vs_ncs2elp6_p0.05_shrink_rmdup.txt",row.names=F,sep="\t",col.names=T,quote=F)

## Need to normalize to mRNA levels. 

# Read raw counts from RNA-seq
mRNA=read.table("/Users/jwu/Downloads/all_gene_mRNA_2batch_1.txt",header = 1)
p5_d=read.table("/Users/jwu/Downloads/all_count_rmdup.txt",header = T,row.names = 1)

# Select ncs2elp6 and wt for both data. 
p5_ds=p5_d[,c(1:3,7:9)]
mRNA_s=mRNA[,c(3,4,7,8)]

# Merge both data. 
all=merge(p5_s,mRNA_s,by=0,suffixes = c("_5p","_mRNA"))
all[is.na(all)] <- 0
rownames(all)=all$Row.names
all=all[,-1]

# Generate ColData. 
colData <- data.frame (row.names = colnames(all), condition = c("wt","wt","wt","mutant","mutant","mutant","mutant","mutant","wt","wt"),sampleType=c("5p","5p","5p","5p","5p","5p","mRNA","mRNA","mRNA","mRNA"))
dds <- DESeqDataSetFromMatrix (countData = all,colData = colData, design = ~ sampleType+condition+sampleType:condition)
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds$sampleType <- relevel (dds$sampleType, ref = "mRNA")
dds$condition <- relevel (dds$condition, ref = "wt")
dds <- DESeq(dds, test="LRT", reduced =  ~ sampleType+condition)
res <- results(dds,name="sampleType5p.conditionmutant",alpha=0.05)
res2 <- lfcShrink(dds, coef=4,type = "apeglm",res=res)
res2 <- res2[order(res2$padj), ]
res2=data.frame(gene_name=rownames(res2),res2)
write.table(res2,"/Users/jwu/Downloads/all_gene_5p_vs_mRNA_ncs2elp6_vs_wt_p0.05_shrink_rmdup.txt",row.names=F,sep="\t",col.names=T,quote=F)

# Similar analysis for ncs2elp6hel2 and wt
p5_ds=p5_d[,c(1:3,10:12)]
mRNA_s=mRNA[,c(5,6,7,8)]
all=merge(p5_ds,mRNA_s,by=0,suffixes = c("_5p","_mRNA"))
head(all)
all[is.na(all)] <- 0
rownames(all)=all$Row.names
all=all[,-1]
colData <- data.frame (row.names = colnames(all), condition = c("wt","wt","wt","mutant","mutant","mutant","mutant","mutant","wt","wt"),sampleType=c("5p","5p","5p","5p","5p","5p","mRNA","mRNA","mRNA","mRNA"))
dds <- DESeqDataSetFromMatrix (countData = all,colData = colData, design = ~ sampleType+condition+sampleType:condition)
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds$sampleType <- relevel (dds$sampleType, ref = "mRNA")
dds$condition <- relevel (dds$condition, ref = "wt")
dds <- DESeq(dds, test="LRT", reduced =  ~ sampleType+condition)
res <- results(dds,name="sampleType5p.conditionmutant",alpha=0.05)
res2 <- lfcShrink(dds, coef=4,type = "apeglm",res=res)
res2 <- res2[order(res2$padj), ]
res2=data.frame(gene_name=rownames(res2),res2)
write.table(res2,"/Users/jwu/Downloads/all_gene_5p_vs_mRNA_ncs2elp6hel2_vs_wt_p0.05_shrink_rmdup.txt",row.names=F,sep="\t",col.names=T,quote=F)

# Similar analysis for hel2 and wt
p5_ds=p5_d[,c(1:3,4:6)]
mRNA_s=mRNA[,c(1,2,7,8)]
all=merge(p5_ds,mRNA_s,by=0,suffixes = c("_5p","_mRNA"))
all[is.na(all)] <- 0
rownames(all)=all$Row.names
all=all[,-1]
colData <- data.frame (row.names = colnames(all), condition = c("wt","wt","wt","mutant","mutant","mutant","mutant","mutant","wt","wt"),sampleType=c("5p","5p","5p","5p","5p","5p","mRNA","mRNA","mRNA","mRNA"))
dds <- DESeqDataSetFromMatrix (countData = all,colData = colData, design = ~ sampleType+condition+sampleType:condition)
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds$sampleType <- relevel (dds$sampleType, ref = "mRNA")
dds$condition <- relevel (dds$condition, ref = "wt")
dds <- DESeq(dds, test="LRT", reduced =  ~ sampleType+condition)
res <- results(dds,name="sampleType5p.conditionmutant",alpha=0.05)
res2 <- lfcShrink(dds, coef=4,type = "apeglm",res=res)
res2 <- res2[order(res2$padj), ]
res2=data.frame(gene_name=rownames(res2),res2)
write.table(res2,"/Users/jwu/Downloads/all_gene_5p_vs_mRNA_hel2_vs_wt_p0.05_shrink_rmdup.txt",row.names=F,sep="\t",col.names=T,quote=F)

# Similar analysis for ncs2elp6hel2 and ncs2elp6
p5_ds=p5_d[,c(7:9,10:12)]
mRNA_s=mRNA[,c(3,4,5,6)]
all=merge(p5_ds,mRNA_s,by=0,suffixes = c("_5p","_mRNA"))
all[is.na(all)] <- 0
rownames(all)=all$Row.names
all=all[,-1]
head(all)
colData <- data.frame (row.names = colnames(all), condition = c("wt","wt","wt","mutant","mutant","mutant","wt","wt","mutant","mutant"),sampleType=c("5p","5p","5p","5p","5p","5p","mRNA","mRNA","mRNA","mRNA"))
dds <- DESeqDataSetFromMatrix (countData = all,colData = colData, design = ~ sampleType+condition+sampleType:condition)
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds$sampleType <- relevel (dds$sampleType, ref = "mRNA")
dds$condition <- relevel (dds$condition, ref = "wt")
dds <- DESeq(dds, test="LRT", reduced =  ~ sampleType+condition)
res <- results(dds,name="sampleType5p.conditionmutant",alpha=0.05)
res2 <- lfcShrink(dds, coef=4,type = "apeglm",res=res)
res2 <- res2[order(res2$padj), ]
res2=data.frame(gene_name=rownames(res2),res2)
write.table(res2,"/Users/jwu/Downloads/all_gene_5p_vs_mRNA_ncs2elp6hel2_vs_ncs2elp6_p0.05_shrink_rmdup.txt",row.names=F,sep="\t",col.names=T,quote=F)

