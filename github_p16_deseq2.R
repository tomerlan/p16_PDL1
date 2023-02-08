#############Libraries#####
library(DESeq2)
library(pcaMethods)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)


#############Setup####
counts = read.delim("P16_exp2_UTAP.txt", sep = "\t")
rownames(counts) = counts$Gene
counts = counts[, -1]
coldata = data.frame(Treatment = substr(colnames(counts), 12, 15), p16 = substr(colnames(counts), 8, 10), Combined = substr(colnames(counts), 4, 15))
rownames(coldata) = colnames(counts)


############DESEQ2 design######
outliers =  NA
design = "controlled" # with control variable

deseq = DESeqDataSetFromMatrix(countData = counts[,setdiff(colnames(counts), outliers)], colData = coldata[setdiff(colnames(counts), outliers),], design = ~ Treatment + p16); 

deseq$Treatment = relevel(deseq$Treatment, ref = "LPS")
deseq$p16 = relevel(deseq$p16, ref = "neg")


###########Gene filtering####
keep = rowSums(counts(deseq)) >= 30
deseq = deseq[keep,]
deseq = DESeq(deseq)


############QC######
ntd <- normTransform(deseq)
vsd <- vst(deseq, blind=FALSE)

select <- order(rowMeans(counts(deseq,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(deseq)[,c("Treatment","p16")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$p16, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


pcaData <- plotPCA(vsd, intgroup=c("Treatment", "p16"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=p16)) +
  geom_point(size = 6) +
  geom_text(aes(label = substr(name,2,2)), color = "black") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

outliers = c("s1_p16_neg_PBS", "s1_p16_pos_PBS")

############DESEQ2 design remove outliers######
outlies_removed = "T"
design = "controlled" # with control variable

deseq = DESeqDataSetFromMatrix(countData = counts[,setdiff(colnames(counts), outliers)], colData = coldata[setdiff(colnames(counts), outliers),], design = ~ Treatment + p16);

deseq$Treatment = relevel(deseq$Treatment, ref = "LPS")
deseq$p16 = relevel(deseq$p16, ref = "neg")

keep = rowSums(counts(deseq)) >= 30
deseq = deseq[keep,]
deseq = DESeq(deseq)


#############DGE##########
logFC_cut = log2(1.25)
alpha = 0.05
factor = "p16"
g1 = "pos"
g2 = "neg"
DEG = results(deseq, contrast = c(factor, g1, g2), alpha = alpha, lfcThreshold = logFC_cut) 

summary(DEG)
mcols(DEG)$description

DESeq2::plotMA(DEG, ylim = c(-3, 3))

plotCounts(deseq, gene = which.min(DEG$padj), intgroup = factor)

d = plotCounts(deseq, gene = which.min(DEG$padj), intgroup=factor, returnData=TRUE)
ggplot(d, aes(x=p16, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle(rownames(DEG)[which.min(DEG$padj)])

DEG[order(DEG$pvalue),]

##volcano plot
EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = alpha,
                FCcutoff = 1.25,
                pointSize = 2,
                labSize = 4,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                max.overlaps = 100,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)
ggsave(paste0("volcano_",g1, "_",g2, "_",design,".png"), width = 200, height = 200, units = "mm", dpi = 200)



###########GSEA fgsea#####
logFC = DEG$log2FoldChange; names(logFC) = rownames(DEG)
gsea_list = sort(logFC, decreasing = T); head(gsea_list, 30)
egSYMBOL = AnnotationDbi::toTable(org.Mm.egSYMBOL)
egSYMBOL = egSYMBOL[-which(duplicated(egSYMBOL$symbol)),]
rownames(egSYMBOL) = egSYMBOL$symbol
names(gsea_list) = egSYMBOL[names(gsea_list),1]
gsea_list = gsea_list[-which(is.na(names(gsea_list)))]

gsea_pathways = get(load("mouse_c5_v5p2.RData"))
#"http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2022.1.Mm/m5.all.v2022.1.Mm.symbols.gmt"

fgsea = fgsea(gsea_pathways, gsea_list, minSize = 5  , maxSize = 500)

gene_key = egSYMBOL
rownames(gene_key) = gene_key$gene_id

gse_browes = as.matrix(fgsea[,-8])
rownames(gse_browes) = gse_browes[,1]
gse_browes = gse_browes[,-1]
class(gse_browes) <- "numeric"
gse_browes = data.frame(gse_browes)
gse_browes$leading_edge = rep(NA, nrow(gse_browes))

for(i in 1:nrow(fgsea)){
  gse_browes[i, "leading_edge"] <- paste(gene_key[fgsea$leadingEdge[[i]], "symbol"], collapse = " ")
}

fgsea_selected = c("GO_EXTRACELLULAR_MATRIX","GO_RESPONSE_TO_RETINOIC_ACID","GO_DNA_REPLICATION" , "GO_MYELOID_DENDRITIC_CELL_ACTIVATION","GO_REGULATION_OF_CHEMOKINE_SECRETION","GO_NEGATIVE_REGULATION_OF_RESPONSE_TO_WOUNDING")
fgsea_selected = gse_browes[fgsea_selected,]
fgsea_selected$pathway = rownames(fgsea_selected)
gsea_order = rev(fgsea_selected$pathway[order(fgsea_selected$NES)])

fgsea_selected = fgsea_selected[,c("pathway", "padj","NES","leading_edge")]
colnames(fgsea_selected) = c("term", "padj","NES","leading_edge")
ggplot(fgsea_selected) +
  geom_bar(aes(y = factor(term, levels = gsea_order), x = NES, fill = -log10(padj)), stat = "identity") +
  scale_y_discrete(position = "right") + 
  theme_classic() +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red") +
  labs(x = "NES", title = "P16 pov vs. neg") +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color = "black"))
ggsave(filename = "GSEA_selected_T2.pdf", width = 160, height = 50, units = "mm", dpi = 200)
