#############Libraries#####
library(DESeq2)
library(pcaMethods)
library(ggplot2)
library(pheatmap)

#############Setup####
counts = read.delim("exp3_UTAP.txt", sep = "\t")[,c(1:9)]
rownames(counts) = counts$Gene
counts = counts[, -1]
coldata = data.frame(Treatment = c("PBS", "PBS","PBS", "PBS","PBS",  "LPS", "LPS", "LPS"))
rownames(coldata) = colnames(counts)

############DESEQ2 design######
outliers = NA
outlies_removed = "T"

deseq = DESeqDataSetFromMatrix(countData = counts[,setdiff(colnames(counts), outliers)], colData = data.frame(Treatment = coldata[setdiff(colnames(counts), outliers),]), design = ~ Treatment)

###########Gene filtering####
keep = rowSums(counts(deseq)) >= 30
deseq = deseq[keep,]
deseq = DESeq(deseq)

############QC######
ntd <- normTransform(deseq)
vsd <- vst(deseq, blind=FALSE)

select <- order(rowMeans(counts(deseq,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(deseq))
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


pcaData <- plotPCA(vsd, intgroup=c("Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Treatment)) +
  geom_point(size = 6) +
  geom_text(aes(label = substr(name,2,3)), color = "black") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#############DGE##########
logFC_cut = log2(1.25)
alpha = 0.05
factor = "Treatment"
g1 = "LPS"
g2 = "PBS"
DEG = results(deseq, contrast = c(factor, g1, g2), alpha = alpha, lfcThreshold = 0) 

summary(DEG)
mcols(DEG)$description

DESeq2::plotMA(DEG, ylim = c(-3, 3))

plotCounts(deseq, gene = which.min(DEG$padj), intgroup = factor)

d = plotCounts(deseq, gene = which.min(DEG$padj), intgroup=factor, returnData=TRUE)
ggplot(d, aes(x=Treatment, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle(rownames(DEG)[which.min(DEG$padj)])

DEG[order(DEG$pvalue),]

DEG_volcano = DEG[-which(is.na(DEG$padj)),]
ggde_table = data.frame(gene = rownames(DEG_volcano), x = DEG_volcano$log2FoldChange, y = DEG_volcano$padj)
ggde_table[is.na(ggde_table$y), "y"] = 1

selected_genes = c("Trem2", "Ch25h", "Orm1", "Cd177", "Clec7a", "Ecrg4", "Cd200r4", "Cd200r1", "Fcamr", "Cfb", "Wfdc17")

nudge = 0.4
force = 10
size = 3
ggplot(ggde_table, aes(x = x, y = -log10(y))) +
  geom_point(color = "grey", size = 2, shape = 19) +
  geom_point(data = ggde_table[ggde_table$gene %in% selected_genes,], size = 3, color = "red", shape = 19) +
  geom_text_repel(data = ggde_table[ggde_table$gene %in% selected_genes,], color = "black", aes(label = gene), segment.size  = 0.15, nudge_x = -nudge, segment.color = "grey", size = size, max.overlaps = 100, force = force) +
  ylim(c(0, quantile(-log10(ggde_table$y), 0.999))) +
  xlab(paste0("log2 ", g1, "/", g2)) +
  ylab("-log10 FDR") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") + 
  theme_classic() +
  theme(legend.position = "none",
        title = element_text(size = 30),
        axis.title = element_text(size = 20))
ggsave(filename = "volcano.pdf", width = 150, height = 150, units = "mm", dpi = 200)

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

fgsea = fgsea(gsea_pathways, gsea_list, minSize = 10  , maxSize = 500)

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

fgsea_selected = c("GO_CILIUM","GO_INNATE_IMMUNE_RESPONSE","GO_ADAPTIVE_IMMUNE_RESPONSE","GO_INFLAMMATORY_RESPONSE", "GO_REGULATION_OF_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","GO_REGULATION_OF_INTERLEUKIN_6_PRODUCTION","GO_REGULATION_OF_CYTOKINE_SECRETION","GO_REGULATION_OF_RESPONSE_TO_WOUNDING","GO_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB" ,"GO_MACROPHAGE_ACTIVATION", "GO_EXTRACELLULAR_MATRIX","GO_MORPHOGENESIS_OF_AN_EPITHELIUM")  #read.xlsx("Epith_GO_indicated_1.xlsx", sheet = 2,rowNames = F)
fgsea_selected = gse_browes[fgsea_selected,]
fgsea_selected$pathway = rownames(fgsea_selected)
gsea_order = fgsea_selected$pathway[order(fgsea_selected$NES)]
fgsea_selected = fgsea_selected[,c("pathway", "padj","NES","leading_edge")]

ggplot(fgsea_selected) +
  geom_bar(aes(y = factor(pathway, levels = gsea_order), x = NES, fill = -log10(padj)), stat = "identity") +
  scale_y_discrete(position = "right") + 
  theme_classic() +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red") +
  labs(x = "NES", title = "Infl pov vs. Ctrl") +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color = "black"))
ggsave(filename = "GSEA_selected.pdf", width = 160, height = 60, units = "mm", dpi = 200)

