#Transcriptomic clustering of COVID-19 patients

#### 0. Setings ####
library(RColorBrewer)
library(ggplot2)
library(critcolors) #To install: > devtools::install_github("cecilomar6/critcolors")

cluster1<-"#486D87"
cluster2<-"#FFB838"
significant<-"blue"
non_significant<-"grey"
up_regulated<-"pink"
down_regulated<-"green"
clusters_colors<-c(cluster1, cluster2)
sig_colors<-critcolors(100)
exp_colors<-colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)


#################################################################################
############################## Clustering #######################################

#### 1. Loading sample info and transcript data ####

load("transcript_data.RData")
samples<-read.csv("clinical_data.csv")
samples<-samples[match(colnames(txi$abundance), samples$id ),] 

library(AnnotationHub)
ah<-AnnotationHub()
edb<-ah[["AH73986"]]


#### 2. Normalizing count matrix to TPM (or FPKM) ####

mat <- txi$abundance

#### 3. Transforming data to log2 ####

mat <- log2(mat +1)

#### 4. Selecting highly variable features ####

library(genefilter)
mat.filter <- varFilter(mat, var.cutoff = 0.95, filterByQuantile = TRUE)
dim(mat)
dim(mat.filter)

genenames<-AnnotationDbi::select(edb, keys=rownames(mat.filter), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(mat.filter)<-genenames$SYMBOL  

#### 5. Distance matrices ####

distance <- dist(t(mat.filter), method = "euclidean")

#### 6. Clustering ####

clusters <- hclust(distance, method = "ward.D") 

library(factoextra)

hc.cut<-hcut(distance, k=2, hc.method="ward.D", hc_func = "hclust")

fviz_dend(hc.cut,
          k=2,
          show_labels = FALSE,
          rect=TRUE,
          rect_fill=TRUE,
          k_colors=rev(clusters_colors),
          main=NULL,
          cex=32,
          ggtheme=theme_classic(base_size = 24))+
  theme(aspect.ratio=1.618)

library(pvclust)
cluster_ps<-pvclust((mat.filter),method.hclust = "ward.D", method.dist="euclidean", nboot=1000) #Cluster probabilities
plot(cluster_ps, labels=FALSE, print.num=FALSE)

clustercut <- cutree(clusters, 2)
samples$clust2 <- as.factor(clustercut)
table(samples$clust2)

library(gplots)

## Heatmap of genes used for clustering
sc.mat.filter <- t(apply(X = mat.filter,MARGIN = 1,FUN = function(y){
  (y - mean(y))/sd(y)}))

samples$colors<-clusters_colors[as.numeric(samples$clust2)]

heatmap.2(sc.mat.filter,
          hclustfun = function(x){hclust(x, method="ward.D")},
          scale = "none",
          col= exp_colors,
          keysize = 1.25,
          #lhei=c(1,5),
          trace = "none",
          density.info = "none",
          ColSideColors = samples$colors,
          labCol = NA,
          labRow = NA,
          margins=c(1,1))

## UMAP plot of clusters
library(umap)

my.umap <- umap::umap(t(mat.filter) )
dat <- as.data.frame(my.umap$layout)
dat$clust <- samples$clust2[match(rownames(dat), samples$id)]

ggplot(dat, aes(x = V1, y = V2,col = clust))+
  geom_point(size = 2.5)+
  scale_color_manual(values = clusters_colors)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio = 1, legend.position = "none")+
  labs(x="UMAP1", y="UMAP2")


#################################################################################
######################## Differential expression ################################

library(DESeq2)
ddsTxi <- DESeqDataSetFromTximport(txi, colData=samples, design= ~clust2)

## Prefiltering
keep<-rowSums(counts(ddsTxi))>=ncol(ddsTxi)
ddsTxi<-ddsTxi[keep,]

dds<-DESeq(ddsTxi)
genenames<-AnnotationDbi::select(edb, keys=rownames(dds), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(dds)<-genenames$SYMBOL  

resultsNames(dds)
res <- results(dds, name = "clust2_2_vs_1")

resOrdered <- res[order(res$log2FoldChange),]
resOrdered
#write.csv(resOrdered, file="suppl_file_1.csv")
summary(res, alpha=0.01)

sig.genes<-resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.01 & abs(resOrdered$log2FoldChange) > 2,]
sig.genes

#Volcano plot (after GSEA)



####################################################################################
######################## Gene Set Enrichment Analysis ##############################

# From 
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library("org.Hs.eg.db", character.only = TRUE)

# formatting data from deseq2
# Rereading dds to work with ENSEMBL IDs, as to avoid duplicates

dds.gsea <- DESeq(ddsTxi)
res.gsea <- results(dds.gsea, contrast = c("clust2", "1", "2"))
res.gsea <- as.data.frame(res.gsea)
sig.genes.gsea<-res.gsea[!is.na(res.gsea$padj) & res.gsea$padj<0.01 & abs(res.gsea$log2FoldChange)>2,]
res.gsea<-res.gsea[rownames(sig.genes.gsea),] # beware! use ENSG, not genename

original_gene_list <- -res.gsea$log2FoldChange #Reverse to represent CTP2 vs CTP1
names(original_gene_list) <- rownames(res.gsea)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

#GSE
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             minGSSize = 2, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")

## OUTPUTS

#dotplot
require(DOSE)

dotplot(gse, showCategory=gse@result$Description[1:55], split=".sign", font.size = 8) + 
  facet_grid(.~.sign)+
  scale_color_gradientn(name = "p.dajust",
                        colours = sig_colors,
                        trans="log",
                        limits= c(1e-8, 0.05), 
                        breaks=c(1e-8,1e-5, 1e-3, 0.05) ) 

dotplot(gse, showCategory=gse@result$Description[56:110], split=".sign", font.size = 8) + 
  facet_grid(.~.sign)+
  scale_color_gradientn(name = "p.adjust",
                        colours = sig_colors,
                        trans="log",
                        limits= c(1e-8, 0.05), 
                        breaks=c(1e-8,1e-5, 1e-3, 0.05) ) 

GOofInterest<-c("GO:0045589", "GO:0045066", "GO:0043331", "GO:0034340", "GO:0002286", "GO:0060337", "GO:0050871", "GO:0002323", "GO:0050853", "GO:0042100")

edox <- setReadable(gse, "org.Hs.eg.db", "ENSEMBL")


categ<-edox@result$Description[gse@result$ID %in% GOofInterest]

dotplot(edox, showCategory=categ, split=".sign", font.size = 24) + 
  facet_grid(.~.sign)+
  scale_color_gradientn(name = "p.dajust",
                        colours = sig_colors,
                        trans="log",
                        limits= c(1e-4, 0.05), 
                        breaks=c(0.0001,0.001, 0.01, 0.05) ) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

cnetplot(edox,showCategory=categ, foldChange = gene_list, 
             circular=FALSE, colorEdge=FALSE,
             layout="fr",
         color_category="#ffffff",
             cex_category=1,
             cex_label_category=0.5,
             cex_gene=0.5,
             cex_label_gene=0.5)+
  scale_color_gradientn(name = "fold change",
                           colours = exp_colors,
                           limits= c(-10, 5), 
                           breaks=c(-10 , 0, 5) ) +
  theme(aspect.ratio = 0.5)

#gene extraction
GenesofInterest<-unique(unlist(strsplit(edox@result$core_enrichment[edox@result$ID %in% GOofInterest], "/")))

#Volcano plot
library(ggrepel)
genes<-as.data.frame(res)
genes$genename<-rownames(genes)
genes$Significant <- ifelse((genes$padj <= 0.01 & abs(genes$log2FoldChange)>2), "FDR < 0.01", "Not Sig")
genes$Significant[genes$Significant=="FDR < 0.01" & genes$genename %in% GenesofInterest]<-"FDR < 0.01, IFN pathway"
ggplot(genes, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant), shape=1, alpha=(0.99-0.2*(genes$padj>0.01))) +
  geom_hline(yintercept = -log10(0.01), linetype=2)+
  scale_color_manual(values=c(sig_colors[100],"red", "grey75"))+
  xlim(-10,5)+
  theme_gray(base_size = 24) + 
  theme(legend.position = "none") +
  geom_text_repel(
    data = genes[genes$genename %in% GenesofInterest & genes$padj<0.05,],
    aes(label = genename),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps=30
  )

#############################################################################################
################ Corrplot from all GO terms related to type I IFN ###########################

# GOs extracted by searching "Type I Interferon" in AmiGO2 
IFN_GOs <- c("GO:0019962","GO:0032606","GO:0005132","GO:0004905","GO:0034340","GO:0038197",
             "GO:0060337","GO:0071357","GO:0032481","GO:0032480","GO:0032479","GO:0039501",
             "GO:0060340","GO:0060339","GO:0060338","GO:0035456","GO:0035455","GO:0039502",
             "GO:0035458","GO:0035457","GO:1990231")


library(biomaRt)
ensembl.hs <- useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = IFN_GOs, mart = ensembl.hs)

length(unique(gene.data$hgnc_symbol))

ifn.sig<-c(unique(gene.data$hgnc_symbol))

ifn.sig.expr<-dds[rownames(dds) %in% ifn.sig,]

#Comparing correlation matrices
library(corrplot)
library(Hmisc)
library(psych)
corr.matrix.CTP1<-cor(t(counts(ifn.sig.expr, normalized=TRUE)[,ifn.sig.expr$clust2==1]))
corr.matrix.CTP2<-cor(t(counts(ifn.sig.expr, normalized=TRUE)[,ifn.sig.expr$clust2==2]))

corr.matrix.test <- cortest(corr.matrix.CTP1, corr.matrix.CTP2, n1 = 42, n2 = 14)
corr.matrix.test

#Plotting only significant correlations
corr.matrix.CTP1<-rcorr(t(counts(ifn.sig.expr, normalized=TRUE)[,ifn.sig.expr$clust2==1]))
corr.matrix.CTP2<-rcorr(t(counts(ifn.sig.expr, normalized=TRUE)[,ifn.sig.expr$clust2==2]))


corrplot(corr.matrix.CTP1$r, p.mat = corr.matrix.CTP1$P, insig = "blank", order="hclust", type="lower", tl.cex = 0.5, tl.col="black")
corrplot1<-corrplot(corr.matrix.CTP1$r, order="hclust", type="lower", tl.cex = 0.5, tl.col="black")
corr.matrix.CTP2$r <- corr.matrix.CTP2$r[unique(corrplot1$corrPos$xName), unique(corrplot1$corrPos$xName)]
corr.matrix.CTP2$P <- corr.matrix.CTP2$P[unique(corrplot1$corrPos$xName), unique(corrplot1$corrPos$xName)]
corrplot(corr.matrix.CTP2$r, p.mat = corr.matrix.CTP2$P,insig = "blank", order="original", type="lower", tl.cex = 0.5, tl.col="black")

##### Pairwise comparisons of correlation coefficients
##### DGCA: package for differential correlation across conditions

library(DGCA)

inputMat <- counts(dds, normalized=TRUE)
inputMat <- as.data.frame(inputMat[rownames(inputMat) %in% ifn.sig,])
design_mat <- model.matrix(~ 0 + dds$clust2)
colnames(design_mat) <- c("CTP1", "CTP2")

cor_res <-  getCors(inputMat = inputMat, design = design_mat) 

dcPairs_res <-  pairwiseDCor(cor_res, compare = c("CTP1", "CTP2"))

ddcor_res <-  ddcorAll(inputMat = inputMat, design = design_mat,
                       compare = c("CTP1", "CTP2"),
                       adjust = "BH", heatmapPlot = TRUE, nPerm = 0, nPairs = "all")
head(ddcor_res)

#Filtering opposite correlations in each cluster
library(dplyr)
sig.ddcor_res <- ddcor_res %>% filter(pValDiff_adj < 0.05) %>% filter(Classes == "-/+" | Classes == "+/-")


### Making a net of our favorite correlations

library(ggnet) #Install with "devtools::install_github("briatte/ggnet")"
library(network)
library(sna)

### Create networks from https://www.jessesadler.com/post/network-analysis-with-r/

nodes <- as.data.frame(unique(c(sig.ddcor_res$Gene1, sig.ddcor_res$Gene2)))
colnames(nodes) <- "label"

edges <- as.data.frame(sig.ddcor_res[,c("Gene1", "Gene2", "CTP1_cor", "CTP2_cor")])

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

color_legend<-map2color(c(edges$CTP1_cor, edges$CTP2_cor), COL2())

edges$CTP1_col <- color_legend[1:round(length(color_legend)/2)]
edges$CTP2_col <- color_legend[round(length(color_legend)/2)+1:round(length(color_legend)/2)]

routes_network <- network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE)
layout.mode<-"fruchtermanreingold"
n1<-ggnet2(routes_network, mode=layout.mode, 
           label = TRUE, label.size = 5, edge.color = "CTP1_col", edge.size = 1)
n2<-ggnet2(routes_network, mode=matrix(c(n1$data$x, n1$data$y), ncol=2), 
           label = TRUE, label.size = 5, edge.color = "CTP2_col", edge.size = 1)
n1
n2

#Corrplots of genes of interest
corr.matrix.CTP1.small<-rcorr(t(counts(ifn.sig.expr, normalized=TRUE)[rownames(counts(ifn.sig.expr)) %in% nodes$label,ifn.sig.expr$clust2==1]))
corr.matrix.CTP2.small<-rcorr(t(counts(ifn.sig.expr, normalized=TRUE)[rownames(counts(ifn.sig.expr)) %in% nodes$label,ifn.sig.expr$clust2==2]))


corrplot(corr.matrix.CTP1.small$r, p.mat = corr.matrix.CTP1.small$P, insig = "blank", order="hclust", type="lower", tl.cex = 1, tl.col="black")
corrplot1.small<-corrplot(corr.matrix.CTP1.small$r, order="hclust", type="lower", tl.cex = 0.5, tl.col="black")
corr.matrix.CTP2.small$r <- corr.matrix.CTP2.small$r[unique(corrplot1.small$corrPos$xName), unique(corrplot1.small$corrPos$xName)]
corr.matrix.CTP2.small$P <- corr.matrix.CTP2.small$P[unique(corrplot1.small$corrPos$xName), unique(corrplot1.small$corrPos$xName)]
corrplot(corr.matrix.CTP2.small$r, p.mat = corr.matrix.CTP2.small$P,insig = "blank", order="original", type="lower", tl.cex = 1, tl.col="black")



######################################################################################
################### Cell Deconvolution using Immunostates ############################


library(cowplot)
library(MetaIntegrator)
library(data.table)
library(dplyr)

#Deconvolution of RNAseq samples using Immunostates
immunoStatesMatrix2<-immunoStatesMatrix[,!colnames(immunoStatesMatrix) %in% c("MAST_cell", "macrophage_m0", "macrophage_m1", "macrophage_m2")]
testRNA<-counts(dds, normalized=TRUE)

sum(rownames(testRNA) %in% rownames(immunoStatesMatrix)) #Out of 318 genes in immunoStatesMatrix
testRNA<-testRNA[rownames(testRNA) %in% rownames(immunoStatesMatrix),]

#collapse unique genes
testRNA<- rowsum(testRNA, rownames(testRNA))

#Immunostates on the expression matrix
outDT <- as.data.table(MetaIntegrator:::iSdeconvolution(immunoStatesMatrix2, 
                                                        testRNA), keep.rownames = T)

outDT[,natural_killer_cell:=CD56bright_natural_killer_cell+CD56dim_natural_killer_cell]
outDT[,           monocyte:=CD14_positive_monocyte+CD16_positive_monocyte]
outDT[,             B_cell:=naive_B_cell+memory_B_cell]
outDT[,             T_cell:=CD8_positive_alpha_beta_T_cell+CD4_positive_alpha_beta_T_cell+gamma_delta_T_cell]
outDT[,        granulocyte:=neutrophil+eosinophil+basophil]
outDT[, lymp:=B_cell+T_cell]

#Cleaning and re-scaling
cytokines<-colnames(outDT)[2:17]
keep<-apply(outDT[,2:17],2, function(x) sum(x>0)>10) #Only cell populations present in more than 10 samples
cytokines<-cytokines[keep]
cytokines<-cytokines[c(1,9,5,6,
                       2,4,7,8,
                       10,3,11)] #Reordering

outDT_scaled<-as.data.frame(outDT)
outDT_scaled<-outDT_scaled[,cytokines]
outDT_scaled<-as.data.frame(t(apply(outDT_scaled,1, FUN=function(x) x/sum(x))))
rowSums(outDT_scaled)


#Stats and Plotting

pcalc<-function(x) {
  wt<-wilcox.test(as.numeric(x)~as.factor(samples$clust2))
  return(round(wt$p.value, 3))
}

p_val<-apply(outDT_scaled, 2, pcalc)
p_val<-paste("p=", p_val, sep="")
p_val[p_val=="p=0"]<-"p<0.001"
p_val<-as.factor(p_val)

cyto_plot<-function(celltype, titulo, subtitulo) {
  cytoquine<-as.numeric(outDT_scaled[[celltype]])
  genot<-as.factor(samples$clust2)
  ggplot(NULL, aes(y=cytoquine, x=genot, col=genot, fill=genot))+
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
    geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
    theme_grey(base_size = 12)+
    theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
          strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
    labs(x=NULL, y=NULL, title=titulo)+
    scale_x_discrete(labels=c("CTP1", "CTP2"))+
    scale_y_continuous(labels=scales::percent)+
    scale_color_manual(values=clusters_colors)+
    scale_fill_manual(values=clusters_colors)+
    facet_wrap(subtitulo)
}

titulos<-cytokines

lista<-list()
for(i in 1:length(cytokines)) {
  lista[[i]]<-cyto_plot(cytokines[[i]], titulos[[i]], p_val[[i]])
  #ggsave(filename=paste(titulos[[i]],".pdf", sep=""), plot=lista[[i]])
}

plot_grid(plotlist=lista, align="hv", nrow=3)

#Major cell groups
ggplot(NULL, aes(x=samples$clust2, y=outDT$granulocyte, col=samples$clust2, fill=samples$clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="Granulocytes")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
wilcox.test((outDT$granulocyte)~samples$clust2)

ggplot(NULL, aes(x=samples$clust2, y=outDT$monocyte, col=samples$clust2, fill=samples$clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="Monocytes")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
wilcox.test((outDT$monocyte)~samples$clust2)

ggplot(NULL, aes(x=samples$clust2, y=outDT$lymp, col=samples$clust2, fill=samples$clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="Lymphocytes")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
wilcox.test((outDT$lymp)~samples$clust2)

ggplot(NULL, aes(x=samples$clust2, y=outDT$natural_killer_cell, col=samples$clust2, fill=samples$clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="Natural killer")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
wilcox.test((outDT$natural_killer_cell)~samples$clust2)

#Lymphocyte subpopulations
ggplot(NULL, aes(x=samples$clust2, y=outDT$CD4_positive_alpha_beta_T_cell/outDT$lymp, col=samples$clust2, fill=samples$clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="CD4+")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
wilcox.test((outDT$CD4_positive_alpha_beta_T_cell/outDT$lymp)~samples$clust2)

ggplot(NULL, aes(x=samples$clust2, y=outDT$CD8_positive_alpha_beta_T_cell/outDT$lymp, col=samples$clust2, fill=samples$clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="CD8+")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
wilcox.test((outDT$CD8_positive_alpha_beta_T_cell/outDT$lymp)~samples$clust2)

ggplot(NULL, aes(x=samples$clust2, y=outDT$naive_B_cell/outDT$lymp, col=samples$clust2, fill=samples$clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="Naive B-cell")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
wilcox.test((outDT$naive_B_cell/outDT$lymp)~samples$clust2)

ggplot(NULL, aes(x=samples$clust2, y=outDT$memory_B_cell/outDT$lymp, col=samples$clust2, fill=samples$clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="Memory B-cell")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
wilcox.test((outDT$memory_B_cell/outDT$lymp)~samples$clust2)


#########################################################################################
###################### miRNAs targets for metascore genes ###############################

# Counts matrix

df_raw<-read.csv("miRNA_counts.csv")
rownames(df_raw)<-df_raw$X
df_raw<-df_raw[,-1]
colnames(df_raw) <- gsub(x = colnames(df_raw), pattern = "X", replacement = "")


# Samples info

mir.samples <- samples[, c("id", "clust2")]

mir.samples <- mir.samples[mir.samples$id %in% colnames(df_raw),]
df_raw <- df_raw[,colnames(df_raw) %in% mir.samples$id]
mir.samples <- mir.samples[match(colnames(df_raw),mir.samples$id),]

all(mir.samples$patient_id %in% colnames(df_raw))
all(mir.samples$patient_id == colnames(df_raw))

## Normalization using DESeq2 
dds.mir <- DESeqDataSetFromMatrix(countData = df_raw,
                                  colData = mir.samples,
                                  design = ~ clust2)

keep<-rowSums(counts(dds.mir))>=ncol(dds.mir)
dds.mir<-dds.mir[keep,]
#dds <- dds[top20_mi_rna,]

dds.mir<-DESeq(dds.mir)

## Read IPA output
library(stringr)
int.miRNAs.up <- readLines("molecules_upnetwork.txt")
int.miRNAs.up <- int.miRNAs.up %>% str_subset("miR") %>% str_replace(" \\(.*$", "")
length(int.miRNAs.up)

full.names.up <- str_subset(rownames(dds.mir),paste(int.miRNAs.up, collapse = "|"))
full.names.up<-full.names.up[!duplicated(str_extract(full.names.up, ".*? "))]
length(full.names.up)

l <- list()
for ( i in full.names.up){
  l[[i]] <- plotCounts(dds.mir, i, intgroup = "clust2", returnData = TRUE)
}

library(reshape2)
l <- melt(l)
l$L1<-gsub(" .*.?", "", l$L1)


library(critcolors)

pcalc_mirna<-function(x) {
  wt<-wilcox.test(x$value~x$clust2)
  return(round(wt$p.value, 3))
}
pvalues<-as.numeric(by(l, l$L1, FUN=pcalc_mirna))
pvalues

ggplot(l, aes(x = clust2, y = value, col = clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 12)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL)+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_y_continuous(labels=scales::percent)+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)+
  scale_y_log10()+
  facet_wrap(~L1, scales = "free")


## int.genes down
int.miRNAs.down <- readLines("molecules_downnetwork.txt")
int.miRNAs.down <- int.miRNAs.down %>% str_subset("miR") %>% str_replace(" \\(.*$", "")
length(int.miRNAs.down)

full.names.down <- str_subset(rownames(dds.mir),paste(int.miRNAs.down, collapse = "|"))
full.names.down<-full.names.down[!duplicated(str_extract(full.names.down, ".*? "))]
length(full.names.down)

l <- list()
for ( i in full.names.down){
  l[[i]] <- plotCounts(dds.mir, i, intgroup = "clust2", returnData = TRUE)
}

library(reshape2)
l <- melt(l)
l$L1<-gsub(" .*.?", "", l$L1)

calc_mirna<-function(x) {
  wt<-wilcox.test(x$value~x$clust2)
  return(round(wt$p.value, 3))
}
pvalues<-as.numeric(by(l, l$L1, FUN=pcalc_mirna))
pvalues


ggplot(l, aes(x = clust2, y = value, col = clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 12)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL)+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_y_continuous(labels=scales::percent)+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)+
  scale_y_log10()+
  scale_y_log10()+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  facet_wrap(~L1, scales = "free")


##########################################################################################
####################### Clinical differences and outcome #################################

samples$n_l_ratio<-samples$neutrophils/samples$lymphs

median_iqr<-function(x) {
  percentiles<-round(quantile(x, probs=c(0.25, 0.75), na.rm=TRUE),2)
  paste(round(median(x, na.rm=TRUE),2), " (", percentiles[1]," - ", percentiles[2],")", sep="")
}

quant<-c("age", "bmi","apache2", "FiO2_1",  "PaFi_1", "DAa_1", 
         "pCO2_1",  "RR_1", 
         "pH_1",  "lact_1", 
         "Vt_1", "vt_ibw_1",  "vei_1", 
         "Pplat_1",   "PEEP_1", 
         "driving_1", 
         "Crs_1", "cr_1", "PCT_ing", 
         "il6_1", "il6_3", "ferritin_1", "ferritin_3",
         "leuc_1", "lymph_1",  "Ddimer_1", "Ddimer_3", "preICU", "vfd28",
         "CK", "LDH", "AST", "ALT", "CRP", "neutrophils", "lymphs", "monos", "n_l_ratio")

quant_results<-data.frame(variable=character(), C1=character(), C2=character(), C3=character(), p_value=numeric())

for(i in quant) {
  values<-aggregate(as.numeric(samples[[i]])~samples$clust2, FUN=median_iqr)[,2]
  p_val<-round(wilcox.test(as.numeric(samples[[i]])~samples$clust2)$p.value,3)
  df.line<-data.frame(variable=i, CTP1=values[1], CTP2=values[2], p_value=p_val)
  quant_results<-rbind(quant_results, df.line)
  
}
View(quant_results)

qual<-c("sex", "hypert", "diabetes", "dislip", "tto_hypert", "tto_diabet", "statin",
        "frailty", "race", "CKD", "COPD", "cirrhosis", "immunosup", "neoplasia",
        "ventmode_1", 
        "vasoactive_1", 
        "early_prone", "ever_prone", "intubated", "nmba", "ecmo", "no_steroids")

qual_results<-data.frame(variable=character(), 
                         values=character(), 
                         C1=numeric(), 
                         C2=numeric(),
                         p.value=character())

for(i in qual) {
  values<-table(samples[[i]], samples$clust2)
  p_val<-round(chisq.test(values)$p.value,3)
  df.chunk<-data.frame(variable=c(i, rep("", nrow(values)-1)),
                       values=rownames(values),
                       C1=as.numeric(values[,"1"]), 
                       C2=as.numeric(values[,"2"]),
                       p.value=c(p_val, rep("", nrow(values)-1)))
  qual_results<-rbind(qual_results, df.chunk)
  
}

View(qual_results)

library(survival)
library(ggfortify)

clust_var <- as.factor(samples$clust2)
samples$icu_mort[is.na(samples$icu_mort)]<-"-"
table(samples$icu_mort, clust_var, useNA = "always")
samples$icu_mort<-factor(samples$icu_mort, levels=c("-", "0", "1"))
table(samples$icu_mort)

covid.surv.icu<-Surv(samples$follow_up.icu, samples$icu_mort)
autoplot(survfit(covid.surv.icu~clust_var)[,2], conf.int = FALSE, surv.size=2)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio = 1/1.618, legend.box.background = element_rect(), legend.position = c(0.85,0.25))+
  labs(x="Time (days)", y="ICU survival", col="Cluster")+
  scale_color_manual(values = clusters_colors, labels=c("CTP1", "CTP2") )

summary(survfit(covid.surv.icu~clust_var), times=c(0,25,50,75,100,125,150))
model<-coxph(covid.surv.icu~clust_var+age+sex+intubated, data=samples, id=id) 
summary(model)

######################################################################################
################################### Validation #######################################

#Definition of a specific gene signature for cluster allocation

genenames <- rownames(res)[rownames(res) %in% rownames(mat.filter) & res$padj < 0.01 & res$log2FoldChange>0]

dds.counts <- counts(dds, normalized=TRUE)
dds.counts <- t(dds.counts[rownames(dds.counts) %in% genenames,])

clust2 <- samples$clust2

library(pROC)

area.num <- apply(X = dds.counts, MARGIN = 2, FUN = function(x){
  stretch.glm <- glm(clust2 ~ x, family = binomial)
  probs <- predict(stretch.glm, type="response")
  roc.obj<-roc(clust2, probs, ci=FALSE)
  as.numeric(roc.obj$auc)
} )

area.num[area.num >0.95]
cluster.signature<-names(area.num[area.num>0.95])
train.up <- counts(dds, normalized = TRUE)
train.up <- train.up[rownames(train.up) %in% cluster.signature,]

gm_mean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[x>0 & !is.na(x)]))
  }
}

score <- apply(train.up, MARGIN = 2, FUN = gm_mean)

samples$score <- score
samples$st.score <- (samples$score-mean(samples$score))/sd(samples$score)

ggplot(samples, aes(x = score, fill = clust2))+
  geom_density(alpha =0.5)

ggplot(samples, aes(y = score, x =clust2, col=clust2, fill=clust2))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y="Transcriptomic score", title=NULL)+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
 
aggregate(score~clust2, samples, FUN=summary)
wilcox.test(score~clust2, data=samples)  #As expected

clust.glm<-glm(clust2~score, data=samples, family=binomial)
probs<-predict(clust.glm, type="response")
roc.obj<-roc(samples$clust2, probs, ci=TRUE)

ggroc<-function(roc.obj, col="black", size=1) {
  datapoints<-data.frame(specificity=roc.obj$specificities, sensitivity=roc.obj$sensitivities)
  ggplot(datapoints, aes(x=(1-specificity), y=sensitivity))+
    geom_path(col=col, size=size)+
    geom_segment(aes(x=0, xend=1, y=0, yend=1), color="grey", linetype=2)+
    theme_grey(base_size = 24)+
    theme(aspect.ratio = 1)+
    labs(x="1 - Specificity", y="Sensitivity")+
    annotate(geom="text", x=0.70, y=0.2, label=paste("AUC =",round(roc.obj$auc, digits=3)))+
    annotate(geom="text", x=0.70, y=0.1,
             label=paste("95% CI:",round(roc.obj$ci, digits=3)[1],"-", round(roc.obj$ci[3], digits=3)))
}

ggroc(roc.obj)

#Transcriptomic score as risk factor
model<-coxph(covid.surv.icu~I(score/100)+age+sex+intubated, data=samples, id=id, model=TRUE) 
summary(model)

#####Survival curves per score quartile #####

dummy_df<-data.frame(
  score=rep(quantile(samples$score)),
  age=mean(samples$age),
  sex=0.5,
  intubated=TRUE)

m<-(survfit(model, newdata=dummy_df))
cox_df<-data.frame(time=m$time, pstate=as.vector(m$pstate[,,2]), score=rep(1:5, each=37))

ggplot(cox_df, aes(x=time, y=pstate, group=(score), col=score)) +
  geom_step(size=2)+
  coord_cartesian(ylim=c(0,1))+
  theme_gray(base_size = 24)+
  theme(aspect.ratio = 1/1.618, legend.text=element_text(size=14), legend.box.background = element_rect(), legend.position = c(0.85,0.25))+
  labs(x="Time (days)", y="ICU survival", col="Score")+
  scale_color_gradient(low=clusters_colors[1], high = clusters_colors[2], labels=round(quantile(samples$score)) )


#Testing neutrophil/lympocyte ratio as marker
n_l_ratio_d<-outDT$neutrophil/outDT$lymp
clust.glm.nl<-glm(clust2~n_l_ratio_d, data=samples, family=binomial)
probs.nl<-predict(clust.glm.nl, type="response")
roc.obj.nl<-roc(samples$clust2, probs.nl, ci=TRUE)
roc.obj.nl

clust.glm.neu.m<-glm(clust2~neutrophils, data=samples, family=binomial)
probs.neu.m<-predict(clust.glm.neu.m, type="response")
roc.obj.neu.m<-roc(samples$clust2, probs.neu.m, ci=TRUE)
roc.obj.neu.m


model.neu<-coxph(covid.surv.icu~neutrophils+age+sex+intubated, data=samples, id=id, model=TRUE) 
summary(model.neu)

clust.glm.nl.m<-glm(clust2~n_l_ratio, data=samples, family=binomial)
probs.nl.m<-predict(clust.glm.nl.m, type="response")
roc.obj.nl.m<-roc(samples$clust2, probs.nl.m, ci=TRUE)
roc.obj.nl.m

model.nl<-coxph(covid.surv.icu~n_l_ratio+age+sex+intubated, data=samples, id=id, model=TRUE) 
summary(model.nl)

#Combination of transcriptomics and n_l ratio
model.3<-coxph(covid.surv.icu~clust_var+n_l_ratio+age+sex+intubated+il6_1, data=samples, id=id) 
summary(model.3)

model.4<-coxph(covid.surv.icu~score+n_l_ratio+age+sex+intubated+il6_1, data=samples, id=id) 
summary(model.4)

model.5<-coxph(covid.surv.icu~score+age+sex+intubated+il6_1, data=samples, id=id) 
summary(model.5)

#Using CRP as marker for clustering
clust.glm.crp<-glm(clust2~CRP, data=samples[!is.na(samples$CRP),], family=binomial)
probs.crp<-predict(clust.glm.crp, type="response")
roc.obj.crp<-roc(samples$clust2[!is.na(samples$CRP)], probs.crp, ci=TRUE)
roc.obj.crp

model.crp<-coxph(covid.surv.icu~CRP+age+sex+intubated, data=samples, id=id, model=TRUE) 
summary(model.crp)
                                    
                                    
##########################################################################################################
################################### Validation data from GSE157103 #######################################

library(GEOquery)
accession <- "GSE157103"
Sys.setenv("VROOM_CONNECTION_SIZE"=300000)
gset <- getGEO(accession, GSEMatrix = TRUE)
eset <- gset[[1]]

#pheno
pheno.1 <- pData(eset)

#counts
filePaths<-getGEOSuppFiles(accession)
count.matrix <- read.delim(rownames(filePaths)[1])
rownames(count.matrix) <- count.matrix$X.symbol
count.matrix <- count.matrix[,-1]
all(colnames(count.matrix) == pheno.1$description)
colnames(count.matrix) <- rownames(pheno.1)

pheno.1 <- pheno.1[pheno.1$characteristics_ch1 =="disease state: COVID-19",]
pheno.1 <- pheno.1[pheno.1$characteristics_ch1.3 =="icu: yes",]

count.matrix <- count.matrix[,colnames(count.matrix) %in% rownames(pheno.1)]
all(colnames(count.matrix) == rownames(pheno.1))

#Deseq to normalize counts
trans.matrix <- mutate_if(count.matrix, is.numeric, round)

dds.val <- DESeqDataSetFromMatrix(trans.matrix, pheno.1, ~1)

keep<-rowSums(counts(dds.val))>=ncol(dds.val)
dds.val <- dds.val[keep,]

dds.val <- DESeq(dds.val)
dds.counts <- counts(dds.val, normalized = TRUE)

#Calculate metascore in validation cohort and clustering
greedy.counts.val <- dds.counts[rownames(dds.counts) %in% cluster.signature,]

pheno.1$score <- apply(greedy.counts.val, MARGIN = 2, FUN = gm_mean)

pheno.1$clust <- factor(pheno.1$score>250, labels=c("CTP1", "CTP2"))
pheno.1$clust<-as.factor(pheno.1$clust)


ggplot(pheno.1, aes(y = score, x =clust, col=clust, fill=clust))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y="Transcriptomic score", title=NULL)+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)


#Outcomes
ggplot(pheno.1, aes(y = as.numeric(pheno.1$`apacheii:ch1`), x =clust, col=clust, fill=clust))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y="APACHE-II score", title=NULL)+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)

aggregate(as.numeric(pheno.1$`apacheii:ch1`)~clust, data=pheno.1, FUN=summary)
wilcox.test(as.numeric(pheno.1$`apacheii:ch1`)~clust, data=pheno.1)

ggplot(pheno.1, aes(y = as.numeric(pheno.1$`sofa:ch1`), x =clust, col=clust, fill=clust))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y="SOFA score", title=NULL)+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)

aggregate(as.numeric(pheno.1$`sofa:ch1`)~clust, data=pheno.1, FUN=summary)
wilcox.test(as.numeric(pheno.1$`sofa:ch1`)~clust, data=pheno.1)


ggplot(pheno.1, aes(y = as.numeric(pheno.1$`ventilator-free days:ch1`), x =clust, col=clust, fill=clust))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y="VFD at day 28", title=NULL)+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)

aggregate(as.numeric(pheno.1$`ventilator-free days:ch1`)~clust, data=pheno.1, FUN=summary)
wilcox.test(as.numeric(pheno.1$`ventilator-free days:ch1`)~clust, data=pheno.1)

quant.val<-c("age (years):ch1", "score", "apacheii:ch1", "sofa:ch1", "ventilator-free days:ch1")

quant_results.val<-data.frame(variable=character(), CTP1=character(), CTP2=character(), p_value=numeric())

for(i in quant.val) {
  values<-aggregate(as.numeric(pheno.1[[i]])~pheno.1$clust, FUN=median_iqr)[,2]
  p_val<-round(wilcox.test(as.numeric(pheno.1[[i]])~pheno.1$clust)$p.value,3)
  df.line<-data.frame(variable=i, CTP1=values[1], CTP2=values[2], p_value=p_val)
  quant_results.val<-rbind(quant_results.val, df.line)
  
}


t<-table(pheno.1$clust, pheno.1$`Sex:ch1`)
print(t)
fisher.test(t)

t<-table(pheno.1$clust, pheno.1$`ventilator-free days:ch1`==0)
print(t)
fisher.test(t)

#Cell populations in validation cohort
testRNA.val<-counts(dds.val, normalized=TRUE)
sum(rownames(testRNA.val) %in% rownames(immunoStatesMatrix2)) #Out of 318 genes in immunoStatesMatrix
testRNA.val<-testRNA.val[rownames(testRNA.val) %in% rownames(immunoStatesMatrix2),]

#collapse unique genes
testRNA.val<- rowsum(testRNA.val, rownames(testRNA.val))

#Immunostates on the expression matrix
outDT.val <- as.data.table(MetaIntegrator:::iSdeconvolution(immunoStatesMatrix2, 
                                                            testRNA.val), keep.rownames = T)

outDT.val[,natural_killer_cell:=CD56bright_natural_killer_cell+CD56dim_natural_killer_cell]
outDT.val[,           monocyte:=CD14_positive_monocyte+CD16_positive_monocyte]
outDT.val[,             B_cell:=naive_B_cell+memory_B_cell]
outDT.val[,             T_cell:=CD8_positive_alpha_beta_T_cell+CD4_positive_alpha_beta_T_cell+gamma_delta_T_cell]
outDT.val[,        granulocyte:=neutrophil+eosinophil+basophil]
outDT.val[, lymp:=B_cell+T_cell]

#Cleaning and re-scaling
cytokines.val<-colnames(outDT.val)[2:17]
keep<-apply(outDT.val[,2:17],2, function(x) sum(x>0)>10) #Only cell populations present in more than 10 samples
cytokines.val<-cytokines.val[keep]
cytokines.val<-cytokines.val[c(1,2, 10, 7,
                               3, 6, 9,
                               4,5,11,12)] #Reordering

outDT_scaled.val<-as.data.frame(outDT.val)
outDT_scaled.val<-outDT_scaled.val[,cytokines.val]
outDT_scaled.val<-as.data.frame(t(apply(outDT_scaled.val,1, FUN=function(x) x/sum(x))))
rowSums(outDT_scaled.val)


#Stats and Plotting

pcalc.val<-function(x) {
  wt<-wilcox.test(as.numeric(x)~as.factor(pheno.1$clust))
  return(round(wt$p.value, 3))
}

p_val<-apply(outDT_scaled.val, 2, pcalc.val)
p_val<-paste("p=", p_val, sep="")
p_val<-as.factor(p_val)

cyto_plot.val<-function(celltype, titulo, subtitulo) {
  cytoquine<-as.numeric(outDT_scaled.val[[celltype]])
  genot<-as.factor(pheno.1$clust)
  ggplot(NULL, aes(y=cytoquine, x=genot, col=genot, fill=genot))+
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
    geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
    theme_grey(base_size = 16)+
    theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
          strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
    labs(x=NULL, y=NULL, title=titulo)+
    scale_x_discrete(labels=c("CTP1", "CTP2"))+
    scale_y_continuous(labels=scales::percent)+
    scale_color_manual(values=clusters_colors)+
    scale_fill_manual(values=clusters_colors)+
    facet_wrap(subtitulo)
}

titulos<-cytokines.val

lista<-list()
for(i in 1:length(cytokines.val)) {
  lista[[i]]<-cyto_plot.val(cytokines.val[[i]], titulos[[i]], p_val[[i]])
 }
plot_grid(plotlist=c(lista[1:7], list(NULL), lista[8:11]), align="hv", nrow=3)


plot_grid(plotlist=lista[c(3,5,6,7)], align="hv", nrow=1)


ggplot(NULL, aes(x=pheno.1$clust, y=outDT.val$CD4_positive_alpha_beta_T_cell/outDT.val$lymp, col=pheno.1$clust, fill=pheno.1$clust))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="CD4+")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)

wilcox.test((outDT.val$CD4_positive_alpha_beta_T_cell/outDT.val$lymp)~pheno.1$clust)

ggplot(NULL, aes(x=pheno.1$clust, y=outDT.val$CD8_positive_alpha_beta_T_cell/outDT.val$lymp, col=pheno.1$clust, fill=pheno.1$clust))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="CD8+")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)

wilcox.test((outDT.val$CD8_positive_alpha_beta_T_cell/outDT.val$lymp)~pheno.1$clust)

ggplot(NULL, aes(x=pheno.1$clust, y=outDT.val$naive_B_cell/outDT.val$lymp, col=pheno.1$clust, fill=pheno.1$clust))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="Naive B-cell")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)

wilcox.test((outDT.val$naive_B_cell/outDT.val$lymp)~pheno.1$clust)

ggplot(NULL, aes(x=pheno.1$clust, y=outDT.val$memory_B_cell/outDT.val$lymp, col=pheno.1$clust, fill=pheno.1$clust))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y=NULL, title="Memory B-cell")+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)
wilcox.test((outDT.val$memory_B_cell/outDT.val$lymp)~pheno.1$clust)


# Validation from PMID: 35216673
# Data downloaded from https://zenodo.org/record/6120249

download.file("https://zenodo.org/record/6120249/files/CBD-KEY-CLINVAR.tar.gz", "clin_data.tar.gz")
untar("clin_data.tar.gz", "CBD-KEY-CLINVAR/COMBAT_CLINVAR_for_processed.txt")
pheno.2<-read.delim("CBD-KEY-CLINVAR/COMBAT_CLINVAR_for_processed.txt")

options(timeout = 360)
download.file("https://zenodo.org/record/6120249/files/CBD-KEY-RNASEQ-WB.tar.gz", "rnaseq.data.tar.gz")
untar("rnaseq.data.tar.gz", "CBD-KEY-RNASEQ-WB/Raw_count_data_143_60683.txt")
count.matrix<-read.delim("CBD-KEY-RNASEQ-WB/Raw_count_data_143_60683.txt")


pheno.2$RNASeq_sample_ID<-gsub("-", ".", pheno.2$RNASeq_sample_ID)
pheno.2<-pheno.2[pheno.2$RNASeq_sample_ID %in% colnames(count.matrix),]
pheno.2<-pheno.2[pheno.2$Source %in% c("COVID_CRIT","COVID_SEV"),] #Filtering severe cases

count.matrix<-count.matrix[, colnames(count.matrix) %in% pheno.2$RNASeq_sample_ID]

dds.val.2 <- DESeqDataSetFromMatrix(count.matrix, pheno.2, ~1)
keep<-rowSums(counts(dds.val.2))>=ncol(dds.val.2)
dds.val.2 <- dds.val.2[keep,]

dds.val.2 <- DESeq(dds.val.2)
dds.counts.2 <- counts(dds.val.2, normalized = TRUE)

symbols <- mapIds(org.Hs.eg.db, keys = rownames(dds.counts.2), keytype = "ENSEMBL", column="SYMBOL")
rownames(dds.counts.2)<-symbols

greedy.counts.val.2 <- dds.counts.2[rownames(dds.counts.2) %in% cluster.signature,]

pheno.2$score <- apply(greedy.counts.val.2, MARGIN = 2, FUN = gm_mean)

#Rescaling 250 (in a range from 36 to 830) to a range from 939 to 3254 results in a cutoff point of 1641

pheno.2$clust <- factor(pheno.2$score>1641, labels=c("CTP1", "CTP2"))
pheno.2$clust<-as.factor(pheno.2$clust)

ggplot(pheno.2, 
       aes(y = score, x =clust, fill=clust))+
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
  geom_jitter(aes(col=as.factor(Death28)), width=0.15, height=0, stroke=0, size=3.5)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
        strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
  labs(x=NULL, y="Transcriptomic score", title=NULL)+
  scale_x_discrete(labels=c("CTP1", "CTP2"))+
  scale_color_manual(values=clusters_colors)+
  scale_fill_manual(values=clusters_colors)

t<-table(pheno.2$clust, pheno.2$Death28)
t
chisq.test(t)

t<-table(as.factor(pheno.2$Age>=5), pheno.2$clust)
t
chisq.test(t)

t<-table(as.factor(pheno.2$Sex), pheno.2$clust)
t
chisq.test(t)

aggregate(score~clust, data=pheno.2, FUN=summary)
wilcox.test(score~clust, data=pheno.2)

#Cell populations in validation cohort 2
testRNA.val<-dds.counts.2
sum(rownames(testRNA.val) %in% rownames(immunoStatesMatrix2)) #Out of 318 genes in immunoStatesMatrix
testRNA.val<-testRNA.val[rownames(testRNA.val) %in% rownames(immunoStatesMatrix2),]

#collapse unique genes
testRNA.val<- rowsum(testRNA.val, rownames(testRNA.val))

#Immunostates on the expression matrix
outDT.val <- as.data.table(MetaIntegrator:::iSdeconvolution(immunoStatesMatrix2, 
                                                            testRNA.val), keep.rownames = T)

outDT.val[,natural_killer_cell:=CD56bright_natural_killer_cell+CD56dim_natural_killer_cell]
outDT.val[,           monocyte:=CD14_positive_monocyte+CD16_positive_monocyte]
outDT.val[,             B_cell:=naive_B_cell+memory_B_cell]
outDT.val[,             T_cell:=CD8_positive_alpha_beta_T_cell+CD4_positive_alpha_beta_T_cell+gamma_delta_T_cell]
outDT.val[,        granulocyte:=neutrophil+eosinophil+basophil]
outDT.val[, lymp:=B_cell+T_cell]

#Cleaning and re-scaling
cytokines.val<-colnames(outDT.val)[2:17]
keep<-apply(outDT.val[,2:17],2, function(x) sum(x>0)>10) #Only cell populations present in more than 10 samples
cytokines.val<-cytokines.val[keep]
cytokines.val<-cytokines.val[c(1,2, 7,
                               3, 5, 6,
                               8, 4,9)] #Reordering

outDT_scaled.val<-as.data.frame(outDT.val)
outDT_scaled.val<-outDT_scaled.val[,cytokines.val]
outDT_scaled.val<-as.data.frame(t(apply(outDT_scaled.val,1, FUN=function(x) x/sum(x))))
rowSums(outDT_scaled.val)


#Stats and Plotting

pcalc.val<-function(x) {
  wt<-wilcox.test(as.numeric(x)~as.factor(pheno.2$clust))
  return(round(wt$p.value, 3))
}

p_val<-apply(outDT_scaled.val, 2, pcalc.val)
p_val<-paste("p=", p_val, sep="")
p_val<-as.factor(p_val)

cyto_plot.val<-function(celltype, titulo, subtitulo) {
  cytoquine<-as.numeric(outDT_scaled.val[[celltype]])
  genot<-as.factor(pheno.2$clust)
  ggplot(NULL, aes(y=cytoquine, x=genot, col=genot, fill=genot))+
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
    geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
    theme_grey(base_size = 16)+
    theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
          strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
    labs(x=NULL, y=NULL, title=titulo)+
    scale_x_discrete(labels=c("CTP1", "CTP2"))+
    scale_y_continuous(labels=scales::percent)+
    scale_color_manual(values=clusters_colors)+
    scale_fill_manual(values=clusters_colors)+
    facet_wrap(subtitulo)
}

titulos<-cytokines.val

lista<-list()
for(i in 1:length(cytokines.val)) {
  lista[[i]]<-cyto_plot.val(cytokines.val[[i]], titulos[[i]], p_val[[i]])
}
plot_grid(plotlist=lista, align="hv", nrow=3)


##################################################################################
################## cluster-specific response to steroids #########################

load("transcript_data_trt.RData")
samples.trt <- read.csv("index_data_trt.csv")

################
### For CTP1 ###
################

library(tidyverse)
samples.trt.CTP1 <- samples.trt %>% filter(CTP == 1)
samples.trt.CTP1$corticosteroid_1 <- as.factor(samples.trt.CTP1$corticosteroid_1)
rownames(samples.trt.CTP1) <- samples.trt.CTP1$sample_id

txi.trt.CTP1 <- txi.trt
txi.trt.CTP1$abundance <- txi.trt.CTP1$abundance[,colnames(txi.trt.CTP1$abundance) %in% samples.trt.CTP1$sample_id]
txi.trt.CTP1$counts <- txi.trt.CTP1$counts[,colnames(txi.trt.CTP1$counts) %in% samples.trt.CTP1$sample_id]
txi.trt.CTP1$length <- txi.trt.CTP1$length[,colnames(txi.trt.CTP1$length) %in% samples.trt.CTP1$sample_id]

all(rownames(samples.trt.CTP1) == colnames(txi.trt.CTP1$abundance))

ddstxi.trt <- DESeqDataSetFromTximport(txi.trt.CTP1, colData=samples.trt.CTP1, design= ~corticosteroid_1)

## Prefiltering
dim(ddstxi.trt)
keep<-rowSums(counts(ddstxi.trt))>=ncol(ddstxi.trt)
ddstxi.trt<-ddstxi.trt[keep,]
dim(ddstxi.trt)

#### GSEA pathways for CTP1

# formatting data from deseq2
# Rereading dds to work with ENSEMBL IDs, to avoid duplicates

dds.gsea <- DESeq(ddstxi.trt)
resultsNames(dds.gsea)
res.gsea <- results(dds.gsea, name = "corticosteroid_1_1_vs_0")
res.gsea <- as.data.frame(res.gsea)

original_gene_list <- res.gsea$log2FoldChange
names(original_gene_list) <- rownames(res.gsea)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
length(gene_list)

ctp1.up<-rownames(res.gsea[res.gsea$pvalue<0.05 & res.gsea$log2FoldChange>0,])
ctp1.down<-rownames(res.gsea[res.gsea$pvalue<0.05 & res.gsea$log2FoldChange<0,])

#GSE
gse.CTP1 <- gseGO(geneList=gene_list, 
                  ont ="BP", 
                  keyType = "ENSEMBL", 
                  minGSSize = 10, 
                  maxGSSize = 800, 
                  pvalueCutoff = 1, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")

################################################################################

################
### For CTP2 ###
################

samples.trt.CTP2 <- samples.trt %>% filter(CTP == 2)
samples.trt.CTP2$corticosteroid_1 <- as.factor(samples.trt.CTP2$corticosteroid_1)
rownames(samples.trt.CTP2) <- samples.trt.CTP2$sample_id

txi.trt.CTP2 <- txi.trt
txi.trt.CTP2$abundance <- txi.trt.CTP2$abundance[,colnames(txi.trt.CTP2$abundance) %in% samples.trt.CTP2$sample_id]
txi.trt.CTP2$counts <- txi.trt.CTP2$counts[,colnames(txi.trt.CTP2$counts) %in% samples.trt.CTP2$sample_id]
txi.trt.CTP2$length <- txi.trt.CTP2$length[,colnames(txi.trt.CTP2$length) %in% samples.trt.CTP2$sample_id]

all(rownames(samples.trt.CTP2) == colnames(txi.trt.CTP2$abundance))

ddstxi.trt <- DESeqDataSetFromTximport(txi.trt.CTP2, colData=samples.trt.CTP2, design= ~corticosteroid_1)

## Prefiltering
dim(ddstxi.trt)
keep<-rowSums(counts(ddstxi.trt))>=ncol(ddstxi.trt)
ddstxi.trt<-ddstxi.trt[keep,]
dim(ddstxi.trt)

#### GSEA pathways for CTP2

# formatting data from deseq2
# Rereading dds to work with ENSEMBL IDs, as to avoid duplicates

dds.gsea <- DESeq(ddstxi.trt)
resultsNames(dds.gsea)
res.gsea <- results(dds.gsea, name = "corticosteroid_1_1_vs_0")
res.gsea <- as.data.frame(res.gsea)

original_gene_list <- res.gsea$log2FoldChange
names(original_gene_list) <- rownames(res.gsea)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
length(gene_list)

ctp2.up<-rownames(res.gsea[res.gsea$pvalue<0.05 & res.gsea$log2FoldChange>0,])
ctp2.down<-rownames(res.gsea[res.gsea$pvalue<0.05 & res.gsea$log2FoldChange<0,])


#GSE
gse.CTP2 <- gseGO(geneList=gene_list, 
                  ont ="BP", 
                  keyType = "ENSEMBL", 
                  minGSSize = 10, 
                  maxGSSize = 800, 
                  pvalueCutoff = 1, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")


#####################################################
## Euler diagram

library(eulerr)
steroids.up<-c("CTP1.up"=length(ctp1.up),
               "CTP2.up"=length(ctp2.up),
               "CTP1.up&CTP2.up"=length(intersect(ctp1.up, ctp2.up)))

steroids.down<-c("CTP1.down"=length(ctp1.down),
                 "CTP2.down"=length(ctp2.down),
                 "CTP1.down&CTP2.down"=length(intersect(ctp1.down, ctp2.down)))

plot(euler(c(steroids.up, steroids.down)),
     quantities=TRUE, 
     fill=rep(c(cluster1, cluster2),2),
     edges=FALSE,
     labels=c("CTP1", "CTP2", "CTP1", "CTP2"))

########################################################
## Plot

gse.CTP1 <- pairwise_termsim(gse.CTP1)
gse.CTP2 <- pairwise_termsim(gse.CTP2)

gse.CTP1.ids <- gse.CTP1@result[,1:10] %>% filter(p.adjust < 0.05) %>% dplyr::select(ID) %>% pull
gse.CTP2.ids <- gse.CTP2@result[,1:10] %>% filter(p.adjust < 0.05) %>% dplyr::select(ID) %>% pull

ids <- c(gse.CTP1.ids, gse.CTP2.ids)
ids <- unique(ids)
length(ids)

gse.CTP1.cut <- gse.CTP1@result[gse.CTP1@result$ID %in% ids, 1:10] %>% 
  arrange(ID) %>% 
  mutate(CTP = "CTP1")
gse.CTP2.cut <- gse.CTP2@result[gse.CTP2@result$ID %in% ids, 1:10] %>% 
  arrange(ID) %>% 
  mutate(CTP = "CTP2")

gse <- rbind(gse.CTP1.cut, gse.CTP2.cut) %>% arrange(CTP, desc(enrichmentScore))
gse$Description <- factor(gse$Description, levels = unique(gse$Description))

divergent_GOs<-tapply(gse$enrichmentScore, gse$ID, prod)
divergent_GOs<-names(divergent_GOs[divergent_GOs<0])
keywords<-"T cell|B cell|interleukin|STAT|JAK|immune"
redundant<-"positive|negative|peptidyl|recombination"

ggplot(gse[gse$ID %in% divergent_GOs & grepl(keywords, gse$Description) & !grepl(redundant, gse$Description),], aes(x = enrichmentScore, y = Description, color = p.adjust, size = setSize))+
  geom_point()+
  geom_vline(xintercept=0, linetype=2)+
  theme_dose(font.size=12)+
  #theme(legend.position="top")+
  facet_wrap(~CTP)+
  scale_color_gradientn(colours = sig_colors,
                        trans="log",
                        limits= c(1e-8, 0.05), 
                        breaks=c(1e-8,1e-5, 1e-3, 0.05) ) +
  labs(size="Number of genes", color="Adjusted p value")





