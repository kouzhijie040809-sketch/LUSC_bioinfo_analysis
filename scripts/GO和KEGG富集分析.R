if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("ggtree")

if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("enrichplot")

install.packages("ggplot2")
install.packages("stringi")
install.packages("GOplot")
library(R.utils)

install.packages("tidytree")
install.packages("ggiraph")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringi)
library(GOplot)
R.utils::setOption("clusterProfiler.download.method", "auto")

setwd("D:\\sdu\\bioinfo\\TCGA_LUSC_analysis\\data\\GDC_data")

# 读入
input_diff = read.table("TCGA.diff.wilcoxon.txt", sep="\t", header=T, check.names=F)
input_gene <- input_diff[,1]
# 去除重复基因
input_gene=unique(as.vector(input_gene))

# 将 gene symbol 转换为基因 id
#org.Hs.eg为人的物种
#https://www.jianshu.com/p/84e70566a6c6
entrezIDs=mget(input_gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)

entrezIDs = unlist(entrezIDs)

entrezIDs=as.character(entrezIDs)
# 去除基因 id 为 NA 的基因
gene=entrezIDs[entrezIDs!="NA"]

#筛选
pvalueFilter=0.05
qvalueFilter=1

if(qvalueFilter>0.05){colorSel="pvalue"}else{colorSel="qvalue"}

# GO
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

# 保存
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

# 显示的数目
showNum=10

#GO柱状图
pdf(file = "GO_barplot.pdf", width = 10, height = 20)

bar <- barplot(
  kk,
  drop = TRUE,
  showCategory = showNum,
  label_format = 15,
  split = "ONTOLOGY",
  color = colorSel
) + facet_grid(ONTOLOGY ~ ., scale = "free")

print(bar)

dev.off()

#GO气泡图
pdf("GO_bubble.pdf", width = 10, height = 20)

bub <- dotplot(
  kk,
  showCategory = showNum,
  orderBy = "GeneRatio",
  label_format = 15,
  split = "ONTOLOGY",
  color = colorSel
) + facet_grid(ONTOLOGY ~ ., scales = "free")

print(bub)

dev.off()


#筛选
pvalueFilter=0.05
qvalueFilter=1

if(qvalueFilter>0.05){colorSel="pvalue"}else{colorSel="qvalue"}

#KEGG

kk <- enrichKEGG(
  gene = gene,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)


KEGG=as.data.frame(kk)

KEGG <- KEGG[(KEGG$pvalue < pvalueFilter & KEGG$qvalue < qvalueFilter), ]

KEGG$geneID <- sapply(strsplit(KEGG$geneID, "/"), function(x) {
  paste(input_gene[match(x, as.character(entrezIDs))], collapse = "/")
})


write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

# 柱状图
pdf(file="KEGGbarplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

# 气泡图
pdf(file="KEGGbubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()







