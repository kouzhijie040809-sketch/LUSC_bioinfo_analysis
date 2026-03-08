install.packages("rjson")
install.packages("tidyverse")

library(tidyverse)
library(rjson)

setwd("D:/sdu/bioinfo/TCGA_LUSC_analysis/data/GDC_data")

#读入meta.data 
metadata <- jsonlite::fromJSON("metadata.cart.2026-03-05.json")

#将样本id从metadata中提取出来  （一个字符向量）
sample_id <- sapply(metadata$associated_entities, function(x) x[,1])

#将文件名与样本id相对应，做成一个数据框
file_sample <- data.frame(sample_id, file_name = metadata$file_name)

#列出TCGA_temp中的tsv文件即获取的数据文件（提取出的count_file是一个文件路径的字符向量）
count_file <- list.files('TCGA_temp', 
                         pattern = '\\.tsv$',          # 严格以 .tsv 结尾
                         recursive = TRUE, 
                         full.names = FALSE)

#将tsv文件名分离并提取出来
count_file_name <- strsplit(count_file, split="/")
count_file_name <- sapply(count_file_name, function(x){x[2]})

#建立一个空数据框
matrix = data.frame(matrix(nrow=60660, ncol=0))


#逐个读取及合并
for (i in 1:length(count_file)){
  path <- paste0('TCGA_temp//',count_file[i])                          #counts文件夹名
  data <- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)  #第一列作行名
  colnames (data)<-data[2,]                                          #第二行作列名
  data <-data[-c(1:6),]                                              #删除前六行   
  #3: unsranded.counts; 4: stranded_first; 5: stranded_second; 6: tpm_unstranded1; 7: fpkm_unstranded; 8: fpkm_uq_unstranded
  #data <- data[3]
   data <- data[6]
  #data <- data[7]
  colnames (data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

  #将gene_name并入
path = paste0('TCGA_temp//',count_file[1])
data <- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <- data[-c(1:6),1]
matrix0<- cbind(gene_name,matrix)
#将gene_type并入
gene_type <- data[-c(1:6),2]
matrix0<- cbind(gene_type,matrix0)

##我们要做基因水平的分析，同一基因会有多个转录本，我们下载的是转录本水平的表达量（transcript counts）所以同一基因会有多行，在接下来的分析中，多行会被当做多个基因，所以我们要删除重复的只保留一个。

#删除重复的基因(这里保留max，min，mean都可以)
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)

#保留mRNA
matrix0 <- subset(x=matrix0,gene_type=="protein_coding")

#将gene_name作为行名，转换为导出格式
rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-c(1,2)]
matrix1 <- data.frame(ID=rownames(matrix0),matrix0)
colnames(matrix1) = gsub('[.]','-',colnames(matrix1))

#导出
write.table(matrix1,'TGCA_LUSC_TPM.txt',sep = "\t",quote = F,row.names = F)















