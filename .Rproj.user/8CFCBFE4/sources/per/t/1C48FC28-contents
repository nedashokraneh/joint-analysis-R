LINE_segway_annotation <- read.csv("../../../data/annotations/segway_LINE.bed", sep = "\t", header = FALSE)
get_compartment <- function(chr, pos){
  a = LINE_segway_annotation[LINE_segway_annotation$V1==chr & LINE_segway_annotation$V2 < pos & LINE_segway_annotation$V3 >= pos,4]
  #print(a)
  return(a)
}
##########################
library(umap)
library(ggplot2)
LINE_embedding <- read.csv("../../../data/embedding/indexed_GM12878_line_embedding_HiC_graph_order_3_samples_25M.embedding", sep = "\t", header= FALSE)
LINE_embedding$chr_num = as.numeric(substring(LINE_embedding$V1,4))
LINE_embedding <- LINE_embedding[order(LINE_embedding$chr_num,
                                       LINE_embedding$V2),]
LINE_embedding[,c(3:102)] = scale(LINE_embedding[,c(3:102)])
LINE_kmeans <- kmeans(LINE_embedding[,c(3:102)], centers = 5, nstart = 25)
LINE_one_feature_kmeans <- kmeans(LINE_embedding[,3], centers = 5, nstart = 25)
LINE_embedding$kmean_compartment = LINE_kmeans$cluster
LINE_embedding$one_feature_kmean_compartment = LINE_one_feature_kmeans$cluster
ggplot(LINE_embedding) + geom_point(aes(x=V3,y=V4,color=kmean_compartment))
ggplot(LINE_embedding) + geom_point(aes(x=V3,y=V4,color=one_feature_kmean_compartment))
LINE.umap = umap(LINE_embedding[,c(3:102)])
LINE_embedding[,c("umap1", "umap2")] = LINE.umap$layout
ggplot(LINE_embedding) + geom_point(aes(x=umap1,y=umap2,color=kmean_compartment))
ggplot(LINE_embedding) + geom_point(aes(x=umap1,y=umap2,color=one_feature_kmean_compartment))
ggplot(LINE_embedding[LINE_embedding$V3<=5,]) + geom_point(aes(x=umap1,y=umap2,color=V3))

LINE_KR_embedding <- read.csv("../../../data/embedding/indexed_GM12878_KR_line_embedding_HiC_graph_order_3_samples_25M.embedding", sep = "\t", header= FALSE)
LINE_KR_embedding$chr_num = as.numeric(substring(LINE_KR_embedding$V1,4))
LINE_KR_embedding <- LINE_KR_embedding[order(LINE_KR_embedding$chr_num,
                                       LINE_KR_embedding$V2),]
LINE_KR_embedding[,c(3:102)] = scale(LINE_KR_embedding[,c(3:102)])
LINE_KR_kmeans <- kmeans(LINE_KR_embedding[,c(3:102)], centers = 5, nstart = 25)
LINE_KR_one_feature_kmeans <- kmeans(LINE_KR_embedding[,3], centers = 5, nstart = 25)
LINE_KR_embedding$kmean_compartment = LINE_KR_kmeans$cluster
LINE_KR_embedding$one_feature_kmean_compartment = LINE_KR_one_feature_kmeans$cluster
ggplot(LINE_KR_embedding) + geom_point(aes(x=V3,y=V4,color=kmean_compartment))
ggplot(LINE_KR_embedding) + geom_point(aes(x=V3,y=V4,color=one_feature_kmean_compartment))
LINE_KR.umap = umap(LINE_KR_embedding[,c(3:102)])
LINE_KR_embedding[,c("umap1", "umap2")] = LINE_KR.umap$layout
ggplot(LINE_KR_embedding) + geom_point(aes(x=umap1,y=umap2,color=kmean_compartment))
ggplot(LINE_KR_embedding) + geom_point(aes(x=umap1,y=umap2,color=one_feature_kmean_compartment))
  


LINE_KR_compartments <- read.csv("../../../data/annotations/GM12878_KR_line_embedding_compartments.txt",
                                 header = FALSE, sep = "\t")
#colnames(LINE_KR_compartments) = LINE_KR_compartments[1,]
LINE_KR_compartments <- LINE_KR_compartments[-1,-2]
LINE_KR_compartments[,2] = as.numeric(as.character(LINE_KR_compartments[,2]))/100000
colnames(LINE_KR_compartments) = c("V1", "V2", "sci_compartment")
LINE_KR_embedding = merge(x = LINE_KR_embedding, y = LINE_KR_compartments,
                          by = c("V1", "V2"), all.x = TRUE, all.y = TRUE)
LINE_KR_embedding2 <- na.omit(LINE_KR_embedding[,-103])
LINE_KR.embedding = LINE_KR_embedding2[,c(3:102)]
LINE_KR.embedding = data.frame(scale(LINE_KR.embedding))
LINE_KR.umap = umap(LINE_KR.embedding)
LINE_KR_embedding2[,c("umap1","umap2")] = LINE_KR.umap$layout
ggplot(LINE_KR_embedding2) + geom_point(aes(x=umap1,y=umap2,color=my_kmean_compartment))
i = 1
current_label = -10
while (i <= nrow(LINE_embedding)){
  
}
###########################
LINE_embedding$V1 = as.character(LINE_embedding$V1)
LINE_segway_annotation$V1 = as.character(LINE_segway_annotation$V1)
LINE_embedding$compartment = mapply(get_compartment, LINE_embedding$V1,
                                    LINE_embedding$V2)
library(ggplot2)
ggplot(LINE_embedding) + geom_point(aes(x = V30, y = V42, color = factor(compartment)))


SNIPER_embedding <- read.csv("../../../data/embedding/SNIPER_embedding.txt",
                             header=FALSE, sep = "\t")
SNIPER_cor = cor(SNIPER_embedding[3:131])
SNIPER_annot <- read.csv("../../../data/annotations/segway_SNIPER_LV.bed",
                         header = FALSE, sep = "\t")
get_compartment2 <- function(chr, pos){
  chr = paste("chr",chr,sep = "")
  a = SNIPER_annot[SNIPER_annot$V1==chr & SNIPER_annot$V2 < pos &
                     SNIPER_annot$V3 >= pos,4]
  #print(a)
  return(a)
}
SNIPER_embedding$compartment = mapply(get_compartment2, SNIPER_embedding$V1,
                                    SNIPER_embedding$V2)
SNIPER_embedding$compartment = as.numeric(SNIPER_embedding$compartment)
ggplot(SNIPER_embedding) + geom_point(aes(x = V3, y = V4, color = factor(compartment)))
SNIPER_kmeans <- kmeans(SNIPER_embedding[3:130], centers = 5, nstart = 25)
SNIPER_embedding$kmean_compartment = SNIPER_kmeans$cluster
ggplot(SNIPER_embedding) + geom_point(aes(x = V3, y = V4, color = factor(kmean_compartment)))


LINE_fa_segway <- read.csv("../../../data/annotations/segway_first_LINE_fa.bed",
                           header = FALSE, sep = "\t")
LINE_fa_segway %>% group_by(V4) %>% summarise(sum(V3 - V2))
fa_segway <- read.csv("../../../data/annotations/segway_func_assays.bed",
                      header = FALSE, sep = "\t")

fa_segway %>% group_by(V4) %>% summarise(sum(V3 - V2))

LINE_KR_embedding %>% group_by(kmean_)


LINE_embedding10 <- read.csv("../../../data/embedding/GM12878_KR_line_d10_embedding", sep = " ", header= FALSE)
cor(LINE_embedding10[,2:11])  
LINE_embedding10_1st_order <- read.csv("../../../data/embedding/GM12878_KR_line_d10_embedding_1st_order", sep = " ", header= FALSE)
cor(LINE_embedding10_1st_order[,2:11]) 
LINE_embedding10_2nd_order <- read.csv("../../../data/embedding/GM12878_KR_line_d10_embedding_2nd_order", sep = " ", header= FALSE)
cor(LINE_embedding10_2nd_order[,2:11]) 

LINE_embedding10_1st_order <- read.csv("../../../data/embedding/_order_1_samples_25M.embedding", sep = " ", header= FALSE)
cor(LINE_embedding10_1st_order[,2:11]) 
