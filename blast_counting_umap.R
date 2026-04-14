
#library(pals)
library(uwot)
library(tidyverse)
library(dplyr)
#library(sva)
library(colorspace)
library(readxl)
library(ggplot2)
library(caret)
library(igraph)
library(RANN)
library(ggprism)
library(Seurat)
#library(leiden)
#library(ggpubr)
args <- commandArgs(trailingOnly = TRUE)
SAMPLE = args[1]
#INPUT <- args[2]
OUTPUT <- args[2]
OUTPUT = '/home/yhsiao/pharmaco_imaging/Cellprofiler/blast_output/pipeline_test4/'
SAMPLE = '64'
setwd('/home/yhsiao/pharmaco_imaging/Cellprofiler/')
uu <- load_uwot("/home/yhsiao/pharmaco_imaging/Cellprofiler/ref_UMAP.rds")
mad = read_csv('/home/yhsiao/pharmaco_imaging/Cellprofiler/mads.csv')
mads = mad$mads
names(mads) = mad$columns
df_ref = as.data.frame(read_csv('/home/yhsiao/pharmaco_imaging/Cellprofiler/blast_output/0518_cropped_classified/blast_IdentifyPrimaryObjects.csv'))
objectID = paste0(df_ref$ImageNumber,'_',df_ref$ObjectNumber)
rownames(df_ref) = df_ref$objectID
df_ref = df_ref[,colnames(df_ref) %in% names(mads)]
df_ref = cbind(objectID,df_ref)
metadata_ref = read_csv('/home/yhsiao/pharmaco_imaging/Cellprofiler/metadata_ref.csv')
metadata_ref$group = 'ref'
removed_ref = metadata_ref[metadata_ref$Classify_Class == 'removed' ,]$objectID
rownames(df_ref) = df_ref$objectID
#df_test_feat = df_test[,-c(2,3,112,113,114,115,116,189,190,191,204,205,302)]



#ref.exprs <- scale(df_ref[!rownames(df_ref) %in% removed_ref, ][,-c(1,2,3,112,113,114,115,116,189,190,191,204,205,302)])
#ref.exprs = ref.exprs[,order(-mads)[1:200] ]
#ref.subtypes <- metadata_ref[metadata_ref$Classify_Class != 'removed' , "Classify_Class"]
#ptrain.overall <- pamr.train(list(x = as.matrix(t(ref.exprs)), y = ref.subtypes), n.threshold = 200)
#saveRDS(ptrain.overall, file = "/home/yhsiao/pharmaco_imaging/Cellprofiler/ptrain.overall.RDS")
#ptrain.overall <- readRDS("/home/yhsiao/pharmaco_imaging/Cellprofiler/ptrain.overall.RDS")


df_test = as.data.frame(read_csv(paste0(OUTPUT,'blast_IdentifyPrimaryObjects.csv')))
#nrow_test = nrow(df_test)
objectID = paste0(df_test$ImageNumber,'_',df_test$ObjectNumber,'_','test')
metadata_test = df_test[,c(1,2,111,203,204)]
metadata_test$sample = metadata_test$ImageNumber
metadata_test$Anno_Shu =''
metadata_test = cbind(objectID,metadata_test)
rownames(metadata_test) = metadata_test$objectID
metadata_test$group = 'target'
df_test = df_test[,colnames(df_test) %in% names(mads)]
df_test = cbind(objectID,df_test)
rownames(df_test) = df_test$objectID

#df_ref_test = df_ref[1:nrow_test,]
df_test = rbind(df_ref,df_test)
#rownames(df_test_feat) = df_test_feat$objectID
#metadata_ref_test = metadata_ref[1:nrow_test,]
metadata_test = rbind(metadata_ref,metadata_test)

removed = metadata_test[metadata_test$Classify_Class == 'removed' ,]$objectID
metadata_test <- metadata_test[metadata_test$Classify_Class != 'removed' , ]
metadata_test$Classify_Class[metadata_test$Classify_Class == "erythrocyte"] <- "Erythrocyte"
metadata_test$Classify_Class[metadata_test$Classify_Class == "positive"] <- "Blast"
metadata_test$Classify_Class[metadata_test$Classify_Class == "negative"] <- "Other"
total_count = nrow(metadata_test[metadata_test$group == 'target',][metadata_test[metadata_test$group == 'target',]$Classify_Class == 'Blast' | metadata_test[metadata_test$group == 'target',]$Classify_Class == 'Other'| metadata_test[metadata_test$group == 'target',]$Classify_Class == 'removed',])
blast_count = nrow(metadata_test[metadata_test$group == 'target',][metadata_test[metadata_test$group == 'target',]$group == 'target',][metadata_test[metadata_test$group == 'target',]$Classify_Class == 'Blast',])
other_count = nrow(metadata_test[metadata_test$group == 'target',][metadata_test[metadata_test$group == 'target',]$Classify_Class == 'Other',])
blast_percentage = round(blast_count / total_count,2)
#removed = metadata_test[!metadata_test$Anno_Shu %in% c('blast','erythrocyte','others','removed') ,]$objectID
#metadata_test <- metadata_test[!metadata_test$objectID %in% removed , ]

df_test <- df_test[!rownames(df_test) %in% removed, ]
df_test <- as.data.frame(scale(df_test[,-1]))
#df_cluster = df_test[,order(-mads)[1:200]]
#set.seed(123)
#k <- 3
#km_result <- kmeans(df_cluster, centers = k)
#metadata_test$cluster <- as.factor(km_result$cluster)
###
#seurat_obj <- CreateSeuratObject(counts = t(df_test[,order(-mads)[1:200]]))  # transpose if rows = samples
#seurat_obj <- FindVariableFeatures(seurat_obj)
#seurat_obj <- ScaleData(seurat_obj)
#seurat_obj <- RunPCA(seurat_obj)
#seurat_obj <- RunPCA(seurat_obj, features = NULL)
#seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
#seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
#df$cluster <- Idents(seurat_obj)
###

#test_id = metadata_test[metadata_test$group == 'target',]$objectID
#pred.overall <- as.data.frame(pamr.predict(ptrain.overall, as.matrix(t(df_test_feat[,order(-mads)[1:200] ])), 8, type  = "posterior"))
#pred.overall = round(pred.overall[rownames(pred.overall) %in% test_id,],3)
#pred.overall$pred <- names(pred.overall)[max.col(pred.overall[, 1:3], ties.method = "first")]
#pred.overall$objectID = rownames(pred.overall)
#metadata_test = dplyr::left_join(metadata_test,pred.overall,by='objectID')
#mads <- apply(df_test,2,mad)
#set.seed(123)
#uu.test <- umap(df_test[,order(-mads)[1:250] ], metric  = "euclidean", spread = 1,min_dist = 0.00001, ret_model = TRUE,  n_neighbors = 15)
#uu.test = as.data.frame(uu.test$embedding)

uu.test <- as.data.frame(umap_transform(df_test[,order(-mads)[1:200] ], model = uu))
uu.test$celltype = metadata_test$Classify_Class
uu.test$group = metadata_test$group
uu.test$PointSize <- ifelse(uu.test$group == "target", 2,1) 
uu.test$Pred = metadata_test$Anno_Shu
uu.test$cluster = 'Erythrocyte'
uu.test[uu.test$V1 < (-3.5),]$cluster = 'Other'
uu.test[uu.test$V1 > (-3.5),][uu.test[uu.test$V1 > (-3.5),]$V2 > (-1),]$cluster = 'Blast'
#uu.test$cluster = as.factor(uu.test$cluster)     
celltype_colors <- c('#C71585','#000080','#00CED1','#B0C4DE')
celltype_colors_named <- setNames(celltype_colors, levels(as.factor(uu.test$celltype)))
cluster_colors <- c('#C71585','#000080','#00CED1','#B0C4DE')
cluster_colors_named <- setNames(cluster_colors, levels(as.factor(uu.test$cluster)))

p = ggplot(uu.test, aes(x = V1, y = V2)) +
  geom_point(aes(shape = group, color = celltype, size = PointSize)) +
  scale_shape_manual(values = c("ref" = 16, "target" = 17)) +
  scale_size_identity() + 
  scale_color_manual(values = celltype_colors_named)+
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  #geom_vline(xintercept=-3.5, color = "black")+
  #geom_segment(aes(x = -3.5, y = -1, xend = 5, yend = -1))+
  ggtitle(paste0('Sample: ',SAMPLE, '\n','Blast Counts: ',blast_count,'\n','Total Counts: ',total_count,'\n','Blast Percentage: ',blast_percentage ))
p
#ggsave(paste0(OUTPUT,SAMPLE,"_BlastCounting_Report.pdf"), plot = p, width = 6, height = 6)


df_test_rank = as.data.frame(apply(df_test, 2, rank))
uu.test.feat = cbind(uu.test,df_test_rank)
p = ggplot(uu.test, aes(x = V1, y = V2)) +
  geom_point(aes(shape = group, color = Pred, size = PointSize)) +
  scale_shape_manual(values = c("ref" = 16, "target" = 17)) +
  scale_size_identity() + 
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  ggtitle(paste0('Sample: ',SAMPLE, '\n','Blast Counts: ',blast_count,'\n','Total Counts: ',total_count,'\n','Blast Percentage: ',blast_percentage ))
p


total_count = nrow(uu.test[uu.test$group == 'target',][uu.test[uu.test$group == 'target',]$cluster == 'Blast' | uu.test[uu.test$group == 'target',]$cluster == 'Other',])
blast_count = nrow(uu.test[uu.test$group == 'target',][uu.test[uu.test$group == 'target',]$cluster == 'Blast',])
other_count = nrow(uu.test[uu.test$group == 'target',][uu.test[uu.test$group == 'target',]$cluster == 'Other',])
blast_percentage = round(blast_count / total_count,2)
p = ggplot(uu.test, aes(x = V1, y = V2)) +
  geom_point(aes(shape = group, color = cluster, size = PointSize)) +
  scale_shape_manual(values = c("ref" = 16, "target" = 17)) +
  scale_size_identity() + 
  scale_color_manual(values = cluster_colors)+
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  ggtitle(paste0('Sample: ',SAMPLE, '\n','Blast Counts: ',blast_count,'\n','Total Counts: ',total_count,'\n','Blast Percentage: ',blast_percentage ))
p
###
blast_list = rownames(uu.test[uu.test$cluster == 'Blast',])
df_blast = df_test[rownames(df_test) %in% blast_list,]
meta_blast = meta





###

pred_colors <- c('#C71585','#000080','#00CED1','#B0C4DE')
pred_colors_named <- setNames(celltype_colors, levels(as.factor(uu.test$Pred)))


total_count = nrow(metadata_test[metadata_test$group == 'target',][metadata_test[metadata_test$group == 'target',]$pred == 'positive' | metadata_test[metadata_test$group == 'target',]$pred == 'negative',])
blast_count = nrow(metadata_test[metadata_test$group == 'target',][metadata_test[metadata_test$group == 'target',]$group == 'target',][metadata_test[metadata_test$group == 'target',]$pred == 'positive',])
other_count = nrow(metadata_test[metadata_test$group == 'target',][metadata_test[metadata_test$group == 'target',]$pred == 'negative',])
blast_percentage = round(blast_count / total_count,2)

cluster_colors <- c('#FF8C00','#8B0000','#6495ED','#008000','#8B4513',
                   '#C71585','#000080','#B0C4DE','#00CED1','#6A5ACD','#FFD700')
cluster_colors_named <- setNames(cluster_colors, levels(as.factor(uu.test$cluster)))

p = ggplot(uu.test, aes(x = V1, y = V2)) +
  geom_point(aes(shape = group, color = cluster)) +
  scale_shape_manual(values = c("ref" = 16, "target" = 17)) +
  scale_size_identity() + 
  scale_color_manual(values = cluster_colors_named)+
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  ggtitle(paste0('Sample: ',SAMPLE, '\n','Blast Counts: ',blast_count,'\n','Total Counts: ',total_count,'\n','Blast Percentage: ',blast_percentage ))
p












###Umap model training

setwd('/home/yhsiao/pharmaco_imaging/Cellprofiler/blast_output/')
shu_anno = as.data.frame(read_csv('annotation_Shu_NN.csv'))
shu_anno$objectID = paste0(shu_anno$ImageNO,'_',shu_anno$Cell)
shu_anno = shu_anno[,c(4,7)]


df_all = as.data.frame(read_csv('/home/yhsiao/pharmaco_imaging/Cellprofiler/blast_output/0518_cropped_classified/blast_IdentifyPrimaryObjects.csv'))
df_all <- df_all[, !grepl("_X", names(df_all))]
df_all <- df_all[, !grepl("_Y", names(df_all))]
objectID = paste0(df_all$ImageNumber,'_',df_all$ObjectNumber)
df_all = cbind(objectID,df_all)
metadata = df_all[,c(1,2,3,106,180,181)]
metadata$sample = metadata$ImageNumber
metadata$Classify_Class[metadata$Classify_Class == "erythrocyte"] <- "Erythrocyte"
metadata$Classify_Class[metadata$Classify_Class == "positive"] <- "Blast"
metadata$Classify_Class[metadata$Classify_Class == "negative"] <- "Other"
metadata = left_join(metadata,shu_anno,by = 'objectID')
df = df_all[,-c(2,3,106,107,108,109,110,170:181)]
rownames(df) = df$objectID
removed = metadata[metadata$Classify_Class == 'removed'  ,]$objectID
df <- df[!rownames(df) %in% removed, ]
df <- as.data.frame(scale(df[,-1]))
mads <- apply(df,2,mad)
set.seed(123)
uu <- umap(df[,order(-mads)[1:200] ], metric  = "euclidean", spread = 1,min_dist = 0.00001, ret_model = TRUE,  n_neighbors = 50)
umap_embedding = as.data.frame(uu$embedding)

umap_embedding$celltype = metadata[metadata$Classify_Class != 'removed' ,]$Classify_Class
#umap_embedding$celltype = metadata$Classify_Class
celltype = unique(metadata$Classify_Class)
celltype_colors <- c('#C71585','#000080','#00CED1','#FFD700')
celltype_colors_named <- setNames(celltype_colors, levels(as.factor(umap_embedding$celltype)))
total_count = nrow(metadata[metadata$Classify_Class == 'Blast' | metadata$Classify_Class == 'Other',])
blast_count = nrow(metadata[metadata$Classify_Class == 'Blast',])
other_count = nrow(metadata[metadata$Classify_Class == 'Other',])
blast_percentage = round(blast_count / total_count,2)


p = ggplot(umap_embedding, aes(x = V1, y = V2, color = cluster)) +
  geom_point(size = 1) +
  scale_color_manual(values = cluster_colors_named)+
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  geom_vline(xintercept=-2, color = "black")+
  geom_hline(yintercept=-2, color = "black")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
p
umap_embedding$cluster = 'Erythrocyte'
umap_embedding[umap_embedding$V1 < (-2),]$cluster = 'Other'
umap_embedding[umap_embedding$V1 > (-2),][umap_embedding[umap_embedding$V1 > (-2),]$V2 > (-2),]$cluster = 'Blast'
blast_list = rownames(umap_embedding[umap_embedding$cluster == 'Blast',])
df_blast = df[rownames(df) %in% blast_list,]
meta_blast = metadata[metadata$objectID %in% blast_list,]
rownames(meta_blast) = meta_blast$objectID
set.seed(123)
mads_blast <- apply(df_blast,2,mad)
uu_blast <- umap(df_blast[,order(-mads_blast)[1:250] ], metric  = "correlation", spread = 0.8,min_dist = 0.005, ret_model = TRUE)
umap_embedding_blast = as.data.frame(uu_blast$embedding)
umap_embedding_blast$celltype = meta_blast$Classify_Class
p = ggplot(umap_embedding_blast, aes(x = V1, y = V2, color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = celltype_colors_named)+
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
p

#save_uwot(uu, file = "/home/yhsiao/pharmaco_imaging/Cellprofiler/ref_UMAP.rds")
mad = as.data.frame(mads)
mad$columns = rownames(mad)
write_csv(mad,'/home/yhsiao/pharmaco_imaging/Cellprofiler/mads.csv')
#ggsave(paste0(OUTPUT,SAMPLE,"_BlastCounting_Report.pdf"), plot = p, width = 8, height = 6)

mad = read_csv('/home/yhsiao/pharmaco_imaging/Cellprofiler/mads.csv')
mads = mad$mads
names(mads) = mad$columns
write_csv(metadata,'/home/yhsiao/pharmaco_imaging/Cellprofiler/metadata_ref.csv')

###
metadata_ref = read_csv('/home/yhsiao/pharmaco_imaging/Cellprofiler/metadata_ref.csv')
metadata_ref$group = 'ref'
removed_ref = metadata_ref[metadata_ref$Classify_Class == 'removed' ,]$objectID
umap_embedding = as.data.frame(uu$embedding)

umap_embedding$celltype = metadata_ref[metadata_ref$Classify_Class != 'removed' ,]$Classify_Class
celltype = unique(metadata_ref$Classify_Class)
celltype_colors <- c('#C71585','#000080','#00CED1','#FFD700')
celltype_colors_named <- setNames(celltype_colors, levels(as.factor(umap_embedding$celltype)))
total_count = nrow(metadata_ref[metadata_ref$Classify_Class == 'Blast' | metadata_ref$Classify_Class == 'Other',])
blast_count = nrow(metadata_ref[metadata_ref$Classify_Class == 'Blast',])
other_count = nrow(metadata_ref[metadata_ref$Classify_Class == 'Other',])
blast_percentage = round(blast_count / total_count,2)


p = ggplot(umap_embedding, aes(x = V1, y = V2, color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = celltype_colors_named)+
  theme_void()+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  ggtitle(paste0('Sample: ','Reference', '\n','Blast Counts: ',blast_count,'\n','Total Counts: ',total_count,'\n','Blast Percentage: ',blast_percentage ))
p


###
k <- 10
df_nn = df_test[,order(-mads)[1:200] ]
nn <- nn2(df_nn, k = k + 1)  # include self
edges <- data.frame(
  from = rep(1:nrow(df_nn), each = k),
  to = as.vector(nn$nn.idx[, 2:(k+1)])  # skip self (1st col)
)

# 🔹 2. Build undirected graph and cluster
g <- graph_from_data_frame(edges, directed = FALSE)
g <- simplify(g)
cl <- cluster_louvain(g)
clusters <- membership(cl)


p = ggplot(uu.test, aes(x = V1, y = V2)) +
  geom_point(aes(shape = group, color = celltype, size = PointSize)) +
  scale_shape_manual(values = c("ref" = 16, "target" = 17)) +
  scale_size_identity() + 
  scale_color_manual(values = celltype_colors_named)+
  guides(color = guide_legend(override.aes = list(size = 5)))
p

test1 = uu.test[uu.test$V1< (-3.75) & uu.test$celltype == 'Blast', ]
test2 = uu.test[uu.test$V1 > (-3.75) & uu.test$V2 > (-1.25) & uu.test$celltype != 'Blast', ]

















