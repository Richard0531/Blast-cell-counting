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
library(tiff)
library(grid)
#library(leiden)
#library(ggpubr)
args <- commandArgs(trailingOnly = TRUE)
SAMPLE1 = args[1]
SAMPLE2 = args[2]
#INPUT <- args[2]
OUTPUT <- args[3]
USER <- args[4]
OUTPUT = '/home/yhsiao/pharmaco_imaging/Cellprofiler/blast_output/validation_manuallyFocus/0813_2/'
SAMPLE1 = '129805'
SAMPLE2 = '129758'
USER = 'Richard'
setwd('/home/yhsiao/pharmaco_imaging/Cellprofiler/')
df_test = as.data.frame(read_csv(paste0(OUTPUT,'blast_Image.csv')))
df_test1 = subset(df_test, grepl(SAMPLE1, FileName_Blue))
df_test2 = subset(df_test, grepl(SAMPLE2, FileName_Blue))

total_count1 = sum(df_test1$Classify_Count_negative) + sum(df_test1$Classify_Count_positive) + sum(df_test1$Classify_Count_removed)
blast_count1 = sum(df_test1$Classify_Count_positive) 
blast_percentage1 = round(blast_count1 / total_count1,2)
total_count2 = sum(df_test2$Classify_Count_negative) + sum(df_test2$Classify_Count_positive) + sum(df_test2$Classify_Count_removed)
blast_count2 = sum(df_test2$Classify_Count_positive) 
blast_percentage2 = round(blast_count2 / total_count2,2)

tiff_files1 <- list.files(path = OUTPUT, pattern = "\\.tiff$", full.names = TRUE)
tiff_files1 = grep(SAMPLE1, tiff_files1, value = TRUE)
tiff_files2 <- list.files(path = OUTPUT, pattern = "\\.tiff$", full.names = TRUE)
tiff_files2 = grep(SAMPLE2, tiff_files2, value = TRUE)
#img1 <- readTIFF(tiff_files[1])
#img2 <- readTIFF(tiff_files[1])
#plot(as.raster(img))
pdf(paste0(OUTPUT,SAMPLE1,'_',SAMPLE2,'_report.pdf'),  width = 7, height = 7)
plot.new()
title( paste0('Blast Counting Report'), line = -2, cex.main = 1.5)
text(0.05, 0.85, paste0('User: ',USER), adj = 0)
text(0.05, 0.8, paste0('Time: ',format(Sys.time(), "%Y-%m-%d ")), adj = 0)

text(0.05, 0.7, paste0("Sample 1 ID: ",SAMPLE1), adj = 0)
text(0.05, 0.65, paste0('Total Cell Counts: ',total_count1), adj = 0)
text(0.05, 0.6, paste0('Blast Cell Counts: ',blast_count1), adj = 0)
text(0.05, 0.55, paste0('Blast Percentage: ',blast_percentage1), adj = 0)

text(0.05, 0.45, paste0("Sample 2 ID: ",SAMPLE2), adj = 0)
text(0.05, 0.4, paste0('Total Cell Counts: ',total_count2), adj = 0)
text(0.05, 0.35, paste0('Blast Cell Counts: ',blast_count2), adj = 0)
text(0.05, 0.3, paste0('Blast Percentage: ',blast_percentage2), adj = 0)

text(0.05, 0.2, paste0('Comments: '), adj = 0)

index_order <- as.vector(rbind(seq(2, 6, by = 2), seq(1, 6 - 1, by = 2)))
for (i in index_order) {
  img <- readTIFF(tiff_files1[i])
  grid.newpage()
  img_raster <- as.raster(img)
  grid.raster(img_raster, width = unit(1, "npc"), height = unit(0.8, "npc"))
  
  # Add title
  grid.text(paste0(SAMPLE1), 
            x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
}
for (i in index_order) {
  img <- readTIFF(tiff_files2[i])
  grid.newpage()
  img_raster <- as.raster(img)
  grid.raster(img_raster, width = unit(1, "npc"), height = unit(0.8, "npc"))
  
  # Add title
  grid.text(paste0(SAMPLE2), 
            x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
}
dev.off()

































