library(TCseq)


Idents(obj) <-"Condition"
health <- subset(obj, idents ="healthy")
Idents(severe) <- 'Time'

each_timepoint_list <- list()

each_timepoint <- unique(severe$Time)


for (RN in each_timepoint) {
  
  each_timepoint_list[[RN]] <- subset(severe , idents = RN)
  
}
names(each_timepoint_list) <- each_timepoint



d3<-AverageExpression(each_timepoint_list[['d3']], assay = 'RNA', layer = 'data')
d1<-AverageExpression(each_timepoint_list[['d1']], assay = 'RNA', layer = 'data')
d7<-AverageExpression(each_timepoint_list[['d7']], assay = 'RNA', layer = 'data')
d1<-d1$RNA
d3<-d3$RNA
d7<-d7$RNA
avg_exp <- AverageExpression(severe, assays = "RNA" , layer = "data")
x <- as.data.frame(avg_exp$RNA)
high_expressed_genes <- x[apply(x, 1, max) >= 0.1, ]
ref_row_names <- rownames(high_expressed_genes)
d1_reordered <- d1[ref_row_names, , drop = TRUE]
d3_reordered <- d3[ref_row_names, , drop = TRUE]
d7_reordered <- d7[ref_row_names, , drop = TRUE]

combined_matrix <- cbind(d1_reordered, d3_reordered, d7_reordered)

# Optionally, you can rename the columns of the combined matrix
colnames(combined_matrix) <- c('d1', 'd3', 'd7')
normalized_matrix <- apply(combined_matrix, 1, function(x) (x - min(x)) / (max(x) - min(x)))

data <- as.matrix(combined_matrix)
data <- log1p(combined_matrix)
kmeans_data <-
  timeclust(x=data,
            algo = 'cm',  #指定聚类方法，选项包括 'km' (kmeans;k均值划分)、'pam'（围绕中心点划分）、'hc'（层次聚类）、'cm' (cmeans)
            k =9, 
            standardize = TRUE , 
  )

pdf('severe.pdf',width = 20, height = 15)
p <-
  timeclustplot(
    kmeans_data, 
    categories="timepoint",  
    value ="z-score",   
    membership.color =colorRampPalette(c("#1a535c","#4ecdc4","#f7fff7","#ffe66d","#ff6b6b"))(50), #自定义 颜色
    cols = 3, 
    cl.color="black"
  )

dev.off()

# within group variance
gene_names_cluster_4 <- names(kmeans_data@cluster[kmeans_data@cluster == 4])
gene_names_cluster_8 <- names(kmeans_data@cluster[kmeans_data@cluster == 8])
change_ranges <- apply(data[gene_names_cluster_8, ], 1, function(x) max(x) - min(x))
change_ranges[change_ranges>0.3]



