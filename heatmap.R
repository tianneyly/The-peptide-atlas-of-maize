library("pheatmap")
library("RColorBrewer")
library("dplyr")

data <- read.csv(file="cp_tf.csv",
                 header=T)
df1 <- data[!duplicated(data),]
head(df1)
rownames(df1) <- df1[,1]
df2 <- df1[,-1]
head(df2)


data1 <- df2[,c(4:16)]

# 热图A
pheatmap(data1,scale='row',show_rownames = F,
         color=colorRampPalette(c("black","#EEE53B"))(100),
         cluster_rows = F, cluster_cols = F,
         legend=F,
         fontsize_col=12,
         )


# # 热图B
# 
data_B <- read.csv(file="cp_tf_p.csv",
                   header=T,row.name=1)
data_B1 <- data_B[,c(1:13)]
# 
pheatmap(data_B1,show_rownames = T,
         # scale='row',
         color=colorRampPalette(c("white","red"))(100),
         cluster_rows = F, cluster_cols = F,
         border_color = 'black',
         legend=T,
         legend_breaks = NA,
         legend_labels = "-log10 P-value",
         fontsize_col=12,
         cellwidth = 10, cellheight = 8
         )

# 热图C

data_C <- read.csv(file="cp_tf_avg.csv",
                   header=T,row.name=1)
data_C1 <- data_C[,c(1:13)]

pheatmap(data_C1,show_rownames = T,
         scale='row',
         color=colorRampPalette(c("blue","white","red"))(100),
         cluster_rows = F, cluster_cols = F,
         border_color = 'black',
         legend=T,
         fontsize_col=12,
         cellwidth = 10,cellheight = 8,
         )


