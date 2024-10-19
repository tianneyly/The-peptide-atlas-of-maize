library(ggplot2)
library(hrbrthemes)
library(ggsci)

#-------------------------------------------------------------------------------
#cp2cp cp2N
data <- read.table("./density1.txt",header = T,check.names = F)

ggplot(data=data, aes(x=R,fill=Class,color=Class)) +
  geom_density(adjust=1.5, alpha=0.5,lwd=0.5,linetype = 1) +
  theme_bw()+theme(panel.grid=element_blank())

#-------------------------------------------------------------------------------
#cp2cp cp2p
#cp2cp cp2N
data <- read.table("./density2.txt",header = T,check.names = F)

ggplot(data=data, aes(x=R,fill=Class,color=Class)) +
  geom_density(adjust=1.5, alpha=0.5,lwd=0.5,linetype = 1) +
  theme_bw()+theme(panel.grid=element_blank())


#-------------------------------------------------------------------------------
#box plot

data_box <- read.table("Boxplot.txt",header = T,check.names = F)

p1 <- ggplot(data_box,aes(x=D,y=R,fill=G))+
  geom_boxplot(width=0.6,alpha=0.8,notch = TRUE)+
  theme_bw()+theme(panel.grid=element_blank())+
  scale_y_continuous(limits = c(-0.5, 1),breaks = seq(-0.5,1,0.5))
p1


#-------------------------------------------------------------------------------
#KEGG Enrich

remotes::install_github("YuLab-SMU/createKEGGdb")
createKEGGdb::create_kegg_db('zma')

file="KEGG_cp2n.txt"

sn=gsub(".txt","",file)
gene <- read.delim(file,header = T,row.names = 1)
kk <- enrichKEGG(gene= gene$GenBank.Gene,organism= 'zma',
                 pvalueCutoff = 0.743,
                 qvalueCutoff = 1,
                 use_internal_data =T)
head(kk)
write.csv(kk@result,file = paste0(sn,".pathway.csv"),row.names = F)
##只展示过阈值的pathway
kkp <- enrichKEGG(gene= gene$GenBank.Gene,organism= 'zma',qvalueCutoff = 0.1,use_internal_data =T)
head(kkp)
p=dotplot(kk);p
ggsave(p,filename = paste0(sn,".pathway.png"),width = 8,height = 7,dpi = 300)
ggsave(p,filename = paste0(sn,".pathway.pdf"),width = 8,height = 7,dpi = 300)

