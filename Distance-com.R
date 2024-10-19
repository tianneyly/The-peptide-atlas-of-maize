library(ggplot2)
library(dplyr)
library(readr)  # 使用readr代替readxl因为您的文件是CSV

# 读取CSV文件
data1 <- read_csv("disstance-en-spe.csv")

pdf("disstance-en-spe.pdf",width=5,height=4)
p <- ggplot(data1, aes(value, fill=variable)) +
  geom_histogram(data = filter(data1, variable == 'enhanced'),
                 aes(y = after_stat(count)), 
                 binwidth = 30) +
  geom_histogram(data = filter(data1, variable == 'specific'),
                 aes(y = -after_stat(count)), 
                 binwidth = 30) +
  scale_x_continuous(limits = c(0, 8000),
                     breaks = c(0, 2000, 4000, 6000, 8000),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-100, 100),
                     breaks = seq(-100, 100, 20),
                     labels = abs(seq(-100, 100, 20))) +
  scale_fill_manual(values = c('enhanced' = '#88c4e8', 'specific' = '#db6968')) +
  labs(x = 'Distance between adjacent CPs', y = 'Number of CPs') +
  theme_classic(base_size = 15) +
  theme(#panel.border = element_rect(linewidth = 1, fill='transparent'),
        legend.position = 'none',
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),  # 调整标签角度
        axis.title.x = element_text(size = 16, color = "black"),  # 设置X轴标题的字体大小和颜色
                      axis.title.y = element_text(size = 16, color = "black")) +  # 设置Y轴标题的字体大小和颜色
  #geom_vline(xintercept = 2000, linetype = 2, linewidth = 1) +
  #annotate('text', x = 9000, y = 70,
           #label = round(median(filter(data1, variable == 'enhanced')$value), digits = 2),
          # size = 4, color = '#88c4e8') +
  #annotate('text', x = 9000, y = -70,
           #label = round(median(filter(data1, variable == 'specific')$value), digits = 2),
          # size = 4, color = '#db6968') +
  annotate('text', x = 6000, y = 90,
           label = 'Tissue enhanced CPs',
           size = 6, color = '#88c4e8') +
  annotate('text', x = 6000, y = -90,
           label = 'Tissue specific CPs',
           size = 6, color = '#db6968')
print(p) 

# 保存图形为PDF文件 
dev.off()