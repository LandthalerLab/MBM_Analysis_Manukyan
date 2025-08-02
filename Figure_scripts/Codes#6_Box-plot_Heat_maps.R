library(ggplot2)
library(ggpubr)
library(devtools)
library(tidyverse)
library(corrplot)
library("gplots")   
library("colorRamps")

## data
my.data<-read.delim(".....txt", sep="\t",h=T)

## Ordered Box plots (ranked by expression level)
ggplot(my.data, aes(x=reorder(BZW2Level,PDCD1LG2), y=PDCD1LG2)) +
  geom_boxplot(fill="gray")+
  labs(title="",x="BZW2Level", y = "PDCD1LG2")+
  theme_classic()

# Change  automatically color by groups
bp <- ggplot(my.data, aes(x=reorder(BZW2Level, PDCD1LG2), y=PDCD1LG2, fill=BZW2Level)) + 
  geom_boxplot()+
  labs(title="",x="BZW2Level", y = "PDCD1LG2")+
  stat_compare_means(method = "t.test", label.y = 0,size=3)       # Add global annova p-value
  #stat_compare_means(label = "p.signif",label.y=-2.5, size=4, method = "t.test",
                     #ref.group = ".all.")  

bp + theme(legend.position="none")
bp + scale_color_gradient(low = "#0091ff", high = "yellow") + theme_minimal()

##Non-ordered Box plots (ranked by given order)
ggplot(my.data, aes(x=Diet, y="GAPDHS")) + 
  geom_boxplot(fill="gray")+
  labs(title="",x="Diet", y = "GAPDHS")+
  theme_classic()
# Change  automatically color by groups
ggplot(my.data, aes(x=Diet , y=GAPDHS, fill=Diet)) + 
  geom_boxplot()+
  labs(title="",x="", y = "GAPDHS")+
  stat_compare_means(method = "t.test", label.y = 7,size=3)+    
  geom_jitter(width = 0.1)+
  geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
scale_fill_manual(values = c("orange","#0091ff","#CC5CC6", "#005199", "#3DFF26", "#BC00FF")) + theme_minimal()  
  
  
########Generation of Heat map

# data 
my.data<-read.delim("...txt",sep="\t",h=T,row.names = 1, stringsAsFactors = T)

Label = c(rep("lightblue",10),rep("orange",250),rep("yellow",250),
          rep("brown",323))

heatmap.2(as.matrix(my.data), col=my.heat<-colorRampPalette(c("blue","lightgrey","yellow")), scale="col", key=T, keysize=1.5,
          density.info="histogram", trace="none", cexRow=0.9,cexCol=0.8, labRow=rownames(my.data), #das hier funtioniert!
          RowSideColors=Label,lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0), distfun = function(my.data) dist(my.data,method = 'euclidean'),
          lwid=c(1.5,0.2,2.5,2.5))

# Correlation
corr_mydata <- cor(my.data)

corrplot(corr_mydata, method = "pie")

Label = c(rep("lightblue",10),rep("orange",250),rep("yellow",250),
          rep("brown",323))

heatmap.2(as.matrix(corr_mydata), col=my.heat<-colorRampPalette(c("blue","lightgrey","yellow")), scale="col", key=T, keysize=1.5,
          density.info="histogram", trace="none", cexRow=0.9,cexCol=0.8, labRow=rownames(corr_mydata), #das hier funtioniert!
          RowSideColors=Label,lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0), distfun = function(my.data) dist(my.data,method = 'euclidean'),
          lwid=c(1.5,0.2,2.5,2.5))
