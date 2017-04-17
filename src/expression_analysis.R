
# R package dependencies -------------------------------------------------------
# Comment this line out if your machine already has these packages installed
#install.packages(c("matrixStats","ggplot2","MASS","scales"),repos='https://watson.nci.nih.gov/cran_mirror/')

## FPKM matrix input -----------------------------------------------------------
data1 <- read.table("./RESULTS/normout/genes.count_table",header=TRUE,sep="\t", stringsAsFactors = FALSE)
attr.table <- read.table("./RESULTS/normout/genes.attr_table",header=TRUE,sep="\t",stringsAsFactors = FALSE)
data1$gene_short_name <- attr.table$gene_short_name
# write out complete gene expression matrix
write.csv(data1, file="./RESULTS/HIO_patterning_FPKM.csv")

png(filename="./RESULTS/PCA.png")
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
## Clear temporary memory
rm(list=ls())

library(ggplot2)
library(RColorBrewer)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Generate FPKM data table
data <- read.table("./RESULTS/normout/genes.count_table",header=TRUE,sep="\t")
data1 <-data
attr.table <- read.table("./RESULTS/normout/genes.attr_table",header=TRUE,sep="\t")
data1$gene_short_name <- attr.table$gene_short_name
write.csv(data1, file="./RESULTS/cufflinks_normalizedFPKM_table.csv")
# PCA analysis
num.data <- data[,sapply(data,is.numeric)]
group <- gsub('.{2}$', '', colnames(num.data))
pca.data <- num.data[apply(num.data, 1, sd, na.rm=TRUE) != 0,]
pca <- prcomp(t(pca.data),scale=TRUE,center=TRUE,na.action=na.omit)
scores <- data.frame(colnames(pca.data), pca$x[,1:ncol(pca$x)],group)

library(ggplot2)
theme <- theme(legend.position="right",
               legend.title=element_blank(),
               legend.background = element_rect(fill="white", size=0.5, linetype="dotted"))
# plot theme
library(ggplot2)
library(grid)
theme1 <-  theme(axis.text.x=element_text(color="black",
                                          size=16,
                                          angle=0,
                                          hjust=0.5,
                                          face="bold"),
                 axis.text.y=element_text(size=16,
                                          face="bold",
                                          hjust=1,
                                          color="black"),
                 legend.position="none",
                 panel.background = element_rect(fill = "white"),
                 panel.grid.major=element_line(size=0.5,
                                               color = "grey40",
                                               linetype = "dashed"),
                 panel.grid.minor=element_line(size=0.5, # element_blank()
                                               color = "grey70",
                                               linetype = "dashed"),
                 axis.title=element_text(size=16,
                                         face="bold"),
                 axis.title.y=element_text(size=16,
                                           face="bold", vjust = 1.5),
                 axis.title.x=element_text(vjust=-0.5),
                 legend.title=element_blank(),
                 panel.border=element_rect(fill=NA,
                                           color="black",size=2),
                 plot.title = element_text(size=16,
                                           face="bold",
                                           hjust = -0.1),
                 legend.text=element_text(size=14,
                                          face="bold"),
                 axis.title.y=element_text(vjust=3,
                                           size=16),
                 plot.margin = unit(c(1,1,1,1), "cm"))
pc1.2 <- qplot(x=PC1, y=PC2, data=scores) + theme1 + theme +
 # scale_fill_brewer(palette="Set3")+
  scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(scores$group)), "Set1"))(length(unique(scores$group)))) +
  geom_point(shape=21,aes(fill=factor(group)), size=8)  

pc1.3 <- qplot(x=PC1, y=PC3, data=scores) + theme1 + theme +
 # scale_fill_brewer(palette="Set3")+
  scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(scores$group)), "Set1"))(length(unique(scores$group)))) +
  geom_point(shape=21,aes(fill=factor(group)), size=8) 

pc2.3 <- qplot(x=PC2, y=PC3, data=scores) + theme1 + theme +
  #scale_fill_brewer(palette="Set3")+
  scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(scores$group)), "Set1"))(length(unique(scores$group)))) +
  geom_point(shape=21,aes(fill=factor(group)), size=8) 
png(file="./RESULTS/PCA.png",width=11,height=8, units="in", res =144)
multiplot(pc1.2,pc1.3,pc2.3,cols=2)
dev.off()
# Spearman Correlation
cor1 <- cor(pca.data,method="spearman")
## Load libraries
library(matrixStats)
library(gplots)
library(RColorBrewer)

# creates a color palette from blue to white to red in "n" increments
my_palette <- colorRampPalette(c("white","white","red"))(n = 299)

# defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,0.8,length=100),  #white
               seq(0.8,0.97,length=100),
               seq(0.95,1,length=100))   #red    

cor.mat <- as.matrix(cor1)
png(file="./RESULTS/spearman_correlation.png",width=8,height=8, units = "in", res = 144)
library(pheatmap)
pheatmap(cor1,clustering_method="average",color=colorRampPalette(rev(brewer.pal(n = 7, name =
"RdYlBu")))(300),main="Spearman Correlation",border_color="black",fontsize=20)
dev.off()
dev.off()
