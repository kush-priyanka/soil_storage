## Load libraries
library(ggplot2)
library(reshape)
library(dplyr)
library(agricolae)

## Import Rarefied file: Bac/Arch rarefied file
bl.asv.tax <- read.table("Data/SBAR_16S_Box_Rarefied_09262022.txt", 
                         sep = "\t", header = T, row.names = 1)
dim(bl.asv.tax)

## Subset asv read and taxonomy info into two variables
# Extract only the sequencing information into asv
bl.asv <- bl.asv.tax[,-c(41:47)]
dim(bl.asv) #9250 40

# Extract taxa info
taxa <- bl.asv.tax[,c(41:47)]
dim(taxa) #9250 7


all(rownames(bl.asv) == rownames(taxa))
bl.asv$taxonomy<-paste(taxa$Kingdom, 
                        taxa$Phylum, 
                        taxa$Class,
                        taxa$Order, 
                        taxa$Family, 
                        taxa$Genus,
                        taxa$Species, sep=";")

## Save the file to use as FAPROTAX input file
write.table(bl.asv, file = "Data/16S_rar_asv_faprotax_09272022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)


## Import faprotax output file with functional groups and read counts
faprotax <- t(read.table("Results/FAPROTAX/sbar_faprotax_09272022.txt",
                         sep = "\t", 
                         header = T, row.names = 1))
dim(faprotax) #40 92
faprotax.filter <- subset(faprotax,
                          select = colSums(faprotax)!= 0)
dim(faprotax.filter) #40 50

# bac.asv.t <- bac.asv[,-20] #remove taxonomy column
# bac.asv.t <- t(bac.asv.t)

## Calculate proportions for each function
faprotax.filter <- faprotax.filter/40000  
faprotax.filter <- subset(faprotax.filter, 
                          select = colSums(faprotax.filter)!= 0) #40 50

## Save a file with function proportion
write.table(faprotax.filter, file = "Results/FAPROTAX/SBAR_16S_faprotax_results_propor_09272022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

# ##Import file with function proportion
# faprotax.filter <- read.table(file = "Results/FAPROTAX/SBAR_16S_faprotax_results_propor_09272022.txt", 
#                               sep = "\t", header = T, row.names = 1)

## import mapping file
md <- read.table("Data/SBAR_Incubation_MappingFile_05142022.txt",
                sep = "\t", header = T, row.names = 1) %>%
      filter(Sample_type == "Box") 
dim(md) #40 14

faprotax.filter <- faprotax.filter[rownames(md),]
all(rownames(faprotax.filter) == rownames(md))

faprotax.filter <- subset(faprotax.filter,
                          select = colSums(faprotax.filter)!= 0)
dim(faprotax.filter) #40 50

## Merge mapping file with function proportion
map.faprotax <- cbind(md, faprotax.filter)


## Save mapping file with function proportion
write.table(map.faprotax, 
            file="Results/FAPROTAX/SBAR_16S_Mapping_faprotax_results_propor_09272022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

## Subset the mapping file to plot selected functions
carbon <- melt(map.faprotax[,c(11,15:17, 27, 31:32, 38, 45:48, 58:62)])
com <- melt(map.faprotax[,c(11,18,19,39,64)])
nitro <- melt(map.faprotax[,c(11,23:26, 29:30,50:52, 63)])
sul <- melt(map.faprotax[,c(11,20:22, 33:36, 57)])
misc <- melt(map.faprotax[,c(11,37, 40:44, 49,53:56)])


c1 <- ggplot(data = carbon, 
       aes(x = variable, y = value, fill = Time_point))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  # geom_point(position = position_jitterdodge()) +
  xlab("Carbon related")+
  ylab("Bac/Arc Functional Potential Proportion")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

ggsave(plot = c1, 
       file = "Plots/Faprotax/faprtotax_carbon_09272022.pdf", 
       width = 12, height = 6)

n1 <- ggplot(data = nitro, 
             aes(x = variable, y = value, fill = Time_point))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  # geom_point(position = position_jitterdodge()) +
  xlab("Nitrogen related")+
  ylab("Bac/Arc Functional Potential Proportion")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

ggsave(plot = n1, 
       file = "Plots/Faprotax/faprtotax_nitro_09272022.pdf", 
       width = 10, height = 6)

s1 <- ggplot(data = sul, 
             aes(x = variable, y = value, fill = Time_point))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  # geom_point(position = position_jitterdodge()) +
  xlab("Sulfur related")+
  ylab("Bac/Arc Functional Potential Proportion")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

ggsave(plot = s1, 
       file = "Plots/Faprotax/faprtotax_sulfur_09272022.pdf", 
       width = 6, height = 4)

c2 <- ggplot(data = com, 
       aes(x = variable, y = value, fill = Time_point))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  # geom_point(position = position_jitterdodge()) +
  xlab("Highest Proportion Functions")+
  ylab("Bac/Arc Functional Potential \nProportion")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
ggsave(plot = c2, 
       file = "Plots/Faprotax/faprtotax_highest_prop_groups_09272022.pdf", 
       width = 6, height = 4)

m <- ggplot(data = misc, 
       aes(x = variable, y = value, fill = Time_point))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  # geom_point(position = position_jitterdodge()) +
  xlab("Misc Functions")+
  ylab("Bac/Arc Functional Potential Proportion")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

ggsave(plot = m, 
       file = "Plots/Faprotax/faprtotax_misc_09272022.pdf", 
       width = 6, height = 4)


#### Matrix for sig letters Soil ####
kw.matrix <- matrix(data = NA, nrow = 4, ncol = 50)
rownames(kw.matrix) = c('T0', 'T1','T2','T3')
colnames(kw.matrix) = colnames(map.faprotax[, 15:64])

for(i in 15:64)
{
  bla <- kruskal(map.faprotax[,i], map.faprotax$Time_point, group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.matrix[, i - 14] <- as.character(bla.group.ordered[,2])
  # print(bla.group.ordered)
}

write.table(t(kw.matrix), file = "Results/FAPROTAX/SBAR_faprotax_KW_09272022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)