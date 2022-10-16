## Load libraries
library(ggplot2)
library(reshape)
library(dplyr)
library(agricolae)

## Import Rarefied file: Bac/Arch rarefied file
asv.tax <- read.table("Data/SBAR_ITS_Box_Rarefied_09232022.txt", 
                         sep = "\t", header = T, row.names = 1)
dim(asv.tax) #1365 45

## Subset asv read and taxonomy info into two variables
# Extract only the sequencing information into asv
asv <- asv.tax[,-c(39:45)]
dim(asv) #1365 38

# Extract taxa info
taxa <- asv.tax[,c(39:45)]
dim(taxa) #1365 7


all(rownames(asv) == rownames(taxa))
asv$taxonomy<-paste(taxa$Kingdom, 
                       taxa$Phylum, 
                       taxa$Class,
                       taxa$Order, 
                       taxa$Family, 
                       taxa$Genus,
                       taxa$Species, sep=";")

## Save the file to use as FAPROTAX input file
write.table(asv, file = "Funguild/ITS_rar_asv_funguild_09272022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

## Import faprotax output file with functional groups and read counts
funguild <- t(read.table("Funguild/ITS_rar_asv_funguild_09272022.guilds_matched.txt",
                         sep = "\t", 
                         header = T, row.names = 1))
dim(funguild) #48 670
funguild <- as.data.frame(funguild)

fun.asv.t <- as.data.frame(funguild[1:38,])
fun.asv.t <- apply(fun.asv.t, 2, function (x) {as.numeric(as.character(x))})
rownames(fun.asv.t) <- rownames(funguild[1:38,])

guilds <- as.data.frame(apply(fun.asv.t, 1, function (x) by(x, as.factor(funguild[43,]), sum)))
guilds <- subset(guilds, 
                          select = colSums(guilds)!= 0)

## Calculate proportion for guilds
## Use the rarefraction value
guilds.proportion <- guilds/25000 

## Save a file with function proportion
write.table(guilds.proportion, file = "Results/SBAR_ITS_funguild_results_propor_09272022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)


#Check for columns with many zeroes
guilds.t <- t(guilds.proportion)
guilds.t <- subset(guilds.t, 
                   select=  colSums(guilds.t)!= 0)

## import mapping file
md <- read.table("Data/SBAR_Incubation_MappingFile_08232022.txt",
                 sep = "\t", header = T, row.names = 1) %>%
      filter(Sample_type == "Box") 
dim(md) #40 13
md <- md[-which(rownames(md)=="S_035"),]
md <- md[-which(rownames(md)=="S_045"),]

all(rownames(guilds.t) == rownames(md))

## Merge mapping file with function proportion
map.guild <- cbind(md, guilds.t)

## Save mapping file with function proportion
write.table(map.guild, 
            file="Results/SBAR_ITS_Mapping_funguild_results_propor_09272022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### Matrix for sig letters Soil ####
kw.matrix <- matrix(data = NA, nrow = 4, ncol = 52)
rownames(kw.matrix) = c('T0', 'T1','T2','T3')
colnames(kw.matrix) = colnames(map.guild[, 14:65])

for(i in 14:65)
{
  bla <- kruskal(map.guild[,i], map.guild$Time_point, group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.matrix[, i - 13] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.matrix), file = "Results/SBAR_funguild_KW_09272022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

###Plot funguild taxa ###
subset <- melt(map.guild[,c(10,16,25, 29,31,41,45,49,57,61,63, 65)])

guild.plot <- ggplot(data = subset, 
             aes(x = variable, y = value, fill = Time_point))+
  geom_boxplot(alpha = 0.7)+
  xlab("")+
  ylab("Fungal Functional \nPotential Proportion")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

ggsave(plot = guild.plot, 
       file = "Plots/Funguild_09272022.pdf", 
       width = 12, height = 7)
