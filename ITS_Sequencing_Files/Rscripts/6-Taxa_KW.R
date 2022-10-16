#load libraries
library(vegan) 
library(dplyr)
library(agricolae)
library(ggplot2)
library(reshape)
library(ggbreak)

# Read metadata/mapping file
metadata <- read.table("Data/SBAR_Incubation_MappingFile_08232022.txt", 
                       sep = "\t", header = T, row.names = 1) %>%
            filter(Sample_type == "Box")
dim(metadata) #40   14
metadata <- metadata[-which(rownames(metadata)=="S_035"),]
metadata <- metadata[-which(rownames(metadata)=="S_045"),]

## Import relative abundance file ####
box.rel <- read.table("Data/SBAR_ITS_Abundance_Rar_Box_09232022.txt", 
                         sep = "\t", header = T, row.names = 1)
dim(box.rel) #1365 45

## Match asv with samples with metadata
all(colnames(box.rel[,1:38]) == rownames(metadata))

#### Aggregate rel abundance at Phylum-Genus level#### 
# Combine with mapping file
taxa_list <- list()
md_taxa <- list()
for (i in 40:44){
  taxa <- as.data.frame(box.rel[, 1:38]) #1:39 are the sample columns
  taxa$Taxa <- box.rel[,i]
  taxa <- aggregate(.~Taxa, 
                        data = taxa, 
                        sum)
  rownames(taxa) <- taxa[,1]
  taxa <- taxa[,-1]
  taxa_list[[i - 39]] <- taxa
  md_taxa[[i-39]] <- as.data.frame(cbind(metadata, t(taxa_list[[i-39]])))
}

#### Phylum Test for differences using KW ####
kw.phylum <- matrix(data = NA, nrow = 4, ncol = 11)
rownames(kw.phylum) = c('T0', 'T1','T2','T3')
colnames(kw.phylum) = colnames(md_taxa[[1]][, 14:24])

for(i in 14:24)
{
  bla <- kruskal(md_taxa[[1]][,i], md_taxa[[1]][,10], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.phylum[, i - 13] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.phylum), file = "Results/RelAbun_KW/ITS_phylum_KW_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### Class Test for differences using KW ####
kw.class <- matrix(data = NA, nrow = 4, ncol = 22)
rownames(kw.class) = c('T0', 'T1','T2','T3')
colnames(kw.class) = colnames(md_taxa[[2]][, 14:35])

for(i in 14:35)
{
  bla <- kruskal(md_taxa[[2]][,i], md_taxa[[2]][,10], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.class[, i - 13] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.class), file = "Results/RelAbun_KW/ITS_class_KW_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### Order Test for differences using KW ####
kw.order <- matrix(data = NA, nrow = 4, ncol = 50)
rownames(kw.order) = c('T0', 'T1','T2','T3')
colnames(kw.order) = colnames(md_taxa[[3]][, 14:63])

for(i in 14:63)
{
  bla <- kruskal(md_taxa[[3]][,i], md_taxa[[3]][,10], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.order[, i - 13] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.order), file = "Results/RelAbun_KW/ITS_order_KW_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### family Test for differences using KW ####
kw.family <- matrix(data = NA, nrow = 4, ncol = 95)
rownames(kw.family) = c('T0', 'T1','T2','T3')
colnames(kw.family) = colnames(md_taxa[[4]][, 14:108])

for(i in 14:108)
{
  bla <- kruskal(md_taxa[[4]][,i], md_taxa[[4]][,10], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.family[, i - 13] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.family), file = "Results/RelAbun_KW/ITS_family_KW_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### genus Test for differences using KW ####
kw.genus <- matrix(data = NA, nrow = 4, ncol = 165)
rownames(kw.genus) = c('T0', 'T1','T2','T3')
colnames(kw.genus) = colnames(md_taxa[[5]][, 14:178])

for(i in 14:178)
{
  bla <- kruskal(md_taxa[[5]][,i], md_taxa[[5]][,10], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.genus[, i - 13] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.genus), file = "Results/RelAbun_KW/genus_KW_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)


write.table(md_taxa[[1]], file = "Results/RelAbundance/Fungi_Phylum_Relabundance_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(md_taxa[[2]], file = "Results/RelAbundance/Fungi_Class_Relabundance_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(md_taxa[[3]], file = "Results/RelAbundance/Fungi_Order_Relabundance_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(md_taxa[[4]], file = "Results/RelAbundance/Fungi_Family_Relabundance_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(md_taxa[[5]], file = "Results/RelAbundance/Fungi_Genus_Relabundance_09232022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)


#### Plot Phylum ####
phyla <- melt(md_taxa[[1]][,c(10,14,15,17,19,21)])%>%
  group_by(Time_point, variable)%>%
  summarize(mean_rel_abund = 100* mean(value),
            sd_rel_abund = sd(value), .groups = "drop")

phyla$Phylum <- phyla$variable

phyla$Time_point <- factor(phyla$Time_point, 
                           levels = c("T0" , "T1", "T2", "T3"))

phyla.plot <- ggplot(data = phyla, 
                 aes(y = Phylum,
                     x = mean_rel_abund, fill = Time_point))+
  geom_bar(position ="dodge", stat = "identity")+
  ylab("Phyla")+
  xlab("Mean Relative Abundance (%)")+
  geom_errorbar(aes(xmin = phyla$mean_rel_abund - phyla$sd_rel_abund, 
                    xmax = phyla$mean_rel_abund + phyla$sd_rel_abund), 
                width = 0.2, position = position_dodge(0.9))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
  
ggsave(plot = phyla.plot, 
         file = "Plots/fungi_sig_phylum_09232022.pdf", 
         width = 7, height = 6)

#### Plot Class ####
class <- melt(md_taxa[[2]][,c(10,17:19, 21:24,28, 32,34 )])%>%
    group_by(Time_point, variable)%>%
    summarize(mean_rel_abund = 100* mean(value),
              sd_rel_abund = sd(value), .groups = "drop")
  
class$Class <- class$variable

class.plot <- ggplot(data = class, 
                       aes(y = Class,
                           x = mean_rel_abund, fill = Time_point))+
    geom_bar(position ="dodge", stat = "identity")+
    ylab("Phyla")+
    xlab("Mean Relative Abundance (%)")+
    geom_errorbar(aes(xmin = class$mean_rel_abund - class$sd_rel_abund, 
                      xmax = class$mean_rel_abund + class$sd_rel_abund), 
                  width = 0.2, position = position_dodge(0.9))+
  theme_bw() +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1))
  
ggsave(plot = class.plot, 
         file = "Plots/fungi_sig_class_09232022.pdf", 
         width = 7, height = 6)
  
#### Plot Family ####
fam <- melt(md_taxa[[4]][,c(10,19, 23,27,29,34,
                            43,47,48,53,55,58,63,65,68,
                            77,81,87,89,96,97,102)])%>%
  group_by(Time_point, variable)%>%
  summarize(mean_rel_abund = 100* mean(value),
            sd_rel_abund = sd(value), .groups = "drop")

fam$Family <- fam$variable

fam.plot <- ggplot(data = fam, 
                     aes(y = Family,
                         x = mean_rel_abund, fill = Time_point))+
  geom_bar(position ="dodge", stat = "identity")+
  ylab("Phyla")+
  xlab("Mean Relative Abundance (%)")+
  geom_errorbar(aes(xmin = fam$mean_rel_abund - fam$sd_rel_abund, 
                    xmax = fam$mean_rel_abund + fam$sd_rel_abund), 
                width = 0.2, position = position_dodge(0.9)) +
# scale_x_break(c(7.5, 45), scales = "free") +
scale_x_break(c(7, 47), scale = 0.3, ticklabels = c(47, 48)) +
scale_x_break(c(48, 63),ticklabels = c(63,64,65) ) +
theme_classic() 
#   theme(axis.text.x = element_text(angle = 45, 
#                                    hjust = 1))

ggsave(plot = fam.plot, 
       file = "Plots/fungi_sig_family_09232022.pdf", 
       width = 7, height = 6)

  