#load libraries
library(vegan) 
library(dplyr)
library(agricolae)
library(ggplot2)
library(reshape)

# Read metadata/mapping file
metadata <- read.table("Data/SBAR_Incubation_MappingFile_05142022.txt", 
                       sep = "\t", header = T, row.names = 1) %>%
            filter(Sample_type == "Box")
dim(metadata) #40   14

## Import realtive abundance file ####
box.rel <- read.table("Data/SBAR_16S_Abundance_Box_09262022.txt", 
                         sep = "\t", header = T, row.names = 1)
dim(box.rel) #9250 47

## Match asv with samples with metadata
all(colnames(box.rel[,1:40]) == rownames(metadata))


#### Aggregrate rel abundance at Phylum-Genus level#### 
# Combine with mapping file

taxa_list <- list()
md_taxa <- list()
for (i in 42:46){
  bac.taxa <- as.data.frame(box.rel[, 1:40]) #1:40 are the sample columns
  bac.taxa$Taxa <- box.rel[,i]
  bac.taxa <- aggregate(.~Taxa, 
                        data = bac.taxa, 
                        sum)
  rownames(bac.taxa) <- bac.taxa[,1]
  bac.taxa <- bac.taxa[,-1]
  taxa_list[[i - 41]] <- bac.taxa
  md_taxa[[i-41]] <- as.data.frame(cbind(metadata, t(taxa_list[[i-41]])))
}

#### Phylum Test for differences using KW ####
kw.phylum <- matrix(data = NA, nrow = 4, ncol = 41)
rownames(kw.phylum) = c('T0', 'T1','T2','T3')
colnames(kw.phylum) = colnames(md_taxa[[1]][, 15:55])

for(i in 15:55)
{
  bla <- kruskal(md_taxa[[1]][,i], md_taxa[[1]][,11], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.phylum[, i - 14] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.phylum), file = "Results/RelAbun_KW/phylum_KW_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### Class Test for differences using KW ####
kw.class <- matrix(data = NA, nrow = 4, ncol = 96)
rownames(kw.class) = c('T0', 'T1','T2','T3')
colnames(kw.class) = colnames(md_taxa[[2]][, 15:110])

for(i in 15:110)
{
  bla <- kruskal(md_taxa[[2]][,i], md_taxa[[2]][,11], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.class[, i - 14] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.class), file = "Results/RelAbun_KW/class_KW_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### Order Test for differences using KW ####
kw.order <- matrix(data = NA, nrow = 4, ncol = 195)
rownames(kw.order) = c('T0', 'T1','T2','T3')
colnames(kw.order) = colnames(md_taxa[[3]][, 15:209])

for(i in 15:209)
{
  bla <- kruskal(md_taxa[[3]][,i], md_taxa[[3]][,11], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.order[, i - 14] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.order), file = "Results/RelAbun_KW/order_KW_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### family Test for differences using KW ####
kw.family <- matrix(data = NA, nrow = 4, ncol = 254)
rownames(kw.family) = c('T0', 'T1','T2','T3')
colnames(kw.family) = colnames(md_taxa[[4]][, 15:268])

for(i in 15:268)
{
  bla <- kruskal(md_taxa[[4]][,i], md_taxa[[4]][,11], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.family[, i - 14] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.family), file = "Results/RelAbun_KW/family_KW_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### genus Test for differences using KW ####
kw.genus <- matrix(data = NA, nrow = 4, ncol = 486)
rownames(kw.genus) = c('T0', 'T1','T2','T3')
colnames(kw.genus) = colnames(md_taxa[[5]][, 15:500])

for(i in 15:500)
{
  bla <- kruskal(md_taxa[[5]][,i], md_taxa[[5]][,11], group = T, p.adj = "BH")
  bla.groups <- bla$groups
  bla.group.ordered <- bla.groups[order(as.character(row.names(bla.groups))),]
  kw.genus[, i - 14] <- as.character(bla.group.ordered[,2])
}

write.table(t(kw.genus), file = "Results/RelAbun_KW/genus_KW_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

write.table(md_taxa[[1]], file = "Results/RelAbundance/Bac_Phylum_Relabundance_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(md_taxa[[2]], file = "Results/RelAbundance/Bac_Class_Relabundance_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(md_taxa[[3]], file = "Results/RelAbundance/Bac_Order_Relabundance_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(md_taxa[[4]], file = "Results/RelAbundance/Bac_Family_Relabundance_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(md_taxa[[5]], file = "Results/RelAbundance/Bac_Genus_Relabundance_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)


#### Plot Phylum ####
phyla <- melt(md_taxa[[1]][,c(11,15:17,19,20,22,23,25,26,28,33,35,42,47,48,52, 53)])%>%
  group_by(Time_point, variable)%>%
  summarize(mean_rel_abund = 100* mean(value),
            sd_rel_abund = sd(value), .groups = "drop")

phyla$Phylum <- phyla$variable

phyla.plot <- ggplot(data = phyla, 
                 aes(y = Phylum,
                     x = mean_rel_abund, fill = Time_point))+
  geom_bar(position ="dodge", stat = "identity")+
  ylab("Phyla")+
  xlab("Mean Relative Abundance (%)")+
  geom_errorbar(aes(xmin = phyla$mean_rel_abund - phyla$sd_rel_abund, 
                    xmax = phyla$mean_rel_abund + phyla$sd_rel_abund), 
                width = 0.2, position = position_dodge(0.9))+
  theme_classic() 
  
  ggsave(plot = phyla.plot, 
         file = "Plots/sig_phylum_09262022.pdf", 
         width = 7, height = 6)

#### Plot Class ####
class <- melt(md_taxa[[2]][,c(11,16,18,20, 23,25,28,29, 35,38,41:43,
                               50,51,61,68,74,76,78,82,85:87, 97,103:106, 109,110)])%>%
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
  theme_classic() 
  
ggsave(plot = class.plot, 
         file = "Plots/sig_class_09262022.pdf", 
         width = 7, height = 8)
  
####Plot Order ####
order <- melt(md_taxa[[3]][,c(11,16,18,28,32,40,42,57,60,64,66,70,72,86,93,
                              94,97,105, 109, 110, 123,126,137, 147,152,
                              154, 155, 159, 163, 164, 166,172,182:184, 186,
                              196, 199, 201, 206)])%>%
  group_by(Time_point, variable)%>%
  summarize(mean_rel_abund = 100* mean(value),
            sd_rel_abund = sd(value), .groups = "drop")

order$Order <- order$variable

order.plot <- ggplot(data = order, 
                     aes(y = Order,
                         x = mean_rel_abund, fill = Time_point))+
  geom_bar(position ="dodge", stat = "identity")+
  ylab("Phyla")+
  xlab("Mean Relative Abundance (%)")+
  geom_errorbar(aes(xmin = order$mean_rel_abund - order$sd_rel_abund, 
                    xmax = order$mean_rel_abund + order$sd_rel_abund), 
                width = 0.2, position = position_dodge(0.9))+
theme_classic() 

ggsave(plot = order.plot, 
       file = "Plots/sig_order_09262022.pdf", 
       width = 8, height = 10)

#### Plot Family ####
fam <- melt(md_taxa[[4]][,c(11,19,23,35,38,40,44,45,50,58,61,68,90,92,
                            100,101,103,105,125,126,136, 147, 151,165,166,
                            168,169,183,190,192,198,203,208,211, 216, 217,227, 
                            228,229, 235,246,251,259,262,266)])%>%
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
                width = 0.2, position = position_dodge(0.9))+
theme_classic() 

ggsave(plot = fam.plot, 
       file = "Plots/sig_family_09262022.pdf", 
       width = 7, height = 10)

  