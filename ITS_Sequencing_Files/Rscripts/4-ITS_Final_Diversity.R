#Install packages
install.packages("vegan")
install.packages("ggplot2")
install.packages("devtools")
install.packages("reshape")
install.packages(dplyr) #for summary
install.packages(ggpubr)
install.packages(gridExtra)
install.packages("ggbiplot")
install.packages("ecodist")
install.packages("ape")
BiocManager::install("metagenomeSeq")

library(devtools)
install.packages("devtools")
install_github("vqv/ggbiplot")


#load libraries
library(vegan) #for multivariate statistics
library(ggplot2) #for visualization and plotting
library(dplyr) #for summary
library(agricolae)
library(metagenomeSeq)

# Import asv table without contaminants
bl.asv.tax <- read.table("Data/SBAR_ITS_ASV_table_Final_woCon_09222022.txt", 
                         sep = "\t", header = T, row.names = 1)
dim(bl.asv.tax) #2227  81

# Extract only the sequencing information into asv 
bl.asv <- t(bl.asv.tax[,-c(75:81)])
dim(bl.asv) #74 2227

# Extract taxa info 
taxa <- t(bl.asv.tax[,c(75:81)])
dim(taxa) #7 2227

# Read metadata/mapping file
metadata <- read.table("Data/SBAR_Incubation_MappingFile_08232022.txt", 
                       sep = "\t", header = T, row.names = 1, check.names = F)
dim(metadata) #89   13

#Remove samples not related to the project by matching samples IDs from mapping file
metadata <- metadata[rownames(bl.asv),] 
dim(metadata) #74 13

# Remove blanks and zymo from mapping file
metadata <- metadata %>%
  filter(Month != "blank") %>%  
  filter(Month != "Seq_Blank") %>% 
  filter(Month != "Zymo_Control") %>%
  filter(Sample_type == "Box")

dim(metadata) #39 13

## Match asv with samples wo blank
bl.asv <- bl.asv[rownames(metadata),] 
dim(bl.asv) #39 2227

bl.asv <- subset(bl.asv, select = colSums(bl.asv)!= 0)
dim(bl.asv) #39 1444

# check number of sequences per sample
hist(rowSums(bl.asv))
sort(rowSums(bl.asv))
summary(rowSums(bl.asv))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#6601   52866   58577   58300   67554   89803

#### Rarefraction ####
rar.asv <- rrarefy(bl.asv, 25000)
rar.asv <- as.data.frame(rar.asv)

# Remove unrarefied samples
rar.asv <- rar.asv[rowSums(rar.asv) == 25000,] #38 1444

rar.asv <- subset(rar.asv, select = colSums(rar.asv)!= 0)
dim(rar.asv) #38 1365

#Completely remove the three samples below rarefraction number
md <- metadata[-which(rownames(metadata)=="S_045"),]
all(rownames(md) == rownames(rar.asv))

## Match taxa with asv
tax <- taxa[,colnames(rar.asv)]
dim(tax) # 7 1365

sbar.rar <- cbind(t(rar.asv), t(tax))
dim(sbar.rar) #1365   45

write.table(sbar.rar,
            file = "Data/SBAR_ITS_Box_Rarefied_09232022.txt",
            sep = "\t",
            quote = F,
            row.names = T,
            col.names = NA)

# Import rarefied file
rar.asv.tax <- read.table("Data/SBAR_ITS_Box_Rarefied_09232022.txt", 
                         sep = "\t", header = T, row.names = 1)
dim(rar.asv.tax) #1365 45

# Extract only the sequencing information into asv 
rar.asv <- t(rar.asv.tax[,-c(39:45)])
dim(rar.asv) #38 1365

# Extract taxa info 
taxa <- t(rar.asv.tax[,c(39:45)])
dim(taxa) #7 1365

## import mapping file
md <- read.table("Data/SBAR_Incubation_MappingFile_08232022.txt",
                 sep = "\t", header = T, row.names = 1) %>%
  filter(Sample_type == "Box") 
dim(md) #40 13

md <- md[-which(rownames(md)=="S_035"),]
md <- md[-which(rownames(md)=="S_045"),]

all(rownames(rar.asv)== rownames(md))

### Calculate Diversity
#richness and shannon
md$Fungi_Richness <- specnumber(rar.asv)
md$Fungi_Shannon <- diversity(rar.asv, index="shannon")

mean(md$Fungi_Richness) #186.5789
mean(md$Fungi_Shannon)  #2.560492

ggplot(md, aes(Time_point, Fungi_Richness))+
  geom_jitter(aes(color = Time_point), 
              position = position_jitter(height = 0.05, width = 0.05), 
              alpha = 0.8, size = 3) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  theme_bw() + 
  theme(axis.text.x =  element_text(size = 8, hjust = 0.5, vjust = 1)) +
  ggtitle("Fungi Richness")

ggplot(md, aes(Time_point, Fungi_Shannon))+
  geom_jitter(aes(color = Time_point), 
              position = position_jitter(height = 0.05, width = 0.05), 
              alpha = 0.8, size = 3) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 0.5)) +
  ggtitle("Fungi Shannon")

###Kruskal-Wallis Richness + Shannon
bla <- kruskal(md$Fungi_Richness, md$Time_point, group = T, p.adj = "BH")
bla.groups <- bla$groups
#       metadata$Fungi_Richness groups
# T1          25.85000      a
# T0          25.33333      a
# T2          14.00000      b
# T3          12.72222      b

bla1 <- kruskal(md$Fungi_Shannon, md$Time_point, group = T, p.adj = "BH")
bla.groups1 <- bla1$groups
#       metadata$Fungi_Shannon groups
# T0         34.00000      a
# T1         17.70000      b
# T2         13.90000      b
# T3         13.22222      b

#NMDS ordination using Bray-Curtis distance
rownames(rar.asv) <- md$Sample
box.dist <- vegdist(rar.asv, method = "bray")
box.nmds <- metaMDS(box.dist, k = 2)
box.nmds$stress #0.08537972

md$Fungi_NMDS1 <- box.nmds$points[,1]
md$Fungi_NMDS2 <- box.nmds$points[,2]

centroid <- md %>%
  group_by(Time_point)%>%
  summarize(Fungi_NMDS1= mean(Fungi_NMDS1),
            Fungi_NMDS2= mean(Fungi_NMDS2))

segs <- merge(md, setNames(centroid, c('Time_point','oNMDS1','oNMDS2')),
              by = 'Time_point', sort = FALSE)

#Save the NMDS plots
fungi.nmds <- ggplot(md, aes(Fungi_NMDS1, Fungi_NMDS2,
                             label= Sample))+
  geom_point(aes(color = Time_point), 
             alpha = 0.8, size = 3) +
  geom_text(vjust = 0, nudge_y = 0.01, size=3) +
  # geom_segment(data = segs,
  #              mapping = aes(xend = oNMDS1, yend = oNMDS2, color =Time_point)) +
  # stat_ellipse(aes(group = Time_point, color = Time_point)) +
  # geom_point(data = centroid, size = 4, shape= 21, color= "black",
  #            aes(fill=Time_point)) +
  scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                     breaks = c("T0","T1", "T2", "T3")) +
  scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                    breaks = c("T0","T1", "T2", "T3")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate( "text", x = 0.01,
            y = -0.25,
            label = "Stress = 0.09") +
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "bottom")

ggsave(file = "Plots/Fungi_NMDS_SampleName_10032022.pdf", 
       plot = fungi.nmds,
       width = 150,
       height = 150,
       units ="mm")

adonis2(box.dist ~ md$Time_point)
#               Df SumOfSqs      R2     F Pr(>F)    
# md$Time_point  3  0.39061 0.29765 4.8029  0.001 ***
# Residual      34  0.92171 0.70235                  
# Total         37  1.31233 1.00000 

#Calculate Beta-dispersivity
box.dispersion <- betadisper(box.dist, group = md$Time_point,type = "centroid")

#Check Average distance to centroid values
box.dispersion

md <- cbind(md, box.dispersion$distances)
colnames(md)[14] <- "Distance_to_centroid"

#Test whether the beta-dispersivity across Time points is different
permutest(box.dispersion)
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     3 0.03009 0.0100299 5.8404    999  0.002 **
# Residuals 34 0.05839 0.0017173

## Kruskal-wallis
bla3 <- kruskal(md$Distance_to_centroid, md$Time_point, group = T, p.adj = "BH")
bla.groups3 <- bla3$groups
#       md$Distance_to_centroid groups
# T0                31.00000      a
# T1                20.40000      b
# T2                14.90000      b
# T3                12.11111      b

#Plot the distance to centroid values
pdf("Plots/Fungi_Beta_dispersivity_10022022.pdf", 
    width = 5, height = 4)
boxplot(box.dispersion, xlab ="Time Point", ylab = "Distance to centroid")
dev.off()

### Plot Bray-curtis distance comparisons ###
dist_com <- box.dist %>%
  as.matrix() %>%
  as_tibble(rownames = "samples") %>%
  pivot_longer(-samples) %>%
  filter(samples < name) %>%
  separate(samples, into = c("time_a", "box_a"), "B") %>%
  separate(name, into = c("time_b", "box_b"), "B") %>%
  mutate(comparison = case_when(
    time_a != time_b & time_a == "T0" & time_b == "T1" ~ "first",
    time_a != time_b & time_a == "T0" & time_b == "T2" ~ "second",
    time_a != time_b & time_a == "T0" & time_b == "T3" ~ "third",
    TRUE ~ NA_character_
  )) %>% 
  drop_na()

dist.comp.plot <- dist_com %>% 
  ggplot(aes(x = comparison, y = value)) +
  geom_jitter (width = 0.1, color ="#4C4E52") +
  stat_summary(fun.data = median_hilow , color ="red",
               size =1, fun.args = list(conf.int = 0.50)) +
  labs( x= NULL, y = "Bray-Curtis distances \n for fungal community") +
  scale_x_discrete(breaks =c("first", "second", "third"),
                   labels =c("T0 vs. T1",
                             "T0 vs. T2",
                             "T0 vs. T3")) +
  scale_y_continuous(limits = c(0.1,0.5), breaks= seq(0.1, 0.5, 0.1)) +
  theme_classic()


bla2 <- kruskal(dist_com$value, dist_com$comparison, group = T, p.adj = "BH")
bla.groups2 <- bla2$groups
#            dist_com$value groups
# second       138.6667      a
# third        136.5000      a
# first        118.3833      a

ggsave(file = "Plots/Bac_BrayCurtis_Dist_Com_10032022.pdf", 
       plot = dist.comp.plot,
       width = 100,
       height = 100,
       units ="mm")

write.table(md,
            file = "Results/SBAR_ITS_Rar_Diversity_09232022.txt",
            sep = "\t",
            quote = F,
            row.names = T,
            col.names = NA)


#### NMDS T0 vs T1 ####
# Comparison 1: T0 vs T1
md1 <- md %>%
          filter(Time_point == "T0" |Time_point == "T1")

asv1 <- rar.asv[rownames(md1),]
asv1 <- subset(asv1, select = colSums(asv1)!= 0)
dim(asv1) #19 971

dist1 <- vegdist(asv1, method = "bray")
nmds1 <- metaMDS(dist1, k = 2)
nmds1$stress #0.07014334

md1$Fungi_NMDS1 <- nmds1$points[,1]
md1$Fungi_NMDS2 <- nmds1$points[,2]

pdf("Plots/Fungi_NMDS_T0_T1_10022022.pdf", 
    width = 5, height = 4)
ggplot(md1, aes(Fungi_NMDS1, Fungi_NMDS2))+
  geom_point(aes(color = Time_point, shape= Box), 
             alpha = 0.8, size = 4) +
  theme_bw() +
  theme(legend.position="right") +
  annotate( "text", x = 0.1,
            y = -0.25,
            label = "Stress = 0.07") +
  annotate( "text", x = -0.18,
            y = 0.2,
            label = bquote(R^2~'(Time)= 0.26***'))
dev.off()

adonis2(dist1 ~ md1$Time_point)
#               Df SumOfSqs      R2     F Pr(>F)    
# md1$Time_point  1  0.18899 0.24527 5.5245  0.001 ***
# Residual       17  0.58156 0.75473                  
# Total          18  0.77055 1.00000 

#Calculate Beta-dispersivity
dispersion1 <- betadisper(dist1, group = md1$Time_point,type = "centroid")

#Check Average distance to centroid values
dispersion1
#Test whether the beta-dispersivity across Microsite is different
permutest(dispersion1)
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.011051 0.0110511 7.0617    999  0.016 *
# Residuals 17 0.026604 0.0015649   

#Plot the distance to centroid values
boxplot(dispersion1, xlab ="Time Point", ylab = "Distance to centroid")


#### NMDS T0 vs T2 ####
# Comparison 2: T0 vs T2
md2 <- md %>%
  filter(Time_point == "T0" |Time_point == "T2")

asv2 <- rar.asv[rownames(md2),]
asv2 <- subset(asv2, select = colSums(asv2)!= 0)
dim(asv2) #19 941

dist2 <- vegdist(asv2, method = "bray")
nmds2 <- metaMDS(dist2, k = 2)
nmds2$stress #0.05797367


md2$Fungi_NMDS1 <- nmds2$points[,1]
md2$Fungi_NMDS2 <- nmds2$points[,2]

pdf("Plots/Fungi_NMDS_T0_T2_10022022.pdf", 
    width = 5, height = 4)
ggplot(md2, aes(Fungi_NMDS1, Fungi_NMDS2))+
  geom_point(aes(color = Time_point, shape= Box), 
             alpha = 0.8, size = 4) +
  theme_bw() +
  theme(legend.position="right") +
  annotate( "text", x = 0.1,
            y = -0.25,
            label = "Stress = 0.06") +
  annotate( "text", x = -0.18,
            y = 0.2,
            label = bquote(R^2~'(Time)= 0.30***'))
dev.off()

adonis2(dist2 ~ md2$Time_point)
#               Df SumOfSqs      R2     F Pr(>F)    
# md2$Time_point  1  0.24785 0.30354 7.4091  0.001 ***
# Residual       17  0.56869 0.69646                  
# Total          18  0.81654 1.00000                  

#Calculate Beta-dispersivity
dispersion2 <- betadisper(dist2, group = md2$Time_point,type = "centroid")

#Check Average distance to centroid values
dispersion2
#Test whether the beta-dispersivity across Microsite is different
permutest(dispersion2)
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.015765 0.0157655 6.646    999  0.015 *
#   Residuals 17 0.040327 0.0023722 

#Plot the distance to centroid values
boxplot(dispersion2, xlab ="Time Point", ylab = "Distance to centroid")

#### NMDS T0 vs T3 ####
# Comparison 3: T0 vs T3
md3 <- md %>%
  filter(Time_point == "T0" |Time_point == "T3")

asv3 <- rar.asv[rownames(md3),]
asv3 <- subset(asv3, select = colSums(asv3)!= 0)
dim(asv3) #18 890

dist3 <- vegdist(asv3, method = "bray")
nmds3 <- metaMDS(dist3, k = 2)
nmds3$stress #0.04680908

md3$Fungi_NMDS1 <- nmds3$points[,1]
md3$Fungi_NMDS2 <- nmds3$points[,2]

pdf("Plots/Fungi_NMDS_T0_T3_10022022.pdf", 
    width = 5, height = 4)
ggplot(md3, aes(Fungi_NMDS1, Fungi_NMDS2))+
  geom_point(aes(color = Time_point, shape= Box), 
             alpha = 0.8, size = 4) +
  theme_bw() +
  theme(legend.position="right") +
  annotate( "text", x = 0.08,
            y = -0.25,
            label = "Stress = 0.05") +
  annotate( "text", x = -0.18,
            y = 0.2,
            label = bquote(R^2~'(Time)= 0.35***'))
dev.off()

adonis2(dist3 ~ md3$Time_point)
#               Df SumOfSqs      R2     F Pr(>F)    
# md3$Time_point  1  0.26064 0.35445 8.7851  0.001 ***
# Residual       16  0.47469 0.64555                  
# Total          17  0.73533 1.00000 

#Calculate Beta-dispersivity
dispersion3 <- betadisper(dist3, group = md3$Time_point,type = "centroid")

#Check Average distance to centroid values
dispersion3
#Test whether the beta-dispersivity across Microsite is different
permutest(dispersion3)
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.027888 0.0278877 34.511    999  0.001 ***
# Residuals 16 0.012929 0.0008081

#Plot the distance to centroid values
boxplot(dispersion3, xlab ="Time Point", ylab = "Distance to centroid")

#### Normalization ####
fungi.MR <- newMRexperiment(t(bl.asv))
p <- cumNormStat(fungi.MR)
fungi.MR <- cumNorm(fungi.MR, p = p)
fungi.norm <- t(MRcounts(fungi.MR, norm = T, log = F))
dim(fungi.norm) #39 1444

# check number of sequences per sample after normalization
hist(rowSums(fungi.norm))
sort(rowSums(fungi.norm))
summary(rowSums(fungi.norm))

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17697   37093   47852   44914   52513   67930

all(rownames(metadata)== rownames(fungi.norm))

## Match taxa with asv
taxa <- taxa[,colnames(fungi.norm)]
dim(taxa) # 7 1444

sbar.rar <- cbind(t(fungi.norm), t(taxa))
dim(sbar.rar) #1444   46

write.table(sbar.rar,
            file = "Data/SBAR_ITS_Box_Normalized_09232022.txt",
            sep = "\t",
            quote = F,
            row.names = T,
            col.names = NA)

### Calculate Diversity
#richness and shannon
metadata$Fungi_Richness <- specnumber(fungi.norm)
metadata$Fungi_Shannon <- diversity(fungi.norm, index="shannon")

mean(metadata$Fungi_Richness) #187.1282
mean(metadata$Fungi_Shannon)  #2.574082

ggplot(metadata, aes(Time_point, Fungi_Richness))+
  geom_jitter(aes(color = Time_point), 
              position = position_jitter(height = 0.05, width = 0.05), 
              alpha = 0.8, size = 3) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  theme_bw() + 
  theme(axis.text.x =  element_text(size = 8, hjust = 0.5, vjust = 1)) +
  ggtitle("Fungi Richness")

ggplot(metadata, aes(Time_point, Fungi_Shannon))+
  geom_jitter(aes(color = Time_point), 
              position = position_jitter(height = 0.05, width = 0.05), 
              alpha = 0.8, size = 3) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 0.5)) +
  ggtitle("Fungi Shannon")

###Kruskal-Wallis Richness + Shannon
bla <- kruskal(metadata$Fungi_Richness, metadata$Time_point, group = T, p.adj = "BH")
bla.groups <- bla$groups
#       metadata$Fungi_Richness groups
# T1                26.80000      a
# T0                22.60000      a
# T2                15.90000      a
# T3                14.11111      a

bla1 <- kruskal(metadata$Fungi_Shannon, metadata$Time_point, group = T, p.adj = "BH")
bla.groups1 <- bla1$groups
#       metadata$Fungi_Shannon groups
# T0               34.50000      a
# T1               17.90000      b
# T2               13.60000      b
# T3               13.33333      b

#NMDS ordination using Bray-Curtis distance
box.dist <- vegdist(fungi.norm, method = "bray")
box.nmds <- metaMDS(box.dist, k = 2)
box.nmds$stress #0.07116345

metadata$Fungi_NMDS1 <- box.nmds$points[,1]
metadata$Fungi_NMDS2 <- box.nmds$points[,2]

#Save the NMDS plots
ggplot(metadata, aes(Fungi_NMDS1, Fungi_NMDS2))+
  geom_point(aes(color = Time_point, shape = Box), 
             alpha = 0.8, size = 4) +
  theme_bw() +
  theme(legend.position="right") +
  annotate( "text", x = 0.35,
            y = -0.25,
            label = "Stress = 0.09") +
  annotate( "text", x = -0.18,
            y = 0.2,
            label = bquote(R^2~'(Time)= 0.30***'))

adonis2(box.dist ~ metadata$Time_point)
#               Df SumOfSqs      R2     F Pr(>F)    
# metadata$Time_point  3  0.85242 0.36111 6.5942  0.001 ***
# Residual            35  1.50813 0.63889                  
# Total               38  2.36055 1.00000

write.table(metadata,
            file = "Results/SBAR_ITS_Norm_Diversity_09232022.txt",
            sep = "\t",
            quote = F,
            row.names = T,
            col.names = NA)
