library(vegan) 
library(dplyr)
library(agricolae)
library(ggplot2)

## Import asv table w/ Feb Field and Box samples####
bl.asv.tax <- read.table("Data/SBAR_16S_ASV_table_Final_05182022.txt", 
                         sep = "\t", header = T, row.names = 1)
dim(bl.asv.tax) #17254    96

# Extract only the sequencing information into asv
bl.asv <- t(bl.asv.tax[,-c(90:96)])
dim(bl.asv) #89 17254

# Extract taxa info 
taxa <- t(bl.asv.tax[,c(90:96)])
dim(taxa) #7 17254


# Read metadata/mapping file
metadata <- read.table("Data/SBAR_Incubation_MappingFile_05142022.txt", 
                       sep = "\t", header = T, row.names = 1, check.names = F) %>%
            filter(Sample_type == "Box")

dim(metadata) #40   14

## Match asv with samples wo blank
box.asv <- bl.asv[rownames(metadata),] 
dim(box.asv) #40 17254

## Match asv with samples with metadata
all(rownames(box.asv) == rownames(metadata))

box.asv <- subset(box.asv, select = colSums(box.asv)!= 0)

# check number of sequences per sample
hist(rowSums(box.asv))
sort(rowSums(box.asv))
summary(rowSums(box.asv))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#41296   75152  102407  118459  149396  360153

#### Rarefraction ####
rar.asv <- rrarefy(box.asv, 40000)
rar.asv <- as.data.frame(rar.asv) #40 10177

rar.asv <- subset(rar.asv, select = colSums(rar.asv)!= 0)
dim(rar.asv) #40 9250

all(rownames(metadata) == rownames(rar.asv))

## Match taxa with asv
tax <- taxa[,colnames(rar.asv)]
dim(tax) # 7 9250

sbar.rar <- cbind(t(rar.asv), t(tax))
dim(sbar.rar) #9250 47

write.table(sbar.rar,
            file = "Data/SBAR_16S_Box_Rarefied_09262022.txt",
            sep = "\t",
            quote = F,
            row.names = T,
            col.names = NA)

## Import rarefied table w/ Feb Field and Box samples####
bl.asv.tax <- read.table("Data/SBAR_16S_Box_Rarefied_09262022.txt", 
                         sep = "\t", header = T, row.names = 1)
dim(bl.asv.tax) #9250 47

# Extract only the sequencing information into asv
bl.asv <- t(bl.asv.tax[,-c(41:47)])
dim(bl.asv) #40 9250


# Read metadata/mapping file
metadata <- read.table("Data/SBAR_Incubation_MappingFile_05142022.txt", 
                       sep = "\t", header = T, row.names = 1, check.names = F) %>%
  filter(Sample_type == "Box")

dim(metadata) #40   14

## Match asv with samples wo blank
box.asv <- bl.asv[rownames(metadata),] 
dim(box.asv) #40 9250

## Match asv with samples with metadata
all(rownames(box.asv) == rownames(metadata))

### Calculate Diversity
#richness and shannon
metadata$Bac_Richness <- specnumber(box.asv)
metadata$Bac_Shannon <- diversity(box.asv, index="shannon")

mean(metadata$Bac_Richness) #1610.5
mean(metadata$Bac_Shannon) # 6.444276

ggplot(metadata, aes(Time_point, Bac_Richness))+
  geom_jitter(aes(color = Time_point, shape = Box), 
              position = position_jitter(height = 0.05, width = 0.05), 
              alpha = 0.8, size = 3) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  theme_bw() + 
  theme(axis.text.x =  element_text(size = 8, hjust = 0.5, vjust = 1)) +
  ggtitle("Bac/Arch Richness")

ggplot(metadata, aes(Time_point, Bac_Shannon))+
  geom_jitter(aes(color = Time_point), 
              position = position_jitter(height = 0.05, width = 0.05), 
              alpha = 0.8, size = 3) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 0.5)) +
  ggtitle("Bac/Arch Shannon")

###Kruskal-Wallis Richness + Shannon
bla <- kruskal(metadata$Bac_Richness, metadata$Time_point, group = T, p.adj = "BH")
bla.groups <- bla$groups
#       metadata$Bac_Richness groups
# T1                  27.9      a
# T0                  21.1      a
# T2                  18.3      a
# T3                  14.7      a

bla1 <- kruskal(metadata$Bac_Shannon, metadata$Time_point, group = T, p.adj = "BH")
bla.groups1 <- bla1$groups
#       metadata$Bac_Shannon groups
# T1                 27.3      a
# T0                 27.0      a
# T2                 15.8      b
# T3                 11.9      b

#NMDS ordination using Bray-Curtis distance
box.dist <- vegdist(box.asv, method = "bray")
box.nmds <- metaMDS(box.dist, k = 2)
box.nmds$stress #0.1353536

metadata$Bac_NMDS1 <- box.nmds$points[,1]
metadata$Bac_NMDS2 <- box.nmds$points[,2]

centroid <- metadata %>%
  group_by(Time_point)%>%
  summarize(Bac_NMDS1= mean(Bac_NMDS1),
            Bac_NMDS2= mean(Bac_NMDS2))

segs <- merge(metadata, setNames(centroid, c('Time_point','oNMDS1','oNMDS2')),
              by = 'Time_point', sort = FALSE)

#Save the NMDS plots
bac.nmds <- ggplot(metadata, aes(Bac_NMDS1, Bac_NMDS2)) +
  geom_point(aes(color = Time_point, shape= Box), size = 3) +
  # geom_text(vjust = 0, nudge_y = 0.01, size=3) +
  # stat_ellipse(aes(group = Time_point, color = Time_point)) +
  # geom_segment(data = segs,
  #              mapping = aes(xend = oNMDS1, yend = oNMDS2, color =Time_point)) +
  # geom_point(data = centroid, size = 4, shape= 21, color= "black",
  #            aes(fill=Time_point)) +
  scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                     breaks = c("T0","T1", "T2", "T3")) +
  scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                     breaks = c("T0","T1", "T2", "T3")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate( "text", x = 0.1,
            y = -0.3,
            label = "Stress = 0.14",
            size = 3) +
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "right")

ggsave(file = "Plots/Bac_NMDS_Box_10032022.pdf", 
       plot = bac.nmds,
       width = 150,
       height = 150,
       units ="mm")

adonis2(box.dist ~ metadata$Time_point)
#                     Df SumOfSqs      R2      F Pr(>F)    
# metadata$Time_point  3   0.3043 0.17511 2.5474  0.001 ***
# Residual            36   1.4335 0.82489                  
# Total               39   1.7378 1.00000

#Calculate Beta-dispersivity
box.dispersion <- betadisper(box.dist, group = metadata$Time_point,type = "centroid")

#Check Average distance to centroid values
box.dispersion

metadata <- cbind(metadata, box.dispersion$distances)
colnames(metadata)[15] <- "Distance_to_centroid"

#Test whether the beta-dispersivity across Time points is different
permutest(box.dispersion)
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     3 0.012429 0.0041429 3.4269    999  0.026 *
# Residuals 36 0.043521 0.0012089 

## Kruskal-wallis
bla3 <- kruskal(metadata$Distance_to_centroid, metadata$Time_point, group = T, p.adj = "BH")
bla.groups3 <- bla3$groups
# metadata$Distance_to_centroid groups
# T0                          30.4      a
# T2                          22.7     ab
# T3                          14.5      b
# T1                          14.4      b

#Plot the distance to centroid values
pdf("Plots/Bac_Beta_dispersivity_10022022.pdf", 
    width = 5, height = 4)
boxplot(box.dispersion, xlab ="Time Point", ylab = "Distance to centroid")
dev.off()


### Plot Bray-curtis distance comparisons ###
rownames(box.asv)== metadata$Sample # re-calculate distance

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
  labs( x= NULL, y = "Bray-Curtis distances \n for bacterial/archaeal community") +
  scale_x_discrete(breaks =c("first", "second", "third"),
                   labels =c("T0 vs. T1",
                             "T0 vs. T2",
                             "T0 vs. T3")) +
  scale_y_continuous(limits = c(0.1,0.5), breaks= seq(0.1, 0.5, 0.1)) +
  theme_classic()


bla2 <- kruskal(dist_com$value, dist_com$comparison, group = T, p.adj = "BH")
bla.groups2 <- bla2$groups
#            dist_com$value groups
# second        163.440      a
# third         159.835      a
# first         128.225      b

ggsave(file = "Plots/Bac_BrayCurtis_Dist_Com_10032022.pdf", 
       plot = dist.comp.plot,
       width = 100,
       height = 100,
       units ="mm")


write.table(metadata, file = "Results/Bac_Box_Final_Diversity_09262022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

