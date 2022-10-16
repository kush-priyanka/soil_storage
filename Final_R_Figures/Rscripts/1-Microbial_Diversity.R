# Load libraries
library(ggplot2)
library(cowplot)
library(tidyverse)

## Import microbial diversity data for plants only
bac <- read.table("16S_Sequencing_Files/Results/Bac_Box_Final_Diversity_09262022.txt", 
                  header = T, row.names = 1, sep = "\t") #80 14


bac$Time_point <-factor(bac$Time_point,
                       levels = c("T0", "T1", "T2", "T3"))

fungi <- read.table("ITS_Sequencing_Files/Results/SBAR_ITS_Rar_Diversity_09232022.txt", 
                  header = T, row.names = 1, sep = "\t") #80 14


fungi$Time_point <- factor(fungi$Time_point,
                        levels = c("T0", "T1", "T2", "T3"))


#### Plots Bac/Arc Richness + Shannon
names <- c("Number of ASVs (phylotypes)","Shannon index")

brich <- ggplot(bac, aes(x = Time_point,
                            y = Bac_Richness,
                         color = Time_point)) +
    geom_boxplot(linetype = "solid", lwd = 0.2, 
                 outlier.shape = NA, 
                 alpha = 0.75) +
    geom_point(aes(color = Time_point), size =3,
                   position = position_jitter(width = 0.01, height = 0.01)) +
    xlab("") +
    ylab(names[1]) + 
    # scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
    #                   breaks = c("T0","T1", "T2", "T3")) +
   scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                    breaks = c("T0","T1", "T2", "T3")) +
    # scale_color_manual(values = c("#E69F00", "#009E73", "#0072B2", "#D55E00"),
    #                 breaks = c("T0","T1", "T2", "T3")) +
    scale_x_discrete(labels = c("T0" = "T0",
                                "T1" = "T1",
                                "T2" = "T2",
                                "T3" = "T3")) +
  # annotate( "text", x = 1,
  #           y = 2300,
  #           label = "a") +
  # annotate( "text", x = 2,
  #           y = 2600,
  #           label = "a") +
  # annotate( "text", x = 3,
  #           y = 2200,
  #           label = "a") +
  # annotate( "text", x = 4,
  #           y = 2000,
  #           label = "a") +
    theme_classic() +
    theme(axis.text = element_text(size = 7,
                                   color = "black"),
          axis.title.y = element_text(size = 9),
          legend.position = "none") 

bshannon <- ggplot(bac, aes(x = Time_point,
                         y = Bac_Shannon,
                         color = Time_point)) +
  geom_boxplot(linetype = "solid", lwd = 0.2, 
               outlier.shape = NA, 
               alpha = 0.75) +
  geom_point(aes(color = Time_point), size =3,
             position = position_jitter(width = 0.01, height = 0.01)) +
  xlab("") +
  ylab(names[2]) + 
  # scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
  #                   breaks = c("T0","T1", "T2", "T3")) +
  scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                     breaks = c("T0","T1", "T2", "T3")) +
  scale_x_discrete(labels = c("T0" = "T0",
                              "T1" = "T1",
                              "T2" = "T2",
                              "T3" = "T3")) +
  # annotate( "text", x = 1,
  #           y = 7.2,
  #           label = "a") +
  # annotate( "text", x = 2,
  #           y = 6.9,
  #           label = "a") +
  # annotate( "text", x = 3,
  #           y = 6.7,
  #           label = "b") +
  # annotate( "text", x = 4,
  #           y = 6.6,
  #           label = "b") +
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "none") 

#### Plots Fungal Richness + Shannon
frich <- ggplot(fungi, aes(x = Time_point,
                         y = Fungi_Richness,
                         color = Time_point)) +
  geom_boxplot(linetype = "solid", lwd = 0.2, 
               outlier.shape = NA, 
               alpha = 0.75) +
  geom_point(aes(color = Time_point), size =3,
             position = position_jitter(width = 0.01, height = 0.01)) +
  xlab("") +
  ylab(names[1]) + 
  scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                    breaks = c("T0","T1", "T2", "T3")) +
  scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                     breaks = c("T0","T1", "T2", "T3")) +
  scale_x_discrete(labels = c("T0" = "T0",
                              "T1" = "T1",
                              "T2" = "T2",
                              "T3" = "T3")) +
  # annotate( "text", x = 1,
  #           y = 260,
  #           label = "a") +
  # annotate( "text", x = 2,
  #           y = 220,
  #           label = "a") +
  # annotate( "text", x = 3,
  #           y = 185,
  #           label = "a") +
  # annotate( "text", x = 4,
  #           y = 190,
  #           label = "a") +
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "none") 

fshannon <- ggplot(fungi, aes(x = Time_point,
                            y = Fungi_Shannon,
                            color = Time_point,
                            label= Sample)) +
  geom_boxplot(linetype = "solid", lwd = 0.2, 
               outlier.shape = NA, 
               alpha = 0.75) +
  geom_point(aes(color = Time_point), size =3,
             position = position_jitter(width = 0.01, height = 0.01)) +
  geom_text(nudge_x = 0.18, size=2) +
  xlab("") +
  ylab(names[2]) + 
  # scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
  #                   breaks = c("T0","T1", "T2", "T3")) +
  scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                     breaks = c("T0","T1", "T2", "T3")) +
  scale_x_discrete(labels = c("T0" = "T0",
                              "T1" = "T1",
                              "T2" = "T2",
                              "T3" = "T3")) +
  # annotate( "text", x = 1,
  #           y = 3.6,
  #           label = "a") +
  # annotate( "text", x = 2,
  #           y = 2.8,
  #           label = "b") +
  # annotate( "text", x = 3,
  #           y = 2.75,
  #           label = "b") +
  # annotate( "text", x = 4,
  #           y = 2.7,
  #           label = "b") +
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "none") 


### Plot NMDS Bac/arc
bac.centroid <- bac %>%
  group_by(Time_point)%>%
  summarize(Bac_NMDS1= mean(Bac_NMDS1),
            Bac_NMDS2= mean(Bac_NMDS2))

bac.segs <- merge(bac, setNames(bac.centroid, c('Time_point','oNMDS1','oNMDS2')),
              by = 'Time_point', sort = FALSE)

bac.nmds <- ggplot(bac, aes(Bac_NMDS1, Bac_NMDS2))+
  geom_point(aes(color = Time_point), size = 1.5)+
  stat_ellipse(aes(group = Time_point, color = Time_point)) +
  # geom_segment(data = bac.segs,
  #                mapping = aes(xend = oNMDS1, yend = oNMDS2, 
  #                              color = Time_point),
  #              linetype =3) +
  geom_point(data = bac.centroid, size = 2.5, shape= 21, color= "black",
             aes(fill = Time_point)) +
  scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                    breaks = c("T0","T1", "T2", "T3")) +
  scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                    breaks = c("T0","T1", "T2", "T3")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate( "text", x = 0.12,
            y = -0.22,
            label = "Stress = 0.14",
            size = 3) +
  # annotate( "text", x = -0.225,
  #           y = 0.17,
  #           label = bquote(R^2~'(Time point) = 0.18***'),
  #           size = 3) +
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "none")

### Plot Fungal NMDS
fun.centroid <- fungi %>%
  group_by(Time_point)%>%
  summarize(Fungi_NMDS1= mean(Fungi_NMDS1),
            Fungi_NMDS2= mean(Fungi_NMDS2))

fun.segs <- merge(fungi, setNames(fun.centroid, c('Time_point','oNMDS1','oNMDS2')),
                  by = 'Time_point', sort = FALSE)

fun.nmds <- ggplot(fungi, aes(Fungi_NMDS1, Fungi_NMDS2))+
  geom_point(aes(color = Time_point), 
             size = 1.5)+
  stat_ellipse(aes(group = Time_point, color = Time_point)) +
  # geom_segment(data = fun.segs,
  #              mapping = aes(xend = oNMDS1, yend = oNMDS2, 
  #                            color = Time_point), linetype=3) +
  geom_point(data = fun.centroid, size = 2.5, shape= 21, color= "black",
             aes(fill = Time_point)) +
  scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                     breaks = c("T0","T1", "T2", "T3")) +
  scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                    breaks = c("T0","T1", "T2", "T3")) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate( "text", x = 0.1,
            y = -0.1,
            label = "Stress = 0.07",
            size =3) +
  # annotate( "text", x = -0.33,
  #           y = 0.25,
  #           label = bquote(R^2~'(Time point) = 0.36***'),
  #           size =3) +
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "none")


fig1a <- plot_grid(brich, bshannon,frich, fshannon, 
                  ncol = 4, align = "h")

fig1b <- plot_grid(bac.nmds, fun.nmds,
                   ncol =2, align = "h")

fig1 <- plot_grid(fig1a, fig1b,
                   ncol=1, align = "v", rel_widths = c(5, 1))


ggsave(file = "Final_R_Figures/Figures/Fig1_Microbial_Diversity_ellipse_10062022.pdf", 
       plot = fig1b,
       width = 130,
       height = 70,
       units ="mm")

