library(ggbreak)
library(dplyr)
library(ggplot2)
library(reshape)
library(cowplot)

#### Bacterial Taxa ####
##Import Bac relative abundance file at family level
bac.fam <- read.table("Results/RelAbundance/Bac_Family_Relabundance_09262022.txt", 
                      sep = "\t", header = T, row.names = 1)
dim(bac.fam) #40 268


bfam <- melt(bac.fam[,c(11,235,266,100,165,
                            44,198,229,40,105,101,103,
                            19,259,227,169,166,151)])%>%
  group_by(Time_point, variable)%>%
  summarize(mean_rel_abund = 100* mean(value),
            sd_rel_abund = sd(100*(value)), .groups = "drop")

colnames(bfam) <- c("Time_point", "Family", "mean_rel_abund","sd_rel_abund") 
bfam$Time_point <- factor(bfam$Time_point, 
                         levels= c("T3", "T2", "T1", "T0"))

bfam.plot <- ggplot(data = bfam, 
                   aes(y = Family,
                       x = mean_rel_abund, fill = Time_point)) +
  geom_col(width = 0.8,    
           position = position_dodge(0.8)) +
  ylab("")+
  xlab("Mean Relative Abundance (%)") +
  geom_errorbar(aes(xmin = bfam$mean_rel_abund, 
                    xmax = bfam$mean_rel_abund + bfam$sd_rel_abund), 
                width = 0.35, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                    breaks = c("T0","T1", "T2", "T3")) +
  scale_x_continuous(limit =c(0,13), breaks=c(0, 2,4, 6, 8,10, 12) ) +
  # # scale_x_break(c(48, 63),ticklabels = c(63,64,65) ) +
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                      color = "black"),
                        axis.title.y = element_text(size = 8),
                        legend.position = "none") 

#### Fungal Taxa ####
fun.fam <- read.table("Results/RelAbundance/Fungi_Family_Relabundance_09232022.txt", 
                      sep = "\t", header = T, row.names = 1)
dim(fun.fam) #38 108


ffam <- melt(fun.fam[,c(10,47,55,23,43,48,
                        68,29,27,96,65,81)])%>%
  group_by(Time_point, variable)%>%
  summarize(mean_rel_abund = 100* mean(value),
            sd_rel_abund = sd(100*(value)), .groups = "drop")

colnames(ffam) <- c("Time_point", "Family", "mean_rel_abund","sd_rel_abund") 
ffam$Time_point <- factor(ffam$Time_point, 
                          levels= c("T3", "T2", "T1", "T0"))

ffam.plot <- ggplot(data = ffam, 
                    aes(y = Family,
                        x = mean_rel_abund, fill = Time_point)) +
  geom_col(width = 0.8,    
           position = position_dodge(0.8)) +
  ylab("")+
  xlab("Mean Relative Abundance (%)") +
  geom_errorbar(aes(xmin = ffam$mean_rel_abund,
                    xmax = ffam$mean_rel_abund + ffam$sd_rel_abund),
                width = 0.35, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                    breaks = c("T0","T1", "T2", "T3")) +
  scale_x_continuous(limit =c(0,70), breaks=c(0, 2,4, 6, 8,10, 12, 50, 60, 70)) +
  scale_x_break(c(13, 47), scale = 0.5) +
  # scale_x_break(c(48, 63),scale = 0.5, ticklabels = c(63,64) ) +
  # scale_x_break(c(8, 47), scale = 0.3, ticklabels = c(47, 48)) +
  # scale_x_break(c(48, 63),scale = 0.5, ticklabels = c(63,64) ) +
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())

ggsave(file = "Fig2a_Bac_Family_10062022.pdf", 
       plot = bfam.plot,
       width = 70,
       height = 127,
       units ="mm")

ggsave(file = "Fig2b_Fungi_Family_10062022.pdf", 
       plot = ffam.plot,
       width = 128,
       height = 108,
       units ="mm")