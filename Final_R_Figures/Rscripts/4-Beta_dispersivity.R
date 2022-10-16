


### Bac/Arc Plot beta-dispersivity
bcentroid <- ggplot(metadata, aes(x = Time_point,
                         y = Distance_to_centroid,
                         color = Time_point)) +
  geom_boxplot(linetype = "solid", lwd = 0.2, 
               outlier.shape = NA, 
               alpha = 0.75) +
  geom_point(aes(color = Time_point), size =1.5,
             position = position_jitter(width = 0.01, height = 0.01)) +
  ylim(0.05, 0.30) +
  xlab("") +
  ylab("Distance to Centroid") + 
  scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                     breaks = c("T0","T1", "T2", "T3")) +
  scale_x_discrete(labels = c("T0" = "T0",
                              "T1" = "T1",
                              "T2" = "T2",
                              "T3" = "T3")) +
theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "none") 


### Fungi Plot beta-dispersivity
fcentroid <- ggplot(md, aes(x = Time_point,
                                  y = Distance_to_centroid,
                                  color = Time_point)) +
  geom_boxplot(linetype = "solid", lwd = 0.2, 
               outlier.shape = NA, 
               alpha = 0.75) +
  geom_point(aes(color = Time_point), size =1.5,
             position = position_jitter(width = 0.01, height = 0.01)) +
  xlab("") +
  ylab("Distance to Centroid") + 
  scale_color_manual(values = c("#60325F", "#C75C8B", "#519A9A", "#EEC154"),
                     breaks = c("T0","T1", "T2", "T3")) +
  scale_x_discrete(labels = c("T0" = "T0",
                              "T1" = "T1",
                              "T2" = "T2",
                              "T3" = "T3")) +
  ylim(0.05, 0.30)+
  theme_classic() +
  theme(axis.text = element_text(size = 7,
                                 color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "none")


fig4 <- plot_grid(bcentroid, fcentroid,
                   ncol =2, align = "h")

ggsave(file = "Final_R_Figures/Figures/Fig_Supp_Beta_Dispersivity_10062022.pdf", 
       plot = fig4,
       width = 130,
       height = 70,
       units ="mm")