## Load library
library(ggplot2)
library(tidyverse)
library(reshape2)


## Import Faprotax file
fap <- read.table("Results/FAPROTAX/SBAR_16S_Mapping_faprotax_results_propor_09272022.txt", 
                         sep = "\t", header = T, row.names = 1)
dim(fap)

fap.sub <- fap[, c(1:14,31:32, 38, 45,46, 48, 39, 64,58, 60:62,  
                        19,26, 29,52, 63)]


## Import Faprotax file
fguild <- read.table("Results/SBAR_ITS_Mapping_funguild_results_propor_09272022.txt", 
                  sep = "\t", header = T, row.names = 1)
dim(fguild)

fguild.sub <- fguild[, c(1:13,29, 31, 41, 45, 49, 57)]

##manipulate rows to match rownames for 38 samples
fap.sub <- fap.sub[c(1:18, 20:24,26:40, 19, 25),]
rownames(fap.sub[1:38,]) == rownames(fguild.sub)

## add 2 rows to fungal data
rownames(fguild.sub[39,]) <- "S_035"
rownames(fguild.sub[40,]) <- "S_045"

rownames(fap.sub) == rownames(fguild.sub)

## Combine the two datasets
com <- cbind(fap.sub, fguild.sub[, c(14:19)])

#Remove NAs
com[is.na(com)] <- 0
com[c(39:40),]

com$Time_point <- factor(com$Time_point,
                                 levels = c("T0", "T1", "T2", "T3"))

### Calculate average for Time-point
com.heatmap <- melt(com[,c(11,15:37)])%>%
  group_by(Time_point, variable) %>%
  summarize(mean_prop = mean(value), .groups = "drop")
 
colnames(com.heatmap) <- c("Time_Point" , "Function", "mean_prop")

## Set the order of labels backwards for correct order output
# com.heatmap$Function <- factor(com.heatmap$Function,
#                          levels = c())     

com.heatmap$round.mean <- round(com.heatmap$mean_prop, digits = 3)

fig3 <- ggplot(data = com.heatmap, aes(x = Time_point , 
                                     y = forcats::fct_rev(variable), 
                                     fill = round.mean)) +
  geom_tile(color = "gray") +
  geom_text(aes(x = Time_point, 
                y = variable, label = round.mean), size = 1.8)+
  coord_fixed(ratio = 0.6) +
  scale_fill_gradientn(colors = c( 
    "#FFFFCC","#FFFF99",
    "#FFFF33", "#F1C232",
    "#F1A832",
    "red"), 
    breaks = c(0,0.04, 0.08,0.12,0.16, 0.2))+
  scale_x_discrete(position = "top") +
  ylab("") +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank()) +
  theme(axis.ticks = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.position = "bottom",
        legend.direction = "horizontal", 
        legend.box = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,-10,-10)) +
  guides(fill = guide_colorbar(title = "Average proportion", 
                               title.position = "top"))

ggsave(file = "Figures/Heatmap_Function_09282022.pdf", 
       plot = fig3,
       width = 120,
       height = 140,
       units ="mm")

