setwd("C:/JuliaWork/ION_Study/Rprojects/dice_split_half_within_subj") #directory with data tables
library(tidyverse)
library(NHANES)

Data.in <- lapply(list.files(pattern = "\\.csv$"),read.table,sep=",",header=F)

#combines all of the tables by column into one 
Data.in <- do.call(cbind,Data.in)
#pbsublist=data.frame(SUB=c('PB0007','PB0009', 'PB0010', 'PB0011', 'PB0015', 'PB0016', 'PB0017', 'PB0018', 'PB0019', 'PB0020'))

colnames(Data.in) <- c('PB0007','PB0009', 'PB0010', 'PB0011', 'PB0015', 'PB0016', 'PB0017', 'PB0018', 'PB0019', 'PB0020')

data <-gather(Data.in, SUB)

colors <- c("#2B4231", "#228833", "#939D5C", "#DC9B41", "#CA5B48", "#E19790", "#AA3377", "#382585", "#56B4E9", "#BBBBBB")

data %>% 
  group_by(SUB) %>% 
  summarise(med=mean(value)) 

data %>% 
  summarise(med=mean(value)) 

data %>% 
  ggplot(aes(
    x=SUB, 
    y=value,
  )) +
  geom_violin(aes(fill = SUB), trim = TRUE) + 
  geom_boxplot(width = 0.2)+
  ylab("dice coefficient top 25% values") +
  xlab("participant") +
  ylim(-0.4,0.8) +
  theme_minimal() +
  theme(
    #strip.text.y = element_text(size = 12),
    #text = element_text(size = 12),
    #legend.title = element_blank(),
    axis.text.x=element_text(angle = 90, size=14, color="black", vjust=0.5),
    axis.text.y=element_text(size=14, color="black"),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    #legend.position=c(0.8, 0.93),
    legend.position="none",
    #panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  scale_y_continuous(n.breaks=5) +
  scale_fill_manual(values=colors)

