setwd("C:/JuliaWork/ION_Study/Rprojects/") #directory with data tables
library(tidyverse)
library(NHANES)


data1 <- read.csv( "dice_split_half_within_subj/sub-0011_dice_vector_top_25pct.csv", header=FALSE, sep=",") 
data2 <- read.csv( "dice_split_half_between_subj/sub-0011_dice_overlap_vector_all_subj.csv", header=FALSE, sep=",")


#colnames() <- c('PB0007','PB0009', 'PB0010', 'PB0011', 'PB0015', 'PB0016', 'PB0017', 'PB0018', 'PB0019', 'PB0020')
colnames(data1) <- 'half2'
colnames(data2) <- c('PB0007', 'PB0009', 'PB0010', 'PB0011', 'PB0015', 'PB0016', 'PB0017', 'PB0018', 'PB0019')

data_short=cbind(data1,data2)

data <-gather(data_short, SUB) %>% 
  mutate(condition=ifelse(SUB=='half2', 'self', 'other'))

model <- lm(value ~ condition, data = data) 
sink(file = "dice_split_half_between_subj/sub-0011_lm_output_self_other.txt")
summary(model)
sink(file = NULL)



#colors <- c("#2B4231", "#228833", "#939D5C", "#DC9B41", "#CA5B48", "#E19790", "#AA3377", "#382585", "#56B4E9", "#BBBBBB")
colors <- c("#BBBBBB", "#F0F0F0", "#F0F0F0", "#F0F0F0", "#F0F0F0", "#F0F0F0", "#F0F0F0", "#F0F0F0", "#F0F0F0", "#F0F0F0")



data %>% 
  ggplot(aes(
    x=SUB, 
    y=value,
  )) +
  geom_violin(aes(fill = SUB), trim = TRUE) + 
  geom_boxplot(width = 0.2)+
  ylab("dice coefficient top 25% values") +
  #xlab("participant") +
  ggtitle("PB0020") +
  ylim(0, 0.75) +
  theme_minimal() +
  theme(
    #strip.text.y = element_text(size = 12),
    #text = element_text(size = 12),
    #legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 22, face='bold'),
    axis.text.x=element_text(angle = 90, size=18, color="black", vjust=0.5),
    axis.text.y=element_text(size=18, color="black"),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    #legend.position=c(0.8, 0.93),
    legend.position="none",
    #panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  scale_fill_manual(values=colors)

