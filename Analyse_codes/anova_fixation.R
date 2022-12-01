################# do two way anova for eye-movement data

## set path
setwd("Z:/Semantic/Results/sv/Behav/") 

## load in data
my_data <- read.csv("eye_movements.csv")


### get needed packages
# effet size -- Cohen's d for t-test
install.packages("effsize")  
library(effsize)
install.packages("lsr")
library(lsr)
# install raincloudplots
packages <- c("ggplot2", "dplyr", "lavaan", "plyr", "cowplot", "rmarkdown",
              "readr", "caTools", "bitops")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}
library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)
source("R_rainclouds.R")
source("summarySE.R")
source("simulateData.R")
if (!require(remotes)) {
  install.packages("remotes")
}
remotes::install_github('jorvlan/raincloudplots')
library(raincloudplots)
install.packages('svglite')




################################ first fixation durations #################################
## Convert voriables as factors
my_data$congruency <- factor(my_data$congruency, 
                       levels = c(1, 2),
                       labels = c("low", "high"))
my_data$location <- factor(my_data$location, 
                          levels = c(1, 2),
                          labels = c("pre", "targ"))
my_data$subject <- factor(my_data$subject)


# pre-target: t(33)=0.839, p=0.408, d = 0.144
data_pre <- subset(my_data, location == "pre")
t_pre <- t.test(duration ~ congruency, data = data_pre, paired = TRUE, alternative = "two.sided")
t_pre
# cohen's d 
pre_low <- subset(data_pre, congruency == "low")
pre_high <- subset(data_pre, congruency == "high")
cohensD(pre_low$duration, pre_high$duration, method = "paired")

# target: t(33)=5.994, p=9.828e-07, d = 1.028
data_targ <- subset(my_data, location == "targ")
t_targ <- t.test(duration ~ congruency, data = data_targ, paired = TRUE, alternative = "two.sided")
t_targ
# cohen's d
targ_low <- subset(data_targ, congruency == "low")
targ_high <- subset(data_targ, congruency == "high")
cohensD(targ_low$duration, targ_high$duration, method = "paired")

###### raincloud plot ######
# prepare the data
n = max(as.numeric(my_data$subject)) # number of subject
df_2x2 <- data_2x2(
  array_1 = my_data$duration[1:n],
  array_2 = my_data$duration[(n+1):(2*n)],
  array_3 = my_data$duration[(2*n+1):(3*n)],
  array_4 = my_data$duration[(3*n+1):(4*n)],
  labels = (c('incongruent','congruent')),
  jit_distance = .09,
  jit_seed = 321,
  spread_x_ticks = TRUE)
# plot
raincloud_2 <- raincloud_2x2_repmes(
  data = df_2x2,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  line_color = 'gray',
  line_alpha = .3,
  size = 1,
  alpha = .6,
  spread_x_ticks = TRUE) +
scale_x_continuous(breaks=c(1,2,3,4), labels=c("", "","", ""), limits=c(0,5)) +
  xlab("Pre-target Target") +
  ylab("First fixation durations (ms)") +
  theme_classic(base_size=7, base_family='Arial')
raincloud_2
ggsave('Z:/Semantic/Results/sv/Behav/FirstFixationDuration.svg', width=5, height=3.6)



############################### run regression probability #######################################
# pre-target: t(33)=7.8294, p=5.042e-09, d=1.342728
data_pre <- subset(my_data, location == "pre")
t_pre <- t.test(regression ~ congruency, data = data_pre, paired = TRUE, alternative = "two.sided")
t_pre
# cohen's d 
pre_low <- subset(data_pre, congruency == "low")
pre_high <- subset(data_pre, congruency == "high")
cohensD(pre_low$regression, pre_high$regression, method = "paired")

# target: t(33)=9.1324, p=1.494e-10, d=1.566194
data_targ <- subset(my_data, location == "targ")
t_targ <- t.test(regression ~ congruency, data = data_targ, paired = TRUE, alternative = "two.sided")
t_targ
# cohen's d
targ_low <- subset(data_targ, congruency == "low")
targ_high <- subset(data_targ, congruency == "high")
cohensD(targ_low$regression, targ_high$regression, method = "paired")


### plot
# prepare the data
n = max(as.numeric(my_data$subject)) # number of subject
df_2x2 <- data_2x2(
  array_1 = 100*(my_data$regression[1:n]),
  array_2 = 100*(my_data$regression[(n+1):(2*n)]),
  array_3 = 100*(my_data$regression[(2*n+1):(3*n)]),
  array_4 = 100*(my_data$regression[(3*n+1):(4*n)]),
  labels = (c('incongruent','congruent')),
  jit_distance = .09,
  jit_seed = 321,
  spread_x_ticks = TRUE)
# plot
raincloud_2 <- raincloud_2x2_repmes(
  data = df_2x2,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  line_color = 'gray',
  line_alpha = .3,
  size = 1,
  alpha = .6,
  spread_x_ticks = TRUE) +
  scale_x_continuous(breaks=c(1,2,3,4), labels=c("", "","", ""), limits=c(0,5)) +
  xlab("Pre-target Target") +
  ylab("Regression probability (%)") +
  theme_classic(base_size=11, base_family='Arial')
raincloud_2
ggsave('Z:/Semantic/Results/sv/Behav/RegressionProb.svg', width=5, height=3.6)



################################## run total gaze ####################################
# pre-target: t(33)=5.7769, p=1.861e-06, d=0.9907328
data_pre <- subset(my_data, location == "pre")
t_pre <- t.test(totalgaze ~ congruency, data = data_pre, paired = TRUE, alternative = "two.sided")
t_pre
# cohen's d 
pre_low <- subset(data_pre, congruency == "low")
pre_high <- subset(data_pre, congruency == "high")
cohensD(pre_low$totalgaze, pre_high$totalgaze, method = "paired")

# target: t(33)=10.546, p=4.202e-12, d=1.808681
data_targ <- subset(my_data, location == "targ")
t_targ <- t.test(totalgaze ~ congruency, data = data_targ, paired = TRUE, alternative = "two.sided")
t_targ
# cohen's d
targ_low <- subset(data_targ, congruency == "low")
targ_high <- subset(data_targ, congruency == "high")
cohensD(targ_low$totalgaze, targ_high$totalgaze, method = "paired")


### plot
# prepare the data
n = max(as.numeric(my_data$subject)) # number of subject
df_2x2 <- data_2x2(
  array_1 = my_data$totalgaze[1:n],
  array_2 = my_data$totalgaze[(n+1):(2*n)],
  array_3 = my_data$totalgaze[(2*n+1):(3*n)],
  array_4 = my_data$totalgaze[(3*n+1):(4*n)],
  labels = (c('incongruent','congruent')),
  jit_distance = .09,
  jit_seed = 321,
  spread_x_ticks = TRUE)
# plot
raincloud_2 <- raincloud_2x2_repmes(
  data = df_2x2,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  line_color = 'gray',
  line_alpha = .3,
  size = 1,
  alpha = .6,
  spread_x_ticks = TRUE) +
  scale_x_continuous(breaks=c(1,2,3,4), labels=c("", "","", ""), limits=c(0,5)) +
  xlab("Pre-target Target") +
  ylab("Total gaze duration (ms)") +
  theme_classic(base_size=11, base_family='Arial')
raincloud_2
ggsave('Z:/Semantic/Results/sv/Behav/TotalGaze.svg', width=5, height=3.6)





############# effect size for coherence t-test ###############################
## set path
setwd("Z:/Semantic/Results/sv/Coh/") 

my_data <- read.csv("coherence.csv")
my_data$congruency <- factor(my_data$congruency, 
                          levels = c(1, 2),
                          labels = c("low", "high"))
my_data$location <- factor(my_data$location, 
                          levels = c(1, 2),
                          labels = c("pre", "targ"))
my_data$subject <- factor(my_data$subject)
# pre-target: t(28)=-2.561, p=0.016, d = 0.476
data_pre <- subset(my_data, location == "pre")
t_pre <- t.test(coherence ~ congruency, data = data_pre, paired = TRUE, alternative = "two.sided")
t_pre
# cohen's d
pre_low <- subset(data_pre, congruency == "low")
pre_high <- subset(data_pre, congruency == "high")
cohensD(pre_low$coherence, pre_high$coherence, method = "paired")

### target: t(28)=0.499, p=0.622, d = 0.093
data_targ <- subset(my_data, location == "targ")
t_targ <- t.test(coherence ~ congruency, data = data_targ, paired = TRUE, alternative = "two.sided")
t_targ
# cohen's d
targ_low <- subset(data_targ, congruency == "low")
targ_high <- subset(data_targ, congruency == "high")
cohensD(targ_low$coherence, targ_high$coherence, method = "paired")






