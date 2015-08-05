########################################
# plot diagnosis of HeLa_input samples #
########################################

# directory for this analysis
setwd("/TEMP_DDN/users/gfilion/rlim/VirusAlignment/Ranalysis/H1hESC_input")

## get the threshold (max threshold) ##
############################################################################################
# Idea: by sampling (with replacement) a.k.a Bootstrapping the number of viral infections  #
############################################################################################

# for reproducibility of sampling procedure
set.seed(123)

# obtain the vector of virus size for all viruses in our viral genome library 
virus_size <- scan("../virus_size.lst", what=0)
head(virus_size)
length(virus_size)

# get the infection scores
infectionScores_001 <- read.delim("sortedNormalizedInfection-unmapped-hg18-H1-008.fastq.map.fastq.map", header=F)

# obtain the number of viral infections
virus_infected_001 <- read.delim("countedInfection-unmapped-hg18-H1-008.fastq.map.fastq.map", header=F)
head (virus_infected_001)
viruses_001 <- sum(virus_infected_001$V3)

# boostrapping to estimate the max_threshold
thresholds <- NA
for (i in 1:1000) { thresholds[i] = max(tabulate(sample(size=viruses_001, 1:length(virus_size), 
                                                        prob=virus_size, replace=TRUE))/virus_size) }

# round up the threshold
max_threshold <- signif(max(thresholds), digits=2)
max_threshold

fix
# generate the diagnosis plot
plot(sample(infectionScores_001$V2[1:300]), type = 'h')
abline(h=max_threshold, col=2, lty=2)

library(ggplot2)
plotSample_001 <- infectionScores_001 [1:30,]
head(plotSample_001)
nrow(plotSample_001)
colnames(plotSample_001) <- c("Viruses", "Scores")

annot_virus <-  subset(plotSample_001, plotSample_001[,2]> max_threshold)
anno <- data.frame(x = annot_virus$Viruses, y = annot_virus$Scores, lab = annot_virus$Viruses)
anno

gp <- ggplot(plotSample_001) + geom_bar(aes(x=Viruses, y=Scores, group=Viruses, fill=Viruses, alpha=Scores, position="dodge")) +
  opts(axis.text.x = theme_blank(),axis.ticks = theme_blank())+
  opts(panel.background = theme_rect(fill = 'white'))+
  geom_hline(yintercept=max_threshold, color="red", lty=2, lwd=2)+
  geom_text(data = anno, aes(x,y,label = lab, angle = 30, family="mono"))


# Here is your plots
gp
annot_threshold <- paste("Threshold:", max_threshold)
gp + guides(fill=FALSE)+coord_cartesian(ylim=c(0.0, max(plotSample_001$Scores)+0.01)) +
  annotate("text", label = annot_threshold, x = 27, y = max_threshold+0.001, size = 5, colour = "red")


#custom color palette
grad <- colorRampPalette(c('red', 'yellow'))(256)
my.col=grad[plotSample_001$Scores*2560]