library(ggplot2)
library(plyr)
library(RColorBrewer)

#define the color pallet
cpal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#FF8B8B", "#0072B2", "#D55E00", "#CC79A7")

#you can also use RColorBrewer to make pallets
#cpal <- (brewer.pal(n = 8, name = "Dark2"))

#cds <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_cds_lengths.txt", header = FALSE)
#exon <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_exon_lengths.txt", header = FALSE)
#intron <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_intron_lengths.txt", header = FALSE)
#fp_utr <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_five_prime_utr_lengths.txt", header = FALSE)
#tp_utr <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_three_prime_utr_lengths.txt", header = FALSE)
#trans <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_transcript_lengths.txt", header = FALSE)

#make a data frame for each feature type
cds <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_longest_cds_cds_lengths.txt", header = FALSE)
exon <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_longest_cds_exon_lengths.txt", header = FALSE)
intron <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_longest_cds_intron_lengths.txt", header = FALSE)
fp_utr <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_longest_cds_five_prime_utr_lengths.txt", header = FALSE)
tp_utr <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_longest_cds_three_prime_utr_lengths.txt", header = FALSE)
trans <- read.table("/Users/mcampbel/projects/CSHL_projects/W22/gene_feature_distributions/W22_round1_longest_cds_transcript_lengths.txt", header = FALSE)

#add a feature lable column to each data frame so they can be merged later
cds$Feature <- 'cds'
exon$Feature <- 'exon'
intron$Feature <- 'intron'
fp_utr$Feature <- 'fp_utr'
tp_utr$Feature <- 'tp_utr'
trans$Feature <- 'trans'

#merge the individual data frames into one data frame for plotting
FeatureLengths <- rbind(cds, exon, intron, fp_utr, tp_utr, trans)

#calculate the means and medians of the features using ddply
#ddply will take a data frame spllitt it up and work on the peices and retrun another dataframe
#with the results
FeatureMeans <- ddply(FeatureLengths, .variables = .(Feature), summarise, Feature.mean=mean(V1))
FeatureMedians <- ddply(FeatureLengths, .variables = .(Feature), summarise, Feature.median=median(V1))

#make the plot
ggplot(FeatureLengths, aes(FeatureLengths$V1, fill = Feature)) + 
  #plot type. alpha controls the tranperency  
  geom_density(alpha = 0.7) + 
  #define the x and y axes
  xlim(0, 2000) +
  ylim(0,0.009) +
  #add the median lines  
  geom_vline(data=FeatureMedians, aes(xintercept=Feature.median, colour=Feature), linetype="12345678",size=0.5) +
  #add data lables
  labs(x = "Length in base pairs", y = "Density", title = "Feature Length Distribution") +
  #chage the look of the plot
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5, size = 20)) +
  #use the color pallet defined abov
  scale_fill_manual(values=cpal) +
  scale_colour_manual(values=cpal)
  

