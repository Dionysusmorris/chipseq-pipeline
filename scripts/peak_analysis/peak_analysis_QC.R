# peak_analysis.R
# Written by the Harvard Chan Bioinformatics Core on November 15th, 2024
# This script was written as a demo for the Peak Analysis Workshop
# In order to use this script, the user will need to have downloaded and uncompressed the R project from https://www.dropbox.com/scl/fi/s9mxwd7ttqgjt040m6bm2/Peak_analysis.zip?rlkey=ceqbv4pyx59jxsoa0xoh9l6kb&st=q7rlclil&dl=1

# Load libraries
library(tidyverse)

# Load QC metrics file
#metrics <- read.csv("meta/metrics.csv")
#View(metrics)

# Generate QC metrics file
#path <- "C:/Users/Diony/Documents/Masters/RBIF109/Course Project/Peak_analysis/results/"
#bcbio_templates(type="peakseq", outpath="C:/Users/Diony/Documents/Masters/RBIF109/Course Project/Peak_analysis/results")


# Plot total reads
metrics %>%
  ggplot(aes(x = sample,
             y = total_reads/1e6, 
             fill = antibody)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  coord_flip() +
  scale_y_continuous(name = "Millions of reads") +
  scale_x_discrete(limits = rev) +
  xlab("") + 
  ggtitle("Total reads") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "Antibody") +
  theme(legend.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept=20, color = "black", linetype = "dashed", linewidth=.5)

# Plot mapping rate
metrics %>%
  ggplot(aes(x = sample,
             y = mapped_reads_pct, 
             fill = antibody)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  coord_flip() +
  scale_y_continuous(name = "Percent of Reads Mapped") +
  scale_x_discrete(limits = rev) +
  xlab("") +
  ggtitle("Mapping Rate") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "Antibody") +
  theme(legend.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept=70, color = "black", linetype = "dashed", linewidth=.5)

# Plot Normalized Strand Correlation
metrics %>%
  filter(antibody != "input") %>% 
  ggplot(aes(x = sample,
             y = nsc, 
             fill = antibody)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  coord_flip() +
  scale_y_continuous(name = "NSC coefficient") +
  scale_x_discrete(limits = rev) +
  xlab("") +
  ggtitle("Normalized Strand Cross-Correlation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

# Plot Relative Strand Cross-correlation
metrics %>%
  filter(antibody != "input") %>% 
  ggplot(aes(x = sample,
             y = rsc, 
             fill = antibody)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  coord_flip() +
  scale_y_continuous(name = "RSC Coefficient") +
  scale_x_discrete(limits = rev) +
  xlab("") + 
  ggtitle("Relative Strand Cross-Correlation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

# Plot Fraction in read Peaks FRiP
metrics %>% 
  filter(antibody != "input") %>% 
  ggplot(aes(x = sample,
             y = frip, 
             fill = antibody)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  coord_flip() +
  scale_y_continuous(name = "FRiP") +
  scale_x_discrete(limits = rev) +
  xlab("") + 
  ggtitle("Fraction of Reads in Peaks") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

# Plot Non Redundant Fraction
metrics %>% 
  ggplot(aes(x = sample,
             y = nrf, 
             fill = antibody)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  coord_flip() +
  scale_y_continuous(name = "Non-Redundant Fraction") +
  scale_x_discrete(limits = rev) +
  xlab("") + 
  ggtitle("Non-Redundant Fraction") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "Antibody") +
  theme(legend.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0.9, linetype = "dashed", color="green") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color="orange") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color="red")

# Plot number of peaks
metrics %>% 
  filter(antibody != "input") %>% 
  ggplot(aes(x = sample,
             y = peak_count, 
             fill = antibody)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  coord_flip() +
  scale_y_continuous(name = "Number of Peaks") +
  scale_x_discrete(limits = rev) +
  xlab("") +
  ggtitle("Number of Peaks") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
