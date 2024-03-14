# nf-core/crisprseq benchmarking with spike-in samples

### Inport all R libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

## Benchmarking of samples

"Calculate the detected indel percentage with nf-core/crisprseq vs samples from the CRISPR-Analytics publication and Connelly et.al. publication (analysed with crisp.py)"

### Read nf-core/crisprseq results

indel_crispra_samples <- read.csv("benchmarking_results/crisprseq_spikes/multiqc_edition_bargraph_crispra-samples.txt", sep="\t", stringsAsFactors = FALSE)
indel_crisppy_samples <- read.csv("benchmarking_results/crisprseq_spikes/multiqc_edition_bargraph_crisppy-samples.txt", sep="\t", stringsAsFactors = FALSE)

### Add percentages detected by original publication

# Manually add CRISPR-A analysed and expected percentage
indel_manual_data_crispra <- data.frame("Sample" = c("Spikes-high-x25-a","Spikes-high-x25-b","Spikes-high-x30-a","Spikes-high-x30-b","Spikes-high-x35-a","Spikes-high-x35-b",
                                             "Spikes-low-x25-a","Spikes-low-x30-a","Spikes-low-x30-b","Spikes-low-x35-a","Spikes-low-x35-b", "Spikes-no-mix-x25"),
                                "crispra_indel" = c(31, 32, 31, 32, 33, 34, 30, 31, 32, 33, 33, 99),
                                "expected" = c(20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 100)
)
# Manually add crisp.py percentages
indel_manual_data_crisppy <- data.frame("Sample" = c("PC15_hF9_0-A", "PC15_hF9_0-B", "PC15_hF9_0-C", "PC15_hF9_100-A", "PC15_hF9_100-B", "PC15_hF9_100-C",
                                                     "PC15_hF9_20-A", "PC15_hF9_20-B", "PC15_hF9_20-C", "PC15_hF9_40-A", "PC15_hF9_40-B", "PC15_hF9_40-C", 
                                                     "PC15_hF9_60-A", "PC15_hF9_60-B", "PC15_hF9_60-C", "PC15_hF9_80-A", "PC15_hF9_80-B", "PC15_hF9_80-C", 
                                                     "PC19_hADAR1_0-A", "PC19_hADAR1_0-B", "PC19_hADAR1_0-C", "PC19_hADAR1_100-A", "PC19_hADAR1_100-B", 
                                                     "PC19_hADAR1_100-C", "PC19_hADAR1_20-A", "PC19_hADAR1_20-B", "PC19_hADAR1_20-C", "PC19_hADAR1_40-A", 
                                                     "PC19_hADAR1_40-B", "PC19_hADAR1_40-C", "PC19_hADAR1_60-A", "PC19_hADAR1_60-B", "PC19_hADAR1_60-C", 
                                                     "PC19_hADAR1_80-A", "PC19_hADAR1_80-B", "PC19_hADAR1_80-C"),
                                        "crisppy_indel" = c(3.4,2.8,3.4,94.9,94.4,95,19.6,22,19.7,38.9,39.8,38.7,59.4,56.9,56.5,74.7,75.5,76.2,
                                                            0.3,2.2,4.1,90.6,90.3,91.5,14.8,21.3,21.8,42.7,36.8,40.7,56.1,58.5,52.8,67,68.4,69.3)
                                        )

#indel <- rbind(indel_crispra_samples, indel_crisppy_samples)
#indel <- left_join(indel, indel_manual_data_crispra)
#indel <- left_join(indel, indel_manual_data_crisppy)
indel_crispra <- left_join(indel_crispra_samples, indel_manual_data_crispra)
indel_crisppy <- left_join(indel_crisppy_samples, indel_manual_data_crisppy)

### Format data to plot

indel_crispra_formatted <- indel_crispra %>%
  # remove hyphens
  mutate_all(str_replace_all, "no-mix", "nomix") %>%
  # Mutate values to numeric
  transform(Wt = as.numeric(Wt)) %>%
  transform(Template.based = as.numeric(Template.based)) %>%
  transform(Delins = as.numeric(Delins)) %>%
  transform(Ins_inframe = as.numeric(Ins_inframe)) %>%
  transform(Ins_outframe = as.numeric(Ins_outframe)) %>%
  transform(Dels_inframe = as.numeric(Dels_inframe)) %>%
  transform(Dels_outframe = as.numeric(Dels_outframe)) %>%
  transform(crispra_indel = as.numeric(crispra_indel)) %>%
  transform(expected_percent = as.numeric(expected)) %>%
  # Calculate indel percentage
  transform(Total = Wt + Template.based + Delins + Ins_inframe + Ins_outframe + Dels_inframe + Dels_outframe) %>%
  transform(WT_perc = Wt / Total * 100) %>%
  transform(Indel_perc = (Delins + Ins_inframe + Ins_outframe + Dels_inframe + Dels_outframe) / Total * 100) %>%
  # Add column to specify the experiment
  mutate(origin_run = "crispra") %>%
  # Exclude sample Spikes-low-x25-b which was removed from the CRISPR-A analysis due to too low concentration
  filter(Sample != "Spikes-low-x25-b") %>%
  # Select relevant columns & format for plot
  select(origin_run, crispra_indel, Indel_perc, expected_percent) %>%
  pivot_longer(cols=c(Indel_perc, crispra_indel), names_to="analysis", values_to="percent")

indel_crisppy_formated <- indel_crisppy %>%
  # Separate by sample and expected percentage
  separate(Sample, c('Sample', 'reference', 'percentage', "replicate"), sep='_|-', fill="right") %>%
  # Calcualte indel percentage
  transform(Total = Wt + Template.based + Delins + Ins_inframe + Ins_outframe + Dels_inframe + Dels_outframe) %>%
  transform(WT_perc = Wt / Total * 100) %>%
  transform(Indel_perc = (Delins + Ins_inframe + Ins_outframe + Dels_inframe + Dels_outframe) / Total * 100) %>%
  # Mutate expected percentage to numeric
  transform(expected_percent = as.numeric(percentage)) %>%
  # Add column to specify the experiment
  mutate(origin_run = ifelse(Sample == "PC15", "crisppy_A", "crisppy_B")) %>%
  # Select relevant columns  & format for plot
  select(origin_run, crisppy_indel, Indel_perc, expected_percent) %>%
  pivot_longer(cols=c(Indel_perc, crisppy_indel), names_to="analysis", values_to="percent")

indel_all <- rbind(indel_crispra_formatted, indel_crisppy_formated)
indel_all <- indel_all %>% mutate(expected_percent = as.factor(expected_percent))

### Plot

indel_boxplot <- ggplot(indel_all, aes(x=expected_percent, y=percent, color=analysis)) + 
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  theme(legend.position="top") +
  labs(x="Expected indel percentage", y="Detected indel percentage") +
  scale_color_manual(name = "Analysis pipeline", labels = c("crisp.py", "crispr-a", "crisprseq"), values = c("#5ab4ac", "#d582f5", "#fcc555")) +
  facet_grid(~origin_run)
indel_boxplot

ggsave(plot=indel_boxplot, filename = "expected_vs_detected_analysis_pipeline.png", device="png",width = 20, height = 10, units="cm")
