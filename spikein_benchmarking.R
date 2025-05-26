# nf-core/crisprseq benchmarking with spike-in samples

### Inport all R libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)

## Benchmarking of samples

"Calculate the detected indel percentage with nf-core/crisprseq vs CRISPR-Analytics using samples from the CRISPR-Analytics publication (Sanvicente-García et.al.) and Connelly et.al. publication"

### Read nf-core/crisprseq results

indel_sanvi_samples_crisprseq_analysis <- read.csv("benchmarking_results/crisprseq_spikes/multiqc_edition_bargraph_crispra-samples.txt", sep="\t", stringsAsFactors = FALSE)
indel_connelly_samples_crisprseq_analysis <- read.csv("benchmarking_results/crisprseq_spikes/multiqc_edition_bargraph_crisppy-samples.txt", sep="\t", stringsAsFactors = FALSE)

# Process crispr-A results (Connelly samples analysed by crispra)
files <- Sys.glob("benchmarking_results/crisprseq_spikes/crispy_spikein_crispra/*")
csv_list <- lapply(files, function(x) read.csv(x, , sep=",", stringsAsFactors = FALSE, ))
final_df <- bind_rows(lapply(seq_along(csv_list), function(index) {
  csv_list[[index]] %>%
    mutate(sample = paste0("sample_", index)) # Add a 'sample' column with the list index
}), .id = "sample")
final_df <- final_df %>%
  select(sample, classes, counts) %>% # Select relevant columns
  pivot_wider(names_from = classes, values_from = counts, values_fill = 0)
sample_names <- c("PC19_hADAR1_100-A", "PC19_hADAR1_100-B", "PC19_hADAR1_100-C", "PC19_hADAR1_80-A", "PC19_hADAR1_80-B", "PC19_hADAR1_80-C", "PC19_hADAR1_60-A",
                  "PC19_hADAR1_60-B", "PC19_hADAR1_60-C", "PC19_hADAR1_40-A", "PC19_hADAR1_40-B", "PC19_hADAR1_40-C", "PC19_hADAR1_20-A", "PC19_hADAR1_20-B",
                  "PC19_hADAR1_20-C", "PC19_hADAR1_0-A", "PC19_hADAR1_0-B", "PC19_hADAR1_0-C", "PC15_hF9_100-A", "PC15_hF9_100-B", "PC15_hF9_100-C", "PC15_hF9_80-A",
                  "PC15_hF9_80-B", "PC15_hF9_80-C", "PC15_hF9_60-A", "PC15_hF9_60-B", "PC15_hF9_60-C", "PC15_hF9_40-A", "PC15_hF9_40-B", "PC15_hF9_40-C",
                  "PC15_hF9_20-A", "PC15_hF9_20-B", "PC15_hF9_20-C", "PC15_hF9_0-A", "PC15_hF9_0-B", "PC15_hF9_0-C")

indel_connelli_samples_crispra_analysis <- as.data.frame(final_df) %>% 
  mutate(sample = sample_names) %>%
  select(sample, Wt, Indels)
# Percentage of Wt reads and reads with an indel (insertion, deletion or delin)

# Format (origin_run, expected_percent, analysis, percent)
indel_connelli_samples_crispra_analysis_formatted <- indel_connelli_samples_crispra_analysis %>% 
  mutate(analysis = rep("crispr-a", 36)) %>%
  # Separate sample name to obtain expected percentage
  separate(sample, c('Sample', 'reference', 'percentage', "replicate"), sep='_|-', fill="right") %>%
  transform(expected_percent = as.numeric(percentage)) %>%
  # Add column to specify the experiment
  mutate(origin_run = ifelse(Sample == "PC15", "PC15_edited", "PC19_edited")) %>%
  # Select relevant columns  & format for plot
  mutate(percent = Indels) %>%
  select(origin_run, expected_percent, analysis, percent)


### Add percentages detected by original publication

# Manually add CRISPR-A analysed and expected percentage
indel_manual_data_crispra <- data.frame("Sample" = c("Spikes-high-x25-a","Spikes-high-x25-b","Spikes-high-x30-a","Spikes-high-x30-b","Spikes-high-x35-a","Spikes-high-x35-b",
                                             "Spikes-low-x25-a","Spikes-low-x30-a","Spikes-low-x30-b","Spikes-low-x35-a","Spikes-low-x35-b", "Spikes-no-mix-x25"),
                                "crispra_indel" = c(31, 32, 31, 32, 33, 34, 30, 31, 32, 33, 33, 99),
                                "expected" = c(20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 100)
)

indel_sanvi <- left_join(indel_sanvi_samples_crisprseq_analysis, indel_manual_data_crispra)
indel_connelli <- indel_connelly_samples_crisprseq_analysis

### Format data to plot

indel_sanvi_formatted <- indel_sanvi %>%
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
  mutate(origin_run = "SDM") %>%
  # Exclude sample Spikes-low-x25-b which was removed from the CRISPR-A analysis due to too low concentration
  filter(Sample != "Spikes-low-x25-b") %>%
  # Select relevant columns & format for plot
  select(origin_run, crispra_indel, Indel_perc, expected_percent) %>%
  pivot_longer(cols=c(Indel_perc, crispra_indel), names_to="analysis", values_to="percent") %>%
  # Specify analysis method
  mutate(analysis = ifelse(analysis == "Indel_perc", "crisprseq", "crispr-a"))
  

indel_connelli_formated <- indel_connelli %>%
  # Separate by sample and expected percentage
  separate(Sample, c('Sample', 'reference', 'percentage', "replicate"), sep='_|-', fill="right") %>%
  # Calcualte indel percentage
  transform(Total = Wt + Template.based + Delins + Ins_inframe + Ins_outframe + Dels_inframe + Dels_outframe) %>%
  transform(WT_perc = Wt / Total * 100) %>%
  transform(Indel_perc = (Delins + Ins_inframe + Ins_outframe + Dels_inframe + Dels_outframe) / Total * 100) %>%
  # Mutate expected percentage to numeric
  transform(expected_percent = as.numeric(percentage)) %>%
  # Add column to specify the experiment
  mutate(origin_run = ifelse(Sample == "PC15", "PC15_edited", "PC19_edited")) %>%
  # Select relevant columns  & format for plot
  select(origin_run, Indel_perc, expected_percent) %>%
  pivot_longer(cols=c(Indel_perc), names_to="analysis", values_to="percent") %>%
  # Analysis method
  mutate(analysis = rep("crisprseq", 36))

indel_all <- rbind(indel_sanvi_formatted, indel_connelli_formated, indel_connelli_samples_crispra_analysis_formatted)
indel_all <- indel_all %>% mutate(expected_percent = as.factor(expected_percent))

### Plot

indel_boxplot <- ggplot(indel_all, aes(x=expected_percent, y=percent, color=analysis)) + 
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  theme(legend.position="top") +
  labs(x="Expected indel percentage", y="Detected indel percentage") +
  scale_color_manual(name = "Analysis pipeline", labels = c("CRISPR-A", "nf-core/crisprseq"), values = c("#d582f5", "#fcc555")) +
  facet_grid(~origin_run, scales = "fixed") +
  theme(
    panel.spacing = unit(1, "lines"), # Adjust spacing between panels
    axis.line.y = element_line() # replicate
  )
indel_boxplot

ggsave(plot=indel_boxplot, filename = "expected_vs_detected_analysis_pipeline.png", device="png",width = 20, height = 10, units="cm")

# Don't separate by experiment

indel_boxplot_all <- ggplot(indel_all, aes(x=expected_percent, y=percent, color=analysis)) + 
  geom_boxplot() +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  theme(legend.position="top") +
  labs(x="Expected indel percentage", y="Detected indel percentage") +
  scale_color_manual(name = "Analysis pipeline", labels = c("crispr-a", "crisprseq"), values = c("#d582f5", "#fcc555")) 
indel_boxplot_all

ggsave(plot=indel_boxplot_all, filename = "expected_vs_detected_analysis_pipeline_allsamples.png", device="png",width = 20, height = 10, units="cm")


