# Analysis of CRISPR publicly available data

## Load libraries

library(XML)
library(dplyr)
library(ggplot2)

## Load data

projects <- xmlToDataFrame("public_data/ena_project_20240111-1415.xml")
# headers: "IDENTIFIERS" "NAME" "TITLE" "DESCRIPTION" "SUBMISSION_PROJECT" "RELATED_PROJECTS" "PROJECT_LINKS" "PROJECT_ATTRIBUTES" "UMBRELLA_PROJECT"

find_class <- function(description) {
  # possible classes
  # screening,singlecell,KO,KI,RNAseq,targeted,template,CRISPRi,CRISPRa
  class <- c()
  if (any(grepl("screening", description))) {
    class <- append(class, c("screening"))
  }
  if (any(grepl("functional", description))) {
    class <- append(class, c("screening"))
  }
  if (any(grepl("library", description))) {
    class <- append(class, c("screening"))
  }
  if (any(grepl("libraries", description))) {
    class <- append(class, c("screening"))
  }
  if (any(grepl("pooled", description))) {
    class <- append(class, c("screening"))
  }
  if (any(grepl("screen", description))) {
    class <- append(class, c("screening"))
  }
  if (any(grepl("single cell", description))) {
    class <- append(class, c("singlecell"))
  }
  if (any(grepl("loss of function", description))) {
    class <- append(class, c("KO"))
  }
  if (any(grepl("KO", description))) {
    class <- append(class, c("KO"))
  }
  if (any(grepl("knock-out", description))) {
    class <- append(class, c("KO"))
  }
  if (any(grepl("RNAseq", description))) {
    class <- append(class, c("RNAseq"))
  }
  if (any(grepl("RNA-seq", description))) {
    class <- append(class, c("RNAseq"))
  }
  if (any(grepl("targeted", description))) {
    class <- append(class, c("targeted"))
  }
  if (any(grepl("amplicon", description))) {
    class <- append(class, c("targeted"))
  }
  if (any(grepl("KI", description))) {
    class <- append(class, c("KI"))
  }
  if (any(grepl("template", description))) {
    class <- append(class, c("template"))
  }
  if (any(grepl("CRISPRi", description))) {
    class <- append(class, c("CRISPRi"))
  }
  if (any(grepl("CRISPRa", description))) {
    class <- append(class, c("CRISPRa"))
  }
  return(paste(unique(class), collapse=","))
}

## Data processing

projects_classified <- projects %>%
  rowwise() %>%
  mutate(class = find_class(DESCRIPTION))

nrow(projects_classified) # 6599
table(projects_classified$class)

projects_classified_formatted <- projects_classified %>%
  # Reduce the number of different classes
  mutate(class = ifelse(class == "", "unclassified", class)) %>%
  mutate(class = ifelse(grepl("RNAseq", class), "RNAseq", class)) %>% # If RNAseq appears, consider it RNAseq
  mutate(class = ifelse(grepl("singlecell", class), "singlecell", class)) %>% # If singlecell appears, consider it singlecell
  mutate(class = ifelse(class == "screening,KO,KI", "screening", class)) %>%
  mutate(class = ifelse(class == "screening,KO,KI", "screening", class)) %>%
  mutate(class = ifelse(class == "KO,KI", "targeted", class)) %>%
  mutate(class = ifelse(class == "KO,template", "targeted", class)) %>%
  mutate(class = ifelse(class == "targeted,KI", "targeted", class)) %>%
  mutate(class = ifelse(class == "KI,template", "targeted", class)) %>%
  mutate(class = ifelse(class == "KO,targeted", "targeted", class)) %>%
  mutate(class = ifelse(class == "targeted,template", "targeted", class)) %>%
  mutate(class = ifelse(class == "KO,targeted,template", "targeted", class)) %>%
  mutate(class = ifelse(class == "KO,KI,template", "targeted", class)) %>%
  mutate(class = ifelse(grepl("screening", class), "screening", class)) %>% # If they contain screening, give it priority and consider them screening
  # Add potential capability of analysing with crisprseq
  mutate(analysis = ifelse(class == "CRISPRa", "screening",
                    ifelse(class == "CRISPRi", "screening",
                    ifelse(class == "targeted,CRISPRa", "screening",
                    ifelse(class == "targeted,CRISPRi", "screening",
                    ifelse(class == "targeted,CRISPRi,CRISPRa", "screening",
                    ifelse(class == "CRISPRi,CRISPRa", "screening",
                    ifelse(class == "KI", "targeted",
                    ifelse(class == "KO", "targeted",
                    ifelse(class == "template", "targeted",
                    ifelse(class == "KO,CRISPRi", "screening",
                    ifelse(class == "KO,CRISPRi,CRISPRa", "screening",
                    ifelse(class == "KO,targeted,CRISPRi", "screening",
                    ifelse(class == "RNAseq", "rna",
                    ifelse(class == "screening", "screening",
                    ifelse(class == "singlecell", "singlecell",
                    ifelse(class == "targeted", "targeted",
                    ifelse(class == "unclassified", "unclassified", class))))))))))))))))))
  
### Summarise by class

projects_classified_summary <- projects_classified_formatted %>%
  group_by(class) %>% 
  summarise(class_count = n()) %>%
  mutate(class_percentage = class_count / sum(class_count) * 100)

### Summarise by analysis

projects_classified_summary_analysis <- projects_classified_formatted %>%
  group_by(analysis) %>% 
  summarise(analysis_count = n()) %>%
  mutate(analysis_percentage = analysis_count / sum(analysis_count) * 100)
projects_classified_summary_analysis
# rna                    24.3%
# screening              28.6%
# singlecell             1.99%
# targeted               11.0%
# unclassified           34.0%

projects_classified_summary_analysis_dna <- projects_classified_formatted %>% 
  # Exclude RNA data and unclassified projects from the calculation
  filter(analysis != "rna") %>% 
  filter(analysis != "unclassified") %>%
  group_by(analysis) %>% 
  summarise(analysis_count = n()) %>%
  mutate(analysis_percentage = analysis_count / sum(analysis_count) * 100)
projects_classified_summary_analysis_dna
# screening              68.8%
# singlecell             4.77%
# targeted               26.5%

## Plot

analysis_percentage_plot <- ggplot(projects_classified_summary_analysis, aes(fill=analysis, x="", y=analysis_percentage)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  labs(x="", y="Percentage of datasets") +
  scale_fill_manual(name = "Type of analysis", values = c("#fc5858", "#8aeba4", "#e3db9d", "#4abd68", "#a486bf"))
analysis_percentage_plot
ggsave(plot=analysis_percentage_plot, filename = "analysis_percentage.png", device="png",width = 20, height = 15, units="cm")



analysis_percentage_plot_dna <- ggplot(projects_classified_summary_analysis_dna, aes(fill=analysis, x="", y=analysis_percentage)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  labs(x="", y="Percentage of datasets") +
  scale_fill_manual(name = "Type of analysis", values = c("#8aeba4", "#84b0d1", "#d582f5"))
analysis_percentage_plot_dna
ggsave(plot=analysis_percentage_plot_dna, filename = "analysis_percentage_dna-only.png", device="png",width = 20, height = 15, units="cm")
