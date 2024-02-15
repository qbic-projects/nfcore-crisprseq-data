# nf-core/crisprseq benchmarking

### Import all R libraries
library(dplyr)
library(tidyr)

## Benchmarking of resources

"Calculate average runtime and CPU hours of nf-core/crisprseq vs CRISPR-Analytics workflows."

### Calculate total runtime and CPU hours

data_directory = "benchmarking_results/"
crisprseq_totals <- read.table(file = paste0(data_directory, "crisprseq_full/totals_crisprseq.tsv"), sep='\t', header=TRUE)
crispra_totals <- read.table(file = paste0(data_directory, "crispra_full/totals_crispra.tsv"), sep='\t', header=TRUE)

totals <- rbind(crisprseq_totals, crispra_totals)

totals_formatted <- totals %>%
  separate(time, c("hours", "minutes", "seconds"), "h|m|s", remove=FALSE) %>%
  mutate(total_minutes = as.numeric(hours) * 60 + as.numeric(minutes) + as.numeric(seconds) / 60) %>%
  #mutate(minutes = as.numeric(strsplit(time, 'h')[1]) * 60 + as.numeric(strsplit(time, 'h|m')[[1]][2]) + as.numeric(strsplit(time, 'h|m|s')[[1]][3]) / 60) %>%
  mutate(run = as.factor(run)) %>%
  mutate(subrun = as.factor(subrun)) %>%
  mutate(workflow = as.factor(workflow))
  
totals_summary <- totals_formatted %>%
  # Sum all crispra subruns
  group_by(workflow, run) %>%
  mutate(sum_cpu_subruns = sum(CPU_hours)) %>%
  mutate(sum_time_subruns = sum(total_minutes)) %>%
  ungroup() %>%
  # Collapse
  select(workflow, run, sum_cpu_subruns, sum_time_subruns) %>%
  distinct(workflow, run, .keep_all = TRUE) %>%
  # Calculate the average of all runs
  group_by(workflow) %>%
  mutate(avg_minutes = mean(sum_time_subruns)) %>%
  mutate(avg_cpu = mean(sum_cpu_subruns)) %>%
  ungroup() %>%
  # Collapse
  select(workflow, avg_minutes, avg_cpu) %>%
  distinct(workflow, .keep_all = TRUE)
  
totals_summary %>%
  mutate(hours = as.integer(avg_minutes / 60)) %>%
  mutate(minutes = as.integer((avg_minutes / 60)%%1 * 60)) %>%
  mutate(seconds = as.integer(((avg_minutes / 60)%%1 * 60)%%1 * 60))

# crisprseq 14h 12m 45s 222CPUh
# crispra   31h 59m 8s  1919CPUh

