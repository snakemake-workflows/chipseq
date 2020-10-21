log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")

data <- lapply(snakemake@input, read.table, header=F, stringsAsFactors = F)
counts <- tibble()
for (i in 1:length(data)) {
  counts <- rbind(counts, data[[i]])
}
names(counts) <- c("sample_control", "count")

peaks_counts <- ggplot(counts, aes(x = count, y = sample_control, fill=sample_control)) +
  geom_bar(stat="Identity", color="black") +
  theme_minimal() +
  labs(x="Peak count", y="") +
  theme(legend.position = "right") +
  guides(fill=guide_legend("samples with controls")) +
  ggtitle("Total peak count")

ggsave(snakemake@output[[1]], peaks_counts)
