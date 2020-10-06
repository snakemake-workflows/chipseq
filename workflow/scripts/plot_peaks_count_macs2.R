log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")

data <- lapply(snakemake@input, read.table, header=F, stringsAsFactors = F)
data2 <- tibble()
names(data2) <- c("sample_control", "count")
for (i in 1:length(data)) {
  data2 <- rbind(data2, data[[i]])
}
names(data2) <- c("sample_control", "count")

peaks_counts <- ggplot(data2, aes(x = count, y = sample_control, fill=sample_control)) +
  geom_bar(stat="Identity", color="black") +
  theme_minimal() +
  labs(x="Peak count", y="") +
  theme(legend.position = "none") +
  ggtitle("Total peak count")

ggsave(snakemake@output[[1]], peaks_counts)
