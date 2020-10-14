log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")

data <- lapply(snakemake@input, read.table, header=F, stringsAsFactors = F)
data2 <- tibble()
for (i in 1:length(data)) {
  data2 <- rbind(data2, data[[i]])
}
names(data2) <- c("sample_control", "frip")

frip <- ggplot(data2, aes(x = sample_control, y = frip, fill = sample_control)) +
  geom_bar(stat="Identity", color="black") +
  theme_minimal() +
  labs(x="", y="FRiP score") +
  theme(legend.position = "none") +
  ggtitle("FRiP score")

ggsave(snakemake@output[[1]], frip)
