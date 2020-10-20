log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")

data <- lapply(snakemake@input, read.table, header=F, stringsAsFactors = F)
frip_scores <- tibble()
for (i in 1:length(data)) {
  frip_scores <- rbind(frip_scores, data[[i]])
}
names(frip_scores) <- c("sample_control", "frip")

frip <- ggplot(frip_scores, aes(x = sample_control, y = frip, fill = sample_control)) +
  geom_bar(stat="Identity", color="black") +
  theme_minimal() +
  labs(x="", y="FRiP score") +
  theme(legend.position = "right") +
  guides(fill=guide_legend("samples with controls")) +
  ggtitle("FRiP score")

ggsave(snakemake@output[[1]], frip)
