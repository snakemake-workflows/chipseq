log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")

homer_data <- read_tsv(snakemake@input[[1]])
homer_data <- homer_data %>% gather(`exon`, `Intergenic`, `intron`, `promoter-TSS`, `TTS`, key="sequence_element", value="counts")

peaks_sum <- ggplot(homer_data, aes(x = counts, y = sample, fill = sequence_element)) +
  geom_bar(position="fill", stat="Identity") +
  theme_minimal() +
  labs(x="", y="Peak count") +
  theme(legend.position = "right") +
  guides(fill=guide_legend("sequence element")) +
  ggtitle("Peak to feature proportion")

ggsave(snakemake@output[[1]], peaks_sum)

