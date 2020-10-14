#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")
#
#library("tidyverse")
#
#data <- lapply(snakemake@input, read.table, header=T, stringsAsFactors = F)
#data2 <- tibble()
#for (i in 1:length(data)) {
#  data2 <- rbind(data2, data[[i]])
#}
#names(data2) <-  c("sample_control", "exon", "intergenic", "intron", "promoterTSS", "TTS")
#
#data2 <- data2 %>%
#   mutate(countsum=exon+intergenic+intron+promoterTSS+TTS)
#
#data2 <- data2 %>%
#   #mutate(exon=100*exon/countsum,
#          intergenic=100*intergenic/countsum,
#          intron=100*intron/countsum,
#          promoterTSS=100*promoterTSS/countsum,
#          TTS=100*TTS/countsum) %>%
#  select(-countsum)
#
#keys <- data2 %>%
#  select(sample_control)
#vals <- data2 %>%
#  select(-sample_control)
#forfill <- list(rep(colnames(data2[-1]),3))
#
#data2 <- data2 %>%
#  gather(key=sample_control, value=vals)
#
#peaksum <- ggplot(data2, aes(x = vals, y =keys, fill = forfill)) +
#  geom_bar(position="fill", stat="Identity") +
#  theme_minimal() +
#  labs(x="", y="Peak count") +
#  geom_text(aes(label = count), position = position_fill(vjust = 0.5)) +
#  theme(legend.position = "none") +
#  ggtitle("Peak to feature proportion")
#
#ggsave(snakemake@output[[1]], peaksum)
#
