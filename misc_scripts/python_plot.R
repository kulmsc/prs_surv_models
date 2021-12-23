library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
theme_set(theme_cowplot())

#analysis_type <- "disease_labels"
analysis_type <- "adjustment"
#analysis_type <- commandArgs(trailingOnly=TRUE)

res <- readRDS(paste0("fuel_for_plot/relative_risk.", analysis_type, ".RDS"))

source("common.R")


#################################################################

plot_df <- res[["python"]]
plot_df$author <- swap_names(plot_df$author, disease_names)

colnames(plot_df) <- c("conc", "std", "author")
plot_df$author <- factor(plot_df$author, levels = plot_df$author[order(plot_df$conc)])

the_plot <- ggplot(plot_df, aes(conc, author)) + geom_point() +
  labs(x = "Concordance", y = "") +
  geom_errorbar(aes(xmin = conc - std, xmax = conc + std, width = 0)) +
  geom_hline(yintercept = (1:length(unique(plot_df$author)))+0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("out_plots/meta/python/", analysis_type, ".png"), the_plot, width = 4.5, height = 5)
