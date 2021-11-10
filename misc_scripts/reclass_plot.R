library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
theme_set(theme_cowplot())

#analysis_type <- "age_covar"
analysis_type <- commandArgs(trailingOnly=TRUE)

res <- readRDS(paste0("fuel_for_plot/reclass.", analysis_type, ".RDS"))


source("common.R")
cut_names <- data.frame("old" = c("cut_0.9", "cut_0.95", "cut_0.99"),
                        "new" = c("90%", "95%", "99%"), stringsAsFactors = F)
other_model_names <- data.frame("old" = c("base", "diff", "score"),
                                "new" = c("Base", "Base - Score", "Score"), stringsAsFactors = F)


##############################################
all_res <- list()
for(i in 1:length(res)){
  
  plot_df <- melt(res[[i]], id.vars = c("model", "type"))
  plot_df$type <- swap_names(plot_df$type, model_names)
  plot_df$variable <- swap_names(as.character(plot_df$variable), cut_names)
  plot_df$model <- swap_names(plot_df$model, other_model_names)
  plot_df$value <- 1 - plot_df$value
  
  the_plot <- ggplot(plot_df, aes(value, type, color = model)) + geom_point() +
    facet_grid( ~ variable, scales = "free") +
    labs(x = "Reclassification Rate", color = "Model\nPred.", y = "Analysis")
  plot(the_plot)
  ggsave(paste0("out_plots/reclass/reclass.", analysis_type, ".",
                gsub(" ", "_", swap_names(names(res)[i], disease_names)), ".png"),
         the_plot, width = 9, height = 5)
  
  all_res[[i]] <- plot_df
  all_res[[i]]$disease <- swap_names(names(res)[i], disease_names)
}



plot_df <- do.call("rbind", all_res)

the_plot <- ggplot(plot_df[plot_df$variable == "95%",], aes(type, value, color = model)) +
  geom_boxplot() + 
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Analysis", y = "Reclassification Rate", color = "Model\nPred.")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/reclass.", analysis_type, ".box.png"), the_plot, width = 6.8, height = 5)



start_levels <- plot_df$disease[
  plot_df$variable == "95%" & plot_df$model == "Base - Score" & plot_df$type == plot_df$type[1]][order(
    plot_df$value[plot_df$variable == "95%" & plot_df$model == "Base - Score" & plot_df$type == plot_df$type[1]])]
u_disease <- unique(plot_df$disease)
start_levels <- c(start_levels, u_disease[!(u_disease %in% start_levels)])
plot_df$disease <- factor(plot_df$disease, levels = start_levels)



the_plot <- ggplot(plot_df[plot_df$variable == "95%",], aes(value, disease, color = type)) + geom_point() +
  facet_grid( ~ model, scales = "free") + 
  labs(x = "Reclassification Rate", y = "", color = "Analysis")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/reclass.", analysis_type, ".dot_all_rates.png"), the_plot, width = 11, height = 5)



the_plot <- ggplot(plot_df[plot_df$variable == "95%" & plot_df$model == "Score",],
       aes(value, disease, color = type)) + geom_point() +
  labs(x = "Reclassification Rate", y = "", color = "Analysis")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/reclass.", analysis_type, ".dot_95.png"), the_plot, width = 6.8, height = 5)

