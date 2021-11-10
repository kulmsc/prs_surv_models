library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
theme_set(theme_cowplot())


#analysis_type <- "competing_risks"
analysis_type <- commandArgs(trailingOnly=TRUE)
res <- readRDS(paste0("fuel_for_plot/conc.", analysis_type, ".RDS"))


source("common.R")


######################################################
all_res <- list()
for(i in 1:length(res)){
  plot_df <- data.frame("type" = rep(swap_names(rownames(res[[i]]), model_names), 3),
                        "val" = c(res[[i]]$base_val, res[[i]]$score_val, res[[i]]$diff_val),
                        "err" = c(res[[i]]$base_var, res[[i]]$score_var, res[[i]]$diff_var),
                        "model" = rep(c("Base", "Score", "Imp."), each = nrow(res[[i]])))
  
  
  the_plot <- ggplot(plot_df[plot_df$model != "Imp.",], aes(type, val, color = model)) + geom_point() +
    geom_errorbar(aes(ymin = val - 2*err, ymax = val + 2*err, width = 0)) +
    theme(axis.text=element_text(size=10)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    labs(x = "Analysis", y = "Concordance", color = "Model")
  plot(the_plot)
  ggsave(paste0("out_plots/conc/conc.", analysis_type, ".",
                gsub(" ", "_", swap_names(names(res)[i], disease_names)), ".png"),
         the_plot, width = 6, height = 4)
  
  #all_res[[i]] <- plot_df[plot_df$model == "Imp.",]
  all_res[[i]] <- plot_df
  all_res[[i]]$disease <- swap_names(names(res)[i], disease_names)

}


all_df <- do.call("rbind", all_res)
plot_df <- all_df[all_df$model == "Imp.",]

the_plot <- ggplot(all_df[all_df$model != "Imp.",], aes(type, val, color = model)) + geom_boxplot() +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Analysis", y = "Concordance", color = "Model")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".abso.png"), the_plot, width = 6.8, height = 5)


the_plot <- ggplot(plot_df, aes(type, val)) + geom_boxplot() + geom_point() +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Analysis", y = "Concordance Improvement")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".imp.png"), the_plot, width = 6.8, height = 5)


###########################################################################
#Minus Plain


if(!("Plain" %in% all_df$type)){
  plain_add_on <- readRDS("fuel_for_plot/conc.plain.1.RDS")
  plot_df <- all_df[all_df$model == "Score",]
  plot_df$val <- plot_df$val - custom_rep(plain_add_on$val[plain_add_on$model == "Score"], plot_df$disease)
} else {
  plot_df <- all_df[all_df$model == "Score" & all_df$type != "Plain",]
  plot_df$val <- plot_df$val - rep(all_df$val[all_df$model == "Score" & all_df$type == "Plain"],
                                   each =  length(unique(all_df$type))-1)
}

plot_df$disease <- factor(plot_df$disease, levels = plot_df$disease[plot_df$type == plot_df$type[1]][
    order(plot_df$val[plot_df$type == plot_df$type[1]])])


the_plot <- ggplot(plot_df, aes(type, val)) + geom_boxplot() + geom_point() +
  labs(x = "Analysis", y = "Score Concordance:\nAnalysis Minus Plain", color = "Analysis") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".box_minus_plain.png"),
       the_plot, width = 6.8, height = 5)


the_plot <- ggplot(plot_df, aes(val, disease, color = type)) + geom_point() +
  labs(x = "Score Concordance:\nAnalysis Minus Plain", y = "", color = "Analysis")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".disease_minus_plain.png"),
       the_plot, width = 6.8, height = 5)

