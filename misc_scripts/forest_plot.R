library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
theme_set(theme_cowplot())

#analysis_type <- "age_covar"
analysis_type <- commandArgs(trailingOnly=TRUE)

res <- readRDS(paste0("fuel_for_plot/forest.", analysis_type, ".RDS"))
colnames(res) <- c("coef", "exp_coef", "se", "z", "p", "author", "model")

plain_res <- readRDS(paste0("fuel_for_plot/forest.age_covar.RDS"))
colnames(plain_res) <- c("coef", "exp_coef", "se", "z", "p", "author", "model")


source("common.R")

res$author <- swap_names(res$author, disease_names)
res$model <- swap_names(res$model, model_names)

############################################

for_meta_plot <- list()

for(i in 1:length(unique(res$author))){
  sub_res <- res[res$author == unique(res$author)[i],]
  plain_sub_res <- plain_res[plain_res$author == unique(plain_res$author)[i],]
  
  if(sub_res$coef[1] < 0){sub_res$coef <- sub_res$coef * -1}
  graph_data <- as.matrix(data.frame("lower" = c(NA, sub_res$coef - sub_res$se),
                                     "mean" = c(NA, sub_res$coef),
                                     "upper" = c(NA, sub_res$coef + sub_res$se)))
  label_data <- cbind(c("Model", sub_res$model),
                      c("Coefficient", signif(sub_res$coef, 4)),
                      c("P-Value", signif(sub_res$p, 3)))
  summary_bool <- c(TRUE, rep(FALSE, nrow(sub_res)))
  
  png(paste0("out_plots/forest/forest.", analysis_type, ".", gsub(" ", "_", unique(res$author)[i]), ".png"),
      width = 450, height = 180)
  forestplot(label_data, graph_data, is.summary = summary_bool, zero = NA)
  dev.off()
  
  
  
  for_meta_plot[[i]] <- data.frame("model" = sub_res$model[sub_res$model != "Plain"],
                                   "coef_diff" = sub_res$coef[sub_res$model != "Plain"] - 
                                     plain_sub_res$coef[plain_sub_res$model == "unchanged"],
                                    stringsAsFactors = F)
  for_meta_plot[[i]]$perc_diff <- for_meta_plot[[i]]$coef_diff/plain_sub_res$coef[plain_sub_res$model == "unchanged"]
  for_meta_plot[[i]]$disease <- unique(res$author)[i]
}


#########################################################
#should do meta plot
#look at perc diff and coef diff from unchaged to others as a boxplot?

plot_df <- do.call("rbind", for_meta_plot)

the_plot <- ggplot(plot_df, aes(model, coef_diff)) + geom_boxplot() + geom_point() +
  labs(x = "", y = "Coefficient Difference\nFrom Plain Model") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/forest.", analysis_type, ".meta_coef_diff.png"), the_plot, width = 6.3, height = 4.4)



the_plot <- ggplot(plot_df, aes(model, perc_diff)) + geom_boxplot() + geom_point() +
  labs(x = "", y = "Percent Difference in Coefficient\nFrom Plain Model") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/forest.", analysis_type, ".meta_perc_diff.png"), the_plot, width = 6.3, height = 4.4)
