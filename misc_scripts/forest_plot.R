library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
theme_set(theme_cowplot())


analysis_type <- "competing_risks"
#analysis_type <- commandArgs(trailingOnly=TRUE)
conc_limit <- 0.75

res <- readRDS(paste0("fuel_for_plot/forest.", analysis_type, ".RDS"))
colnames(res) <- c("coef", "exp_coef", "se", "z", "p", "author", "model")

if(analysis_type != "age_covar_AGE"){
  plain_res <- readRDS(paste0("fuel_for_plot/forest.age_covar.RDS"))
  plain_res <- plain_res[plain_res$model == "unchanged",]
} else {
  plain_res <- res[res$model == "unchanged",]
}

colnames(plain_res) <- c("coef", "exp_coef", "se", "z", "p", "author", "model")


source("common.R")

#res <- res[res$model %in% model_names$old,]

if(analysis_type == "adjustment" | analysis_type == "disease_labels"){
  python_acc <- readRDS(paste0("fuel_for_plot/python.", analysis_type, ".RDS"))
  
  res <- res[res$author %in% python_acc$author[python_acc$concordance > conc_limit], ]
  plain_res <- plain_res[plain_res$author %in% python_acc$author[python_acc$concordance > conc_limit], ]
}

res$author <- swap_names(res$author, disease_names)
res$model <- swap_names(res$model, model_names)

res <- res[!(res$model %in% c("Left Censor at EHR", "Censor at EHR", "Incl 1st Feat.", "Incl. 2nd Feat.", "Incl. 3rd Feat.")),]

plain_res$author <- swap_names(plain_res$author, disease_names)

############################################

for_meta_plot <- list()
just_plain <- list()

for(i in 1:length(unique(res$author))){
  sub_res <- res[res$author == unique(res$author)[i],]
  plain_sub_res <- plain_res[plain_res$author == unique(res$author)[i],]
  
  if(analysis_type != "age_covar_AGE"){
    if(plain_sub_res$coef[1] < 0){sub_res$coef <- sub_res$coef * -1} #swap so prs inflicts risk
  }
  graph_data <- as.matrix(data.frame("lower" = c(NA, sub_res$coef - sub_res$se),
                                     "mean" = c(NA, sub_res$coef),
                                     "upper" = c(NA, sub_res$coef + sub_res$se)))
  label_data <- cbind(c("Model", sub_res$model),
                      c("Coefficient", signif(sub_res$coef, 4)),
                      c("P-Value", signif(sub_res$p, 3)))
  summary_bool <- c(TRUE, rep(FALSE, nrow(sub_res)))
  
  # png(paste0("out_plots/forest/forest.", analysis_type, ".", gsub(" ", "_", unique(res$author)[i]), ".png"),
  #     width = 450, height = 180)
  # forestplot(label_data, graph_data, is.summary = summary_bool, zero = NA)
  # dev.off()
  
  
  
  for_meta_plot[[i]] <- data.frame("model" = sub_res$model[sub_res$model != "Plain"],
                                   "just_coef" = sub_res$coef[sub_res$model != "Plain"],
                                   "coef_diff" = sub_res$coef[sub_res$model != "Plain"] - 
                                     plain_sub_res$coef,
                                   "p" = sub_res$p[sub_res$model != "Plain"],
                                    stringsAsFactors = F)
  for_meta_plot[[i]]$perc_diff <- for_meta_plot[[i]]$coef_diff/plain_sub_res$coef
  for_meta_plot[[i]]$disease <- unique(res$author)[i]
  
  if(analysis_type == "age_covar"){
  just_plain[[i]] <- data.frame("model" = sub_res$model[sub_res$model == "Plain"],
                                   "just_coef" = sub_res$coef[sub_res$model == "Plain"],
                                   "coef_diff" = sub_res$coef[sub_res$model == "Plain"] - plain_sub_res$coef,
                                   "p" = sub_res$p[sub_res$model == "Plain"], "disease" = unique(res$author)[i],
                                   stringsAsFactors = F)
  }
}


#########################################################
#should do meta plot
#look at perc diff and coef diff from unchaged to others as a boxplot?

plot_df <- do.call("rbind", for_meta_plot)

to_save <- plot_df[,c(1,6,2,3,5,4)]
colnames(to_save) <- c("Model", "Disease", "Coef.", "Coef. Diff. from Plain", "Coef. % Diff. from Plain", "P")
to_save[,3] <- signif(to_save[,3], 3)
to_save[,4] <- signif(to_save[,4], 3)
to_save[,5] <- signif(to_save[,5], 3)
to_save[,6] <- signif(to_save[,6], 3)
write.table(to_save,
            paste0("supp_tables/coef.", analysis_type, ".txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")

if(analysis_type == "age_covar"){
  just_plain <- do.call("rbind", just_plain)
  just_plain <- just_plain[,c(5,2,4)]
  just_plain[,2] <- signif(just_plain[,2], 3)
  just_plain[,3] <- signif(just_plain[,3], 3)
  write.table(just_plain, "supp_tables/coef.plain.txt", row.names = F, col.names = T, quote = F, sep = "\t")
}



if(analysis_type == "adjustment"){
  box_disease <- names(table(plot_df$disease)[table(plot_df$disease) == max(table(plot_df$disease))])
} else {
  box_disease <- unique(plot_df$disease)
}

uvals <- unlist(lapply(unique(plot_df$model), function(x) mean(plot_df$coef_diff[plot_df$model == x])))
plot_df$model <- factor(plot_df$model, levels = unique(plot_df$model)[order(uvals)])

#plot_df$p[plot_df$p < 1e-10] <- 1e-10
the_plot <- ggplot(plot_df[plot_df$disease %in% box_disease,], aes(model, coef_diff)) + geom_boxplot() + geom_point() +
  labs(x = "", y = "Coefficient Difference\nFrom Plain Model") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/forest.", analysis_type, ".meta_coef_diff.png"), the_plot, width = 4.5, height = 4.4)


the_plot <- ggplot(plot_df[plot_df$disease %in% box_disease,], aes(model, -log10(p))) + geom_boxplot() + geom_point() +
  labs(x = "", y = "-log10(P-value) of Coef.") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
plot(the_plot)


the_plot <- ggplot(plot_df[plot_df$disease %in% box_disease,], aes(model, perc_diff*100)) + geom_boxplot() + geom_point() +
  labs(x = "", y = "Percent Difference in\nCoefficient From Plain Model") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/forest.", analysis_type, ".meta_perc_diff.png"), the_plot, width = 4.5, height = 4.4)




the_plot <- ggplot(plot_df[plot_df$disease %in% box_disease,], aes(perc_diff*100, disease, color = model)) + geom_point() +
  labs(x = "Coefficient % Difference\nFrom Plain Model", y = "", color = "Analysis") +
  geom_hline(yintercept = 1:length(unique(plot_df$disease))+0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/forest.", analysis_type, ".meta_all_perc_diff.png"), the_plot, width = 6.3, height = 4.4)



#ggplot(plot_df[plot_df$disease %in% box_disease,], aes(coef_diff, disease)) + geom_point()
