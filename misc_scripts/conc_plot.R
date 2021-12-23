library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
library(stringr)
theme_set(theme_cowplot())


analysis_type <- "adjustment"
#analysis_type <- commandArgs(trailingOnly=TRUE)
res <- readRDS(paste0("fuel_for_plot/conc.", analysis_type, ".RDS"))
conc_limit <- 0.75

source("common.R")

if(analysis_type == "adjustment" | analysis_type == "disease_labels"){
  source("check_acc.R")
}


######################################################
all_res <- list()


for(i in 1:length(res)){

  res[[i]] <- res[[i]][!(rownames(res[[i]]) %in% c("left_ehr_censor", "ehr", "just1_coef", "just2_coef", "just_3_coef")),]
  
  
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
         the_plot, width = 5, height = 4)
  
  #all_res[[i]] <- plot_df[plot_df$model == "Imp.",]
  all_res[[i]] <- plot_df
  all_res[[i]]$disease <- swap_names(names(res)[i], disease_names)

}


################################################################################################333

all_df <- do.call("rbind", all_res)


if(analysis_type == "competing_risks" | analysis_type == "adjustment"){
  box_disease <- names(table(all_df$disease)[table(all_df$disease) == max(table(all_df$disease))])
} else {
  box_disease <- unique(all_df$disease)
}

if(!("Plain" %in% all_df$type)){
  plain_add_on <- readRDS("fuel_for_plot/conc.plain.1.RDS")
  plain_add_on <- plain_add_on[plain_add_on$disease %in% all_df$disease,]
  all_df <- rbind(all_df, plain_add_on) 
}




plot_df <- all_df[all_df$model == "Imp.",]

the_plot <- ggplot(all_df[all_df$model != "Imp." & all_df$disease %in% box_disease,], aes(type, val, color = model)) +
  geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Analysis", y = "Concordance", color = "Model")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".abso.png"), the_plot, width = 6.8, height = 5)


the_plot <- ggplot(all_df[all_df$model == "Score" & all_df$disease %in% box_disease,], aes(type, val)) +
  geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Analysis", y = "Concordance")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".abso.png"), the_plot, width = 4.8, height = 3.8)
save_abs_stat_val <- all_df[all_df$model == "Score",]


the_plot <- ggplot(plot_df[plot_df$disease %in% box_disease,], aes(type, val)) + geom_boxplot() + geom_point() +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Analysis", y = "Concordance\nImprovement")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".imp.png"), the_plot, width = 4.8, height = 3.8)
save_imp_stat_val <- data.frame(plot_df)



new_plot <- plot_df[plot_df$type != "Plain",]
new_plot$val <- new_plot$val - custom_rep(plot_df$val[plot_df$type == "Plain"], new_plot$disease)
new_plot$perc_diff <- new_plot$val/custom_rep(plot_df$val[plot_df$type == "Plain"], new_plot$disease)
#new_plot$type <- factor(new_plot$type, levels = c("No Censoring", "Delay Entry by Age Time", "Delay Entry by Study Time"))
#new_plot$disease <- factor(new_plot$disease, levels = new_plot$disease[new_plot$type == "Delay Entry by Study Time"][
#  order(new_plot$val[new_plot$type == "Delay Entry by Study Time"])])
save_diff_imp_stat_val <- data.frame(new_plot)
if(analysis_type == "censoring"){
  new_plot$fresh_type <- factor(gsub("\n", " ", as.character(new_plot$type)),
                              levels = c("No Censoring", "Delay Entry by Age Time", "Delay Entry by Study Time"))
} else {
  new_plot$fresh_type <- new_plot$type
}
the_plot <- ggplot(new_plot, aes(val, disease, color = fresh_type)) + geom_point() +
  labs(x = "Score Concordance Imp.:\nAnalysis Minus Plain", y = "", color = "Analysis") +
  geom_hline(yintercept = (1:length(unique(plot_df$disease)))+0.5, color = "grey80") 
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".disease_IMP_minus_plain.png"),
       the_plot, width = 7, height = 5)



###########################################################################
#Minus Plain

all_df <- do.call("rbind", all_res)

if(analysis_type == "age_covar"){
  plain_supp <- all_df[all_df$type == "Plain" & all_df$model == "Score",]
  plain_supp <- plain_supp[,c(5,2,3)]
  plain_supp$imp <- plain_supp$val - all_df[all_df$type == "Plain" & all_df$model == "Base",]$val
  plain_supp[,2] <- signif(plain_supp[,2])
  plain_supp[,3] <- signif(plain_supp[,3])
  plain_supp[,4] <- signif(plain_supp[,4])
  write.table(plain_supp, "supp_tables/conc.plain.txt", row.names = F, col.names = T, quote = F, sep = "\t")
}


if(!("Plain" %in% all_df$type)){
   plain_add_on <- readRDS("fuel_for_plot/conc.plain.1.RDS")
   plain_add_on <- plain_add_on[plain_add_on$disease %in% all_df$disease,]
   plot_df <- all_df[all_df$model == "Score",]
   plot_df$val <- plot_df$val - custom_rep(plain_add_on$val[plain_add_on$model == "Score"], plot_df$disease)
   plot_df$perc_diff <- plot_df$val/custom_rep(plain_add_on$val[plain_add_on$model == "Score"], plot_df$disease)
} else {
  plot_df <- all_df[all_df$model == "Score" & all_df$type != "Plain",]
  plot_df$val <- plot_df$val - custom_rep(all_df$val[all_df$model == "Score" & all_df$type == "Plain"],
                                          plot_df$disease)
  #plot_df$val <- plot_df$val - rep(all_df$val[all_df$model == "Score" & all_df$type == "Plain"],
  #                                 each =  length(unique(all_df$type))-1)
  plot_df$perc_diff <- plot_df$val/custom_rep(all_df$val[all_df$model == "Score" & all_df$type == "Plain"],
                                              plot_df$disease)
}

plot_df$disease <- factor(plot_df$disease, levels = plot_df$disease[plot_df$type == plot_df$type[1]][
    order(plot_df$val[plot_df$type == plot_df$type[1]])])
plot_df$type <- as.character(plot_df$type)
#plot_df$type <- factor(plot_df$type, levels = c("No Censoring", "Delay Entry by Age Time", "Delay Entry by Study Time"))


the_plot <- ggplot(plot_df[plot_df$disease %in% box_disease,], aes(type, val)) + geom_boxplot() + geom_point() +
  labs(x = "Analysis", y = "Score Concordance:\nAnalysis Minus Plain", color = "Analysis") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".box_minus_plain.png"),
       the_plot, width = 4.8, height = 3.8)


save_diff_abs_stat_val <- data.frame(plot_df)
if(analysis_type == "censoring"){
  plot_df$type <- factor(gsub("\n", " ", plot_df$type), levels = c("No Censoring", "Delay Entry by Age Time", "Delay Entry by Study Time"))
}
the_plot <- ggplot(plot_df, aes(val, disease, color = type)) + geom_point() +
  labs(x = "Score Concordance:\nAnalysis Minus Plain", y = "", color = "Analysis") +
  geom_hline(yintercept = (1:length(unique(plot_df$disease)))+0.5, color = "grey80") 
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/conc.", analysis_type, ".disease_minus_plain.png"),
       the_plot, width = 7, height = 5)




########
#Create Supp

final_save <- join_save(save_abs_stat_val[,c("type", "disease", "val")],
                        save_imp_stat_val[,c("type", "disease", "val")],
                        save_diff_abs_stat_val[,c("type", "disease", "perc_diff")],
                        save_diff_imp_stat_val[,c("type", "disease", "perc_diff")], "Conc.")
final_save$Analysis <- gsub("\n", " ", as.character(final_save$Analysis))
#saveRDS(join_save, "supp_tables/conc.RDS")
write.table(final_save, paste0("supp_tables/conc.", analysis_type, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
