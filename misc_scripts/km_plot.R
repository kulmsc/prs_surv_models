library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
theme_set(theme_cowplot())


#analysis_type <- "disease_labels"
analysis_type <- commandArgs(trailingOnly=TRUE)
res <- readRDS(paste0("fuel_for_plot/km.", analysis_type, ".RDS"))
conc_limit <- 0.75


source("common.R")

if(analysis_type == "adjustment" | analysis_type == "disease_labels"){
  source("check_acc.R")
}


############################################

if(!("unchanged" %in% res[[1]]$type)){
  plain_data <- readRDS("fuel_for_plot/km.plain.1.RDS")
  full_plain_data <- readRDS("fuel_for_plot/km.plain.2.RDS")
  
  plain_data <- plain_data[plain_data$model == "Score",]
  full_plain_data <- full_plain_data[full_plain_data$model == "Score",]
}



#all_over_time_df <- list()
all_imp_df <- list() #score vals minus base vals - last time
all_diff_df <- list() #score vals minus plain vals - last time
all_diff_perc_df <- list()
all_score_df <- list() #score vals - last time
all_score_perc_df <- list() #score hi-risk perc change from inter - last time
all_imp_perc_df <- list() #imp hi-risk perc change from inter - last time



for(i in 1:length(res)){
  res[[i]]$author <- swap_names(res[[i]]$author, disease_names)
  res[[i]]$type <- swap_names(res[[i]]$type, model_names)
  res[[i]]$model <- stringr::str_to_title(res[[i]]$model)
  
  # if(res[[i]]$author[1] = "Psoriasis"){
  #   exit()
  # }

  res[[i]] <- res[[i]][!(res[[i]]$type %in% c("Left Censor at EHR", "Censor at EHR", "Incl 1st Feat.", "Incl. 2nd Feat.", "Incl. 3rd Feat.")),]


  if(analysis_type == "censoring"){
    res[[i]]$split_up <- "From Assessment (Right Censor)"
    res[[i]]$split_up[grepl("Left", res[[i]]$type)] <- "From Age 60 (Left Censor)"
    res[[i]]$split_up <- as.factor(res[[i]]$split_up)
  } else if(analysis_type == "competing_risks"){
    res[[i]]$split_up <- "From Assessment"
    res[[i]]$split_up[grepl("iCARE", res[[i]]$type)] <- "From Age 60"
    res[[i]]$split_up <- as.factor(res[[i]]$split_up)
    res[[i]]$time[res[[i]]$type == "iCARE"] <- res[[i]]$time[res[[i]]$type == "iCARE"] *365
    
    plain_data$split_up <- "From Assessment"
    full_plain_data$split_up <- "From Assessment"
  }
  

  #all_over_time_df[[i]] <- res[[i]][res[[i]]$type == "Plain",]
  
  plot_df <- res[[i]][res[[i]]$risk_group == "high" & res[[i]]$model == "Score",]

  
  if(!("Plain" %in% plot_df$type)){
    plot_df <- rbind(plot_df, full_plain_data[full_plain_data$author == plot_df$author[1] & 
                                              full_plain_data$risk_group == "high" & full_plain_data$model == "Score",])
  }

  if(analysis_type == "censoring" ){
    plot_df$time[plot_df$type == "Left Censor"] <- plot_df$time[plot_df$type == "Plain"]
    plot_df$time[plot_df$type == "Left Censor at EHR"] <- plot_df$time[plot_df$type == "Plain"]
    plot_df$time[plot_df$type == "Delay Entry by Age Time"] <- plot_df$time[plot_df$type == "Plain"]
    plot_df$time[plot_df$type == "Delay Entry by Study Time"] <- plot_df$time[plot_df$type == "Plain"]
  } 
  

  
  
  ##################################################################################3
  #score over time
  
  the_plot <- ggplot(plot_df, aes(time/365, val, color = type)) +
    geom_point() + geom_line() +
    labs(x = "Year From Assessment", y = "Cumulative Hazard", color = "Analysis")
  plot(the_plot)
  ggsave(paste0("out_plots/km/km.", analysis_type, ".over_time.", gsub(" ", "_", res[[i]]$author[1]), ".png"),
         the_plot, width = 6.3, height = 4)
  
  
 
  
  ##################################################################################3
  #score at top time

  plot_df <- get_top_time(res[[i]])
  
  if(!("Plain" %in% plot_df$type)){
    plot_df <- rbind(plot_df, plain_data[plain_data$author == plot_df$author[1],])
  }
  
  the_plot <- ggplot(plot_df, aes(type, val, color = risk_group, shape = model)) +
    geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
    labs(x = "Analysis", y = "Cumulative Hazard", color = "Risk Group", shape = "Model Type") +
    geom_errorbar(aes(ymin = val - sd, ymax = val + sd, width = 0), position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
    theme(axis.text=element_text(size=10)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    coord_cartesian(ylim = c(0, min(c(0.2, max(get_top_time(res[[i]])$val + get_top_time(res[[i]])$sd)*1.1))))
  plot(the_plot)
  ggsave(paste0("out_plots/km/km.", analysis_type, ".single_time.", gsub(" ", "_", res[[i]]$author[1]), ".png"),
         the_plot, width = 7.3, height = 5)
  all_score_df[[i]] <- plot_df
  
  
  
  ##################################################################################3
  #score perc diff over time
  
  plot_df <- res[[i]][res[[i]]$model == "Score",]
  base_plot_df <- res[[i]][res[[i]]$model == "Base",]
  
  if(!("Plain" %in% base_plot_df$type)){
    plot_df <- rbind(plot_df, full_plain_data[full_plain_data$author == plot_df$author[1] & full_plain_data$model == "Score",])
    base_plot_df <- rbind(base_plot_df, full_plain_data[full_plain_data$author == plot_df$author[1] & full_plain_data$model == "Base",])
  }
  if(analysis_type == "censoring" ){
    plot_df$time[plot_df$type == "Left Censor"] <- plot_df$time[plot_df$type == "Plain"]
    plot_df$time[plot_df$type == "Left Censor at EHR"] <- plot_df$time[plot_df$type == "Plain"]
    plot_df$time[plot_df$type == "Delay Entry by Age Time"] <- plot_df$time[plot_df$type == "Plain"]
    plot_df$time[plot_df$type == "Delay Entry by Study Time"] <- plot_df$time[plot_df$type == "Plain"]
  } 
  
  perc_plot_df <- scale_by_intermediate(plot_df)
  perc_plot_df <- get_top_time(perc_plot_df)
  
  the_plot <- ggplot(perc_plot_df[perc_plot_df$risk_group != "inter",], aes(type, val*100, color = risk_group)) + 
    geom_point() + geom_hline(yintercept = 0) +
    labs(x = "Analysis", y = "% Cumulative Hazard\nDiff. From Inter. Risk", color = "Analysis") +
    theme(axis.text=element_text(size=10)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  plot(the_plot)
  ggsave(paste0("out_plots/km/km.", analysis_type, ".perc_score.", gsub(" ", "_", res[[i]]$author[1]), ".png"),
         the_plot, width = 6.3, height = 4)
  all_score_perc_df[[i]] <- perc_plot_df[perc_plot_df$risk_group != "inter",]
  
  
  
  ##################################################################################3
  #imp over time
  

  plot_df$val <- plot_df$val - base_plot_df$val
  plot_df$group_by <- paste0(plot_df$risk_group, "_", plot_df$type)
  
  the_plot <- ggplot(plot_df[plot_df$risk_group == "high",], aes(time/365, val, color = type, group = group_by)) +
    geom_point() + geom_line() +
    labs(x = "Year From Assessment", y = "Imp. Cumulative Hazard", color = "Analysis")
  plot(the_plot)
  all_imp_df[[i]] <- get_top_time(plot_df)
  
  
  
  
  ##################################################################################3
  #imp perc at final time

  perc_plot_df <- scale_by_intermediate(plot_df)
  perc_plot_df <- get_top_time(perc_plot_df)
  
  the_plot <- ggplot(perc_plot_df[perc_plot_df$risk_group != "inter",], aes(type, val*100, color = risk_group)) + 
    geom_point() + geom_hline(yintercept = 0) +
    labs(x = "Analysis", y = "% Imp. Cumulative Hazard\nDiff. From Inter. Risk", color = "Analysis") +
    theme(axis.text=element_text(size=10)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  plot(the_plot)
  ggsave(paste0("out_plots/km/km.", analysis_type, ".perc_imp.", gsub(" ", "_", res[[i]]$author[1]), ".png"),
         the_plot, width = 6.3, height = 4)
  all_imp_perc_df[[i]] <- perc_plot_df[perc_plot_df$risk_group != "inter",]
  #perc_plot_df <- scale_by_intermediate(plot_df)
  
  

  
  
  
  #################################################################################
  
  
  
  temp_df <- get_top_time(res[[i]])
  imp_df <- temp_df[temp_df$model == "Score",]

  if("Plain" %in% imp_df$type){
    all_diff_df[[i]] <- imp_df[imp_df$type != "Plain",]
    all_diff_df[[i]]$val <- all_diff_df[[i]]$val - rep(imp_df$val[imp_df$type == "Plain"], nrow(all_diff_df[[i]])/3)
    
    all_diff_perc_df[[i]] <- imp_df[imp_df$type != "Plain",]
    all_diff_perc_df[[i]]$val <- (all_diff_df[[i]]$val - rep(imp_df$val[imp_df$type == "Plain"], nrow(all_diff_df[[i]])/3))/
      rep(imp_df$val[imp_df$type == "Plain"], nrow(all_diff_df[[i]])/3)
  } else {
    use_plain_data <- plain_data[plain_data$author %in% imp_df$author,]
    all_diff_df[[i]] <- imp_df
    all_diff_df[[i]]$val <- all_diff_df[[i]]$val - rep(use_plain_data$val[use_plain_data$author == imp_df$author[1]], nrow(all_diff_df[[i]])/3)
    
    all_diff_perc_df[[i]] <- imp_df
    all_diff_perc_df[[i]]$val <- (all_diff_df[[i]]$val - rep(use_plain_data$val[use_plain_data$author == imp_df$author[1]], nrow(all_diff_df[[i]])/3))/
      rep(use_plain_data$val[use_plain_data$author == imp_df$author[1]], nrow(all_diff_df[[i]])/3)
  }
  
}






#####################################################################################################
###
###
###
###
#######################################################################################################






imp_df <- do.call("rbind", all_imp_df)
imp_df$risk_group <- stringr::str_to_title(imp_df$risk_group)

if(analysis_type == "competing_risks" | analysis_type == "adjustment"){
  box_disease <- names(table(imp_df$author)[table(imp_df$author) == max(table(imp_df$author))])
} else {
  box_disease <- unique(imp_df$author)
}

imp_df <- imp_df[imp_df$author %in% box_disease,]

the_plot <- ggplot(imp_df, aes(type, val)) + geom_boxplot() + 
  labs(x = "Analysis", y = "Difference in Cum. Haz.:\nScore Model - Base Model", color = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_cartesian(ylim = c(-0.1, min(c(0.5, max(imp_df$val)*1.1))))
plot(the_plot)
#ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".imp_last_time.png"), the_plot, width = 6.8, height = 5)


############################################


#over_time_df <- do.call("rbind", all_over_time_df)

score_df <- do.call("rbind", all_score_df)
score_df$risk_group <- stringr::str_to_title(score_df$risk_group)

score_df <- score_df[score_df$author %in% box_disease,]

the_plot <- ggplot(score_df, aes(type, val, color = risk_group)) + geom_boxplot() +
  labs(x = "Analysis", y = "Score Model Cum. Haz.", color = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".score_last_time.png"), the_plot, width = 6.8, height = 5)


supp_df <- score_df[score_df$model == "Score",c(5,7,2,3,4)]
supp_df$val <- signif(supp_df$val, 3)
supp_df$sd <- signif(supp_df$sd, 3)
colnames(supp_df) <- c("Disease", "Analysis", "Risk Group", "Avg. Cum. Haz.", "SD Cum. Haz.")
write.table(supp_df, paste0("supp_tables/km.", analysis_type, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

#################################################


diff_df <- do.call("rbind", all_diff_df)
diff_df$risk_group <- stringr::str_to_title(diff_df$risk_group)

diff_df <- diff_df[diff_df$author %in% box_disease,]

the_plot <- ggplot(diff_df, aes(type, val, color = risk_group)) + geom_boxplot() +
  labs(x = "Analysis", y = "Difference in Cum. Haz.:\nScore Model - Plain Model", color = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(the_plot)
#ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".plaindiff_last_time.png"), the_plot, width = 6.8, height = 5)


#################################################

diff_df <- do.call("rbind", all_diff_perc_df)
diff_df$risk_group <- stringr::str_to_title(diff_df$risk_group)

diff_df <- diff_df[diff_df$author %in% box_disease & diff_df$risk_group == "High",]

the_plot <- ggplot(diff_df[diff_df$risk_group == "High",], aes(type, val*100)) + geom_boxplot() +
  labs(x = "Analysis", y = "% Diff. Score Cum. Haz.\nFrom Plain Cum. Haz.", color = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".perc_diff_plain_box.png"), the_plot, width = 6.8, height = 5)


diff_df$author <- factor(diff_df$author, levels = diff_df$author[diff_df$type == diff_df$type[1]][order(diff_df$val[diff_df$type == diff_df$type[1]])])
the_plot <- ggplot(diff_df[diff_df$author %in% box_disease & diff_df$risk_group == "High",], aes(val*100, author, color = type)) + geom_point() +
  labs(x = "% Diff. Score Cum. Haz.\nFrom Plain Cum. Haz.", y = "", color = "Analysis") +
  geom_hline(yintercept = 1:length(unique(diff_df$author)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".perc_diff_plain_dot.png"), the_plot, width = 6.8, height = 5)


#################################################


score_df <- do.call("rbind", all_score_perc_df)
score_df$risk_group <- stringr::str_to_title(score_df$risk_group)

score_df <- score_df[score_df$author %in% box_disease,]

the_plot <- ggplot(score_df[score_df$author %in% box_disease,], aes(type, val*100, color = risk_group)) + geom_boxplot() +
  labs(x = "Analysis", y = "% Diff. Score Model Cum. Haz.\nFrom Intermediate Risk", color = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".perc_diff_inter_box.png"), the_plot, width = 6.8, height = 5)


score_df$author <- factor(score_df$author, levels = score_df$author[score_df$risk_group == "High" & score_df$type == "Plain"][
  order(score_df$val[score_df$risk_group == "High" & score_df$type == "Plain"])])

the_plot <- ggplot(score_df[score_df$risk_group == "High",], aes(val, author, color = type)) + geom_point() +
  labs(x = "% Diff. Score Model Cum. Haz.\nFrom Intermediate Risk", y = "", color = "Analysis") +
  geom_hline(yintercept = 1:length(unique(score_df$author)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".perc_diff_inter_dot.png"), the_plot, width = 6.8, height = 5)


#################################################

imp_df <- do.call("rbind", all_imp_perc_df)
imp_df$risk_group <- stringr::str_to_title(imp_df$risk_group)

imp_df <- imp_df[imp_df$author %in% box_disease,]

the_plot <- ggplot(imp_df[imp_df$author %in% box_disease,], aes(type, val*100, color = risk_group)) + geom_boxplot() +
  labs(x = "Analysis", y = "% Diff. Imp. Cum. Haz.\nFrom Intermediate Risk", color = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_cartesian(ylim = c(min(imp_df$val)*100, quantile(imp_df$val, 0.95)*100))
plot(the_plot)

#meta plot of diff and perc diff from unchanged improvement of low, medium, high - just look at last time point