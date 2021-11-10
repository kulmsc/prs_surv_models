library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
theme_set(theme_cowplot())


#analysis_type <- "competing_risks"
analysis_type <- commandArgs(trailingOnly=TRUE)
res <- readRDS(paste0("fuel_for_plot/km.", analysis_type, ".RDS"))


source("common.R")



############################################

if(!("unchanged" %in% res[[1]]$type)){
  plain_data <- readRDS("fuel_for_plot/km.plain.1.RDS")
}

all_imp_df <- list()
all_diff_df <- list()
all_score_df <- list()


for(i in 1:length(res)){
  res[[i]]$author <- swap_names(res[[i]]$author, disease_names)
  res[[i]]$type <- swap_names(res[[i]]$type, model_names)
  res[[i]]$model <- stringr::str_to_title(res[[i]]$model)


  if(analysis_type == "censoring"){
    res[[i]]$split_up <- "From Assessment (Right Censor)"
    res[[i]]$split_up[grepl("Left", res[[i]]$type)] <- "From Age 60 (Left Censor)"
    res[[i]]$split_up <- as.factor(res[[i]]$split_up)
  } else if(analysis_type == "competing_risks"){
    res[[i]]$split_up <- "From Assessment"
    res[[i]]$split_up[grepl("iCARE", res[[i]]$type)] <- "From Age 60"
    res[[i]]$split_up <- as.factor(res[[i]]$split_up)
    res[[i]]$time[res[[i]]$type == "iCARE"] <- res[[i]]$time[res[[i]]$type == "iCARE"] *365
  }
  

  the_plot <- ggplot(res[[i]][res[[i]]$risk_group == "high" & res[[i]]$model == "Score",],
         aes(time/365, val, color = type)) +
    geom_point() + geom_line() +
    labs(x = "Year From Assessment", y = "Cumulative Hazard", color = "Analysis")
  if(analysis_type == "censoring" | analysis_type == "competing_risks"){
    the_plot <- the_plot + facet_wrap( ~ split_up, scales = "free", nrow = 2) + labs(x = "Year")
  }
  plot(the_plot)
  ggsave(paste0("out_plots/km/km.", analysis_type, ".over_time.", gsub(" ", "_", res[[i]]$author[1]), ".png"),
         the_plot, width = 6.3, height = 4)
  

  
  the_plot <- ggplot(get_top_time(res[[i]]), aes(type, val, color = risk_group, shape = model)) +
    geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
    labs(x = "Analysis", y = "Cumulative Hazard", color = "Risk Group", shape = "Model Type") +
    geom_errorbar(aes(ymin = val - sd, ymax = val + sd, width = 0), position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
    theme(axis.text=element_text(size=10)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    coord_cartesian(ylim = c(0, min(c(0.1, max(get_top_time(res[[i]])$val + get_top_time(res[[i]])$sd)*1.1))))
  plot(the_plot)
  ggsave(paste0("out_plots/km/km.", analysis_type, ".single_time.", gsub(" ", "_", res[[i]]$author[1]), ".png"),
         the_plot, width = 7.3, height = 5)
  
  
  
  temp_df <- get_top_time(res[[i]])
  imp_df <- temp_df[temp_df$model == "Score",]
  all_score_df[[i]] <- imp_df
  
  if("Plain" %in% imp_df$type){
    all_diff_df[[i]] <- imp_df[imp_df$type != "Plain",]
    all_diff_df[[i]]$val <- all_diff_df[[i]]$val - rep(imp_df$val[imp_df$type == "Plain"], nrow(all_diff_df[[i]])/3)
                                                       
  } else {
    all_diff_df[[i]] <- imp_df
    all_diff_df[[i]]$val <- all_diff_df[[i]]$val - rep(plain_data$val[plain_data$author == imp_df$author[1]], nrow(all_diff_df[[i]])/3)
                                                       
  }
  
  imp_df$val <- imp_df$val - temp_df$val[temp_df$model == "Base"]
  imp_df$model <- "diff"
  all_imp_df[[i]] <- imp_df
  
}



imp_df <- do.call("rbind", all_imp_df)
imp_df$risk_group <- stringr::str_to_title(imp_df$risk_group)

the_plot <- ggplot(imp_df, aes(type, val, color = risk_group)) + geom_boxplot() +
  labs(x = "Analysis", y = "Difference in Cum. Haz.\nScore Model - Base Model", color = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_cartesian(ylim = c(-0.1, min(c(0.5, max(imp_df$val)*1.1))))
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".imp_last_time.png"), the_plot, width = 6.8, height = 5)



score_df <- do.call("rbind", all_score_df)
score_df$risk_group <- stringr::str_to_title(score_df$risk_group)

the_plot <- ggplot(score_df, aes(type, val, color = risk_group)) + geom_boxplot() +
  labs(x = "Analysis", y = "Score Model Cum. Haz.", color = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".score_last_time.png"), the_plot, width = 6.8, height = 5)



diff_df <- do.call("rbind", all_diff_df)
diff_df$risk_group <- stringr::str_to_title(diff_df$risk_group)

the_plot <- ggplot(diff_df, aes(type, val, color = risk_group)) + geom_boxplot() +
  labs(x = "Analysis", y = "Score Model Cum. Haz. Minus\nPlain Model Cum. Haz.", color = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(the_plot)
ggsave(paste0("out_plots/meta/", analysis_type, "/km.", analysis_type, ".plaindiff_last_time.png"), the_plot, width = 6.8, height = 5)


#meta plot of diff and perc diff from unchanged improvement of low, medium, high - just look at last time point