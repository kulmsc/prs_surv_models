library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
theme_set(theme_cowplot())

#analysis_type <- "competing_risks"
analysis_type <- commandArgs(trailingOnly=TRUE)

res <- readRDS(paste0("fuel_for_plot/relative_risk.", analysis_type, ".RDS"))

source("common.R")


#################################################################
#diff is the prevalence for individuals unique to the score compared to base group (or vice versa)
stat_types <- c("Prevalence", "Odds Ratio", "Reclassified Group Prevalence")
save_types <- c("prev", "or", "rgp")

all_res <- list("prev" = list(), "or" = list(), "diff" = list())
for(i in 1:length(res[[1]])){
  
  for(j in 1:3){
    
    plot_df <- res[[j]][[i]]
    if(j == 2){
      colnames(plot_df)[colnames(plot_df) == "hi_cutoff"] <- "cutoff"
    } else if(j == 3){
      plot_df <- melt(plot_df, id.vars = c("cutoff", "model"))
      colnames(plot_df)[3:4] <- c("type", "val")
      plot_df$type <- unlist(lapply(strsplit(as.character(plot_df$type), "_"), function(x) x[1]))
    }
    
    plot_df$model <- swap_names(plot_df$model, model_names)
    plot_df$type <- stringr::str_to_title(plot_df$type)
    plot_df$cutoff <- as.factor(plot_df$cutoff)
    
    #fix names
    if(j==1){
      plot_df$hi_val <- plot_df$prev
      plot_df$low_val <- plot_df$prev
      colnames(plot_df)[1] <- "val"
    } else if(j == 2){
      colnames(plot_df)[1:3] <- c("low_val", "val", "hi_val")
    } else if(j == 3){
      plot_df$hi_val <- plot_df$val
      plot_df$low_val <- plot_df$val
    }
    
    the_plot <- ggplot(plot_df, aes(model, val, color = cutoff, shape = type)) + 
      geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.3)) +
      theme(axis.text=element_text(size=10)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      labs(x = "Analysis", y = stat_types[j], color = "Predicted\nGroup\nCut-Off", shape = "Model")
    plot(the_plot)
    ggsave(paste0("out_plots/relative_risk/relative_risk.", analysis_type, ".",
                  gsub(" ", "_", swap_names(names(res[[1]])[i], disease_names)), ".png"),
           the_plot, width = 7.5, height = 5)
    
    plot_df$disease <- swap_names(names(res[[j]])[i], disease_names)
    all_res[[j]][[i]] <- plot_df
  
  }
  
}




stat_types <- c("Prevalence", "Odds Ratio", "Reclassified Group Prevalence Diff.")
for(i in 1:3){
  
  prev_df <- do.call("rbind", all_res[[i]])
  if(analysis_type %in% c("disease_labels", "adjustment", "competing_risks")){
    plain_data <- readRDS(paste0("fuel_for_plot/relative_risk.plain.", i, ".RDS"))
    plain_data <- plain_data[plain_data$model == "Plain",]
    prev_df <- rbind(prev_df, plain_data)
  }
  prev_df$model <- factor(prev_df$model, levels = c(unique(prev_df$model[prev_df$model != "Plain"]), "Plain"))
  
  prev_df <- prev_df[prev_df$cutoff == 0.95,]
  
  
  
  the_plot <- ggplot(prev_df, aes(x = model, y = val, color = type)) + 
    geom_boxplot() + 
    geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width = 0, jitter.height = 0)) +
    theme(axis.text=element_text(size=10)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Analysis", y = stat_types[i], color = "Model")
  plot(the_plot)
  ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".box.png"),
         the_plot, width = 6.8, height = 5)
  

  
  
  sub_prev_df <- prev_df[prev_df$type == "Score",]
  the_plot <- ggplot(sub_prev_df, aes(x = val, y = disease, color = model)) + 
    geom_point() +
    geom_hline(yintercept = (1:length(unique(prev_df$disease)))+0.5, color = "grey80") +
    labs(x = stat_types[i], y = "",  color = "Model") + theme(axis.text=element_text(size=10)) 
  plot(the_plot)
  ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".dot.png"),
         the_plot, width = 11, height = 5)

  #go back and check each is not in rep
  
  
  dubsub_prev_df <- sub_prev_df[sub_prev_df$model != "Plain",]
  dubsub_prev_df$val <- dubsub_prev_df$val - custom_rep(sub_prev_df$val[sub_prev_df$model == "Plain"], dubsub_prev_df$disease)
                                                 
  the_plot <- ggplot(dubsub_prev_df, aes(x = model, y = val)) + 
    geom_point() + geom_boxplot() +
    labs(y = stat_types[i], x = "Analysis") + theme(axis.text=element_text(size=10)) +
    theme(axis.text=element_text(size=10)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  plot(the_plot)
  ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".minus_score.png"),
         the_plot, width = 5, height = 5)
  
  #need to get plain minus
}
