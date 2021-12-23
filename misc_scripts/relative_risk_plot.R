library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
library(stringr)
theme_set(theme_cowplot())

analysis_type <- "censoring"
#analysis_type <- commandArgs(trailingOnly=TRUE)
conc_limit <- 0.75

res <- readRDS(paste0("fuel_for_plot/relative_risk.", analysis_type, ".RDS"))

if(analysis_type == "adjustment" | analysis_type == "disease_labels"){
  source("check_acc.R")
}

source("common.R")


#################################################################
#diff is the prevalence for individuals unique to the score compared to base group (or vice versa)
stat_types <- c("Prevalence", "log(Odds Ratio)", "Reclassified Group Prevalence")
save_types <- c("prev", "or", "rgp")

all_res <- list("prev" = list(), "or" = list(), "diff" = list())
for(i in 1:length(res[[1]])){
  
  for(j in 1:3){
    
    #if(j == 3){
     res[[j]][[i]] <- res[[j]][[i]][!(res[[j]][[i]]$model %in% c("left_ehr_censor", "ehr", "just1_coef", "just2_coef", "just_3_coef")),]
    #} else {
    #  res[[j]][[i]] <- res[[j]][[i]][!(res[[j]][[i]]$type %in% c("left_ehr_censor", "ehr")),]
    #}
    
    plot_df <- res[[j]][[i]]
    
    if(j == 2){
      colnames(plot_df)[colnames(plot_df) == "hi_cutoff"] <- "cutoff"
      plot_df$or <- log(plot_df$or)
      plot_df$low_or <- log(plot_df$low_or)
      plot_df$hi_or <- log(plot_df$hi_or)
      print(min(plot_df$or))
      plot_df$or[is.nan(plot_df$or) | is.infinite(plot_df$or)] <- NA
      plot_df$low_or[is.nan(plot_df$low_or) | is.infinite(plot_df$low_or)] <- NA
      plot_df$hi_or[is.nan(plot_df$hi_or) | is.infinite(plot_df$hi_or)] <- NA
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
      geom_point() +
      #geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.3)) +
      theme(axis.text=element_text(size=10)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      labs(x = "Analysis", y = stat_types[j], color = "Predicted\nGroup\nCut-Off", shape = "Model")
    options(warn=0)
     plot(the_plot)
     
    ggsave(paste0("out_plots/relative_risk/relative_risk.", analysis_type, ".",
                  gsub(" ", "_", swap_names(names(res[[1]])[i], disease_names)), ".png"),
           the_plot, width = 4.7, height = 5)
    #options(warn=2)
    

    plot_df$disease <- swap_names(names(res[[j]])[i], disease_names)
    all_res[[j]][[i]] <- plot_df
  
  }
  
}





stat_types <- c("Prevalence", "log(Odds Ratio)", "Reclassified Group Prevalence")
for(i in 1:3){
  
  
  prev_df <- do.call("rbind", all_res[[i]])
  
  if(analysis_type == "age_covar" & i == 2){
    plain_supp <- prev_df[prev_df$model == "Plain" & prev_df$type == "Score" & prev_df$cutoff == 0.95,][,c(8,2,1,3)]
    plain_supp[,2] <- signif(plain_supp[,2])
    plain_supp[,3] <- signif(plain_supp[,3])
    plain_supp[,4] <- signif(plain_supp[,4])
    plain_supp$imp <- signif(prev_df[prev_df$model == "Plain" & prev_df$type == "Score" & prev_df$cutoff == 0.95,]$val -
                             prev_df[prev_df$model == "Plain" & prev_df$type == "Base" & prev_df$cutoff == 0.95,]$val, 3)
    write.table(plain_supp, "supp_tables/rr_or.plain.txt", row.names = F, col.names = T, quote = F, sep = "\t")
  }
  

  
  if(analysis_type == "competing_risks" | analysis_type == "adjustment"){
    box_disease <- names(table(prev_df$disease)[table(prev_df$disease) == max(table(prev_df$disease))])
  } else {
    box_disease <- unique(prev_df$disease)
  }
  

  if(analysis_type %in% c("disease_labels", "adjustment", "competing_risks")){
    plain_data <- readRDS(paste0("fuel_for_plot/relative_risk.plain.", i, ".RDS"))
    plain_data <- plain_data[plain_data$model == "Plain" & plain_data$disease %in% prev_df$disease ,]
    # if(i == 2){
    #   plain_data$val <- log(plain_data$val)
    #   plain_data$val[is.infinite(plain_data$val) | is.nan(plain_data$val)] <- NA
    #   plain_data$low_val <- log(plain_data$low_val)
    #   plain_data$low_val[is.infinite(plain_data$low_val) | is.nan(plain_data$low_val)] <- NA
    #   plain_data$hi_val <- log(plain_data$hi_val)
    #   plain_data$hi_val[is.infinite(plain_data$hi_val) | is.nan(plain_data$hi_val)] <- NA
    # }
    prev_df <- rbind(prev_df, plain_data)
  }
  prev_df$model <- factor(prev_df$model, levels = c(unique(prev_df$model[prev_df$model != "Plain"]), "Plain"))
  
  prev_df <- prev_df[prev_df$cutoff == 0.95,]
  
  
  
  the_plot <- ggplot(prev_df[prev_df$disease %in% box_disease,], aes(x = model, y = val, color = type)) + 
    geom_boxplot() + 
    geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width = 0, jitter.height = 0)) +
    theme(axis.text=element_text(size=10)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Analysis", y = stat_types[i], color = "Model")
  options(warn=0)
  plot(the_plot)
  ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".box.png"),
         the_plot, width = 5, height = 5)
  #options(warn=2)
  save_abs_stat_val <- data.frame(prev_df[prev_df$type == "Score",])
  
  
  
  if( save_types[i] != "rgp"){
    spec_df <- prev_df[prev_df$type == "Score",]
    spec_df$val <- spec_df$val - prev_df$val[prev_df$type == "Base"]
    plain_spec_df <- spec_df[spec_df$model != "Plain",]
    plain_spec_df$plain_val <- plain_spec_df$val - custom_rep(spec_df$val[spec_df$model == "Plain"], plain_spec_df$disease)
    plain_spec_df$perc_diff <- plain_spec_df$plain_val
    plain_spec_df$perc_diff <- plain_spec_df$perc_diff/custom_rep(spec_df$val[spec_df$model == "Plain"], plain_spec_df$disease)
    
    mean_vals <- unlist(lapply(unique(spec_df$model), function(x) mean(spec_df$val[spec_df$model == x], na.rm=T)))
    spec_df$model <- factor(spec_df$model, levels = unique(spec_df$model)[order(mean_vals)])
    the_plot <- ggplot(spec_df[spec_df$disease %in% box_disease,], aes(x = model, y = val)) + 
      geom_boxplot() + 
      geom_point() +
      theme(axis.text=element_text(size=10)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      labs(x = "Analysis", y = paste("Imp.", stat_types[i]), color = "Model")
    plot(the_plot)
    options(warn=0) 
    ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".IMP_box.png"),
           the_plot, width = 4.8, height = 3.8)
    save_imp_stat_val <- data.frame(spec_df)
    
    
    mean_vals <- unlist(lapply(unique(plain_spec_df$model), function(x) mean(plain_spec_df$val[plain_spec_df$model == x], na.rm=T)))
    plain_spec_df$model <- factor(plain_spec_df$model, levels = unique(plain_spec_df$model)[order(mean_vals)])
    the_plot <- ggplot(plain_spec_df[plain_spec_df$disease %in% box_disease,], aes(x = model, y = perc_diff*100)) + 
      geom_boxplot() + 
      geom_point() +
      theme(axis.text=element_text(size=10)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      labs(x = "Analysis", y = paste("Imp.", stat_types[i], "\n% Diff From Plain Analysis"), color = "Model")
    if(analysis_type == "age_covar"){
      the_plot <- the_plot + coord_cartesian(ylim = c(-150, 250))
    }
    # } else if(analysis_type  == "disease_labels"){
    #   the_plot <- the_plot + coord_cartesian(ylim = c(-175, 225))
    # }
    plot(the_plot)
    options(warn=0) 
    ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".IMP_diff_box.png"),
           the_plot, width = 4.2, height = 4)
    save_diff_imp_stat_val <- data.frame(plain_spec_df)
    save_diff_imp_stat_val$perc_diff <- save_diff_imp_stat_val$perc_diff*100 
    
    
    plain_spec_df$disease <- as.character(plain_spec_df$disease)
    #plain_spec_df$disease <- factor(plain_spec_df$disease, levels = 
    #        plain_spec_df$disease[plain_spec_df$model == "Multi-State"][
    #order(plain_spec_df$perc_diff[plain_spec_df$model == "Multi-State"])])
    the_plot <- ggplot(plain_spec_df, aes(x = perc_diff*100, y = disease, color = model)) + 
      geom_vline(aes(xintercept = 0)) + geom_point() +
      theme(axis.text=element_text(size=10)) +
      geom_hline(yintercept = (1:length(unique(prev_df$disease)))+0.5, color = "grey80") +
      labs(color = "", x = paste("Imp.", stat_types[i], "\n% Diff From Plain Analysis"), y="")
    #the_plot <- the_plot + coord_cartesian(ylim = c(-150, 250))
    plot(the_plot)
    options(warn=0) 
    ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".IMP_diff_dot.png"),
           the_plot, width = 6.8, height = 5)
    
    
    
    
    #options(warn=2) 
    if(analysis_type == "adjustment"){
      spec_df$model = as.character(spec_df$model)
      spec_df$model <- factor(spec_df$model, levels = c(sort(unique(spec_df$model))[c(1,2,4,3)], "Plain"))
    }
    temp <- spec_df$disease[spec_df$model == spec_df$model[1]][order(spec_df$val[spec_df$model == spec_df$model[1]])]
    temp <- c(temp, unique(spec_df$disease)[!(unique(spec_df$disease) %in% temp)])
    spec_df$disease <- factor(spec_df$disease, levels = temp)
    the_plot <- ggplot(spec_df, aes(x = val, y = disease, color = model)) + 
      geom_point() + geom_vline(xintercept = 0) +
      geom_hline(yintercept = (1:length(unique(prev_df$disease)))+0.5, color = "grey80") +
      labs(x = paste("Imp.", stat_types[i]), y = "",  color = "Model") + theme(axis.text=element_text(size=10)) 
    plot(the_plot)
    options(warn=0) 
    ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".IMP_dot.png"),
           the_plot, width = 5.6, height = 5)
    if(analysis_type == "adjustment"){
      the_plot <- the_plot + scale_color_manual(values = gg_color_hue(5))
      ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".IMP_dot.png"),
             the_plot, width = 6, height = 4)
    }
    #options(warn=2) 
  }
  
  
  
  sub_prev_df <- prev_df[prev_df$type == "Score",]
  sub_prev_df$disease <- factor(sub_prev_df$disease, levels = sub_prev_df$disease[sub_prev_df$model == "Plain"][
    order(sub_prev_df$val[sub_prev_df$model == "Plain"])])
  the_plot <- ggplot(sub_prev_df, aes(x = val, y = disease, color = model)) + 
    geom_point() +
    geom_hline(yintercept = (1:length(unique(prev_df$disease)))+0.5, color = "grey80") +
    labs(x = stat_types[i], y = "",  color = "Model") + theme(axis.text=element_text(size=10)) 
  options(warn=0) 
   plot(the_plot)
  ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".dot.png"),
         the_plot, width = 7, height = 5)
  #options(warn=2)

  
  
  dubsub_prev_df <- sub_prev_df[sub_prev_df$model != "Plain",]
  dubsub_prev_df$val <- dubsub_prev_df$val - custom_rep(sub_prev_df$val[sub_prev_df$model == "Plain"], dubsub_prev_df$disease)
  mod_mean <- unlist(lapply(unique(dubsub_prev_df$model), function(x) mean(dubsub_prev_df$val[dubsub_prev_df$model == x])))
  dubsub_prev_df$model <- factor(dubsub_prev_df$model, unique(dubsub_prev_df$model)[order(mod_mean)])
                                                 
  the_plot <- ggplot(dubsub_prev_df[dubsub_prev_df$disease %in% box_disease,], aes(x = model, y = val)) + 
     geom_boxplot() + geom_point() +
    labs(y = paste0(stat_types[i], "\nDiff. From Plain Analysis"), x = "Analysis") + theme(axis.text=element_text(size=10)) +
    theme(axis.text=element_text(size=10)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  options(warn=0)
  plot(the_plot)
  ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".minus_score.png"),
         the_plot, width = 4, height = 5)
  #options(warn=2)
  

  dubsub_prev_df <- sub_prev_df[sub_prev_df$model != "Plain",]
  dubsub_prev_df$val <- (dubsub_prev_df$val - custom_rep(sub_prev_df$val[sub_prev_df$model == "Plain"], dubsub_prev_df$disease))/
    custom_rep(sub_prev_df$val[sub_prev_df$model == "Plain"], dubsub_prev_df$disease)
  dubsub_prev_df$val <- dubsub_prev_df$val * 100
  
  dubsub_prev_df$disease <- as.character(dubsub_prev_df$disease)
  temp <- dubsub_prev_df$disease[dubsub_prev_df$model == dubsub_prev_df$model[1]][
    order(dubsub_prev_df$val[dubsub_prev_df$model == dubsub_prev_df$model[1]])]
  temp <- c(temp, unique(dubsub_prev_df$disease)[!(unique(dubsub_prev_df$disease) %in% temp)])
  dubsub_prev_df$disease <- factor(dubsub_prev_df$disease, levels = temp)
  the_plot <- ggplot(dubsub_prev_df, aes(x = val, y = disease, color = model)) + 
    geom_vline(xintercept = 0) + geom_point() +
    geom_hline(yintercept = (1:length(unique(prev_df$disease)))+0.5, color = "grey80") +
    labs(x = paste0(stat_types[i], " % Diff. From\nPlain Analysis"), y = "",  color = "Model") + theme(axis.text=element_text(size=10)) 
  options(warn=0)
   plot(the_plot)
  ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".dot_perc_diff_plain.png"),
         the_plot, width = 6.8, height = 5)
  save_diff_abs_stat_val <- data.frame(dubsub_prev_df)
  
  mod_mean <- unlist(lapply(unique(dubsub_prev_df$model), function(x) mean(dubsub_prev_df$val[dubsub_prev_df$model == x])))
  dubsub_prev_df$model <- factor(dubsub_prev_df$model, unique(dubsub_prev_df$model)[order(mod_mean)])
  the_plot <- ggplot(dubsub_prev_df, aes(y = val, x = model)) + 
    geom_boxplot() + geom_point() +
    labs(y = paste0(stat_types[i], "\n% Diff. From Plain Analysis"), x = "",  color = "Model") + 
    theme(axis.text=element_text(size=10)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  #the_plot <- the_plot + coord_cartesian(ylim = c(-120, 50))
  plot(the_plot)
  ggsave(paste0("out_plots/meta/", analysis_type, "/relative_risk.", analysis_type, ".", save_types[i], ".box_perc_diff_plain.png"),
         the_plot, width = 4.2, height = 4)
  
  

  final_save <- join_save(save_abs_stat_val[,c("model", "disease", "val")],
            save_imp_stat_val[,c("model", "disease", "val")],
            save_diff_abs_stat_val[,c("model", "disease", "val")],
            save_diff_imp_stat_val[,c("model", "disease", "perc_diff")], toupper(save_types[i]))
  if(analysis_type == "censoring"){
    final_save[,1] <- gsub("\n", " ", as.character(final_save[,1]))
  }
  write.table(final_save, paste0("supp_tables/relative_risk_", save_types[i], ".", analysis_type, ".txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
  #options(warn=2)
}
