

disease_names <- read.table("helper_files/disease_names", stringsAsFactors = F, sep = "\t")

if(analysis_type == "age_covar" | analysis_type == "age_covar_AGE"){
  model_names <- data.frame("old" = c("unchanged", "square", "interaction", "dob", "age_match", "age_small", "age_strat", "sex_strat"),
                            "new" = c("Plain", "Square", "Interaction", "Date of Birth", "Age Match Cohort", "Narrow Age Cohort", "Age Stratified", "Sex Stratified"),
                            stringsAsFactors = F)
}else if(analysis_type == "censoring"){
  model_names <- data.frame("old" = c("no_censor", "unchanged", "left_unchanged_censor", "ehr", "left_ehr_censor", "delay_entry"),
                            "new" = c("No Censoring", "Plain", "Delay Entry by\nAge Time", "Censor at EHR", "Left Censor at EHR", "Delay Entry by\nStudy Time"),
                            stringsAsFactors = F)
}else if(analysis_type == "competing_risks"){
  model_names <- data.frame("old" = c("fg", "ms", "rr", "icare_gbd", "icare_holt", "icare_forecast", "icare_median", "gbd_icare", "holt_icare", "forecast_icare", "median_icare"),
                            "new" = c("Fine-Gray", "Multi-State", "Risk Regression", "iCARE: GBD", "iCARE: Holt", "iCARE: Forecast", "iCARE: Median", "iCARE: GBD", "iCARE: Holt", "iCARE: Forecast", "iCARE: Median"),
                            stringsAsFactors = F)
}else if(analysis_type == "disease_labels"){
  model_names <- data.frame("old" = c("reclassify_950", "reclassify_990", "reclassify_995",
                                      "remove_950", "remove_990", "remove_995"),
                            "new" = c("Reclassify: 5%", "Reclassify: 1%", "Reclassify: 0.5%",
                                      "Remove: 5%", "Remove: 1%", "Remove: 0.5%"),
                            stringsAsFactors = F)
}else if(analysis_type == "adjustment"){
  model_names <- data.frame("old" = c("top3_coef", "top5_coef", "top10_coef", "just1_coef", "just2_coef", "just_3_coef", "big_pred"),
                            "new" = c("Incl. Top 3 Feats.", "Incl. Top 5 Feats.", "Incl. Top 10 Feats.",
				      "Incl 1st Feat.", "Incl. 2nd Feat.", "Incl. 3rd Feat.", "Incl. CoxNet Pred."),
                            stringsAsFactors = F)
}


swap_names <- function(old, key){
  for(x in unique(old)){
    old[old == x] <- key[key[,1] == x, 2]
  }
  return(old)
}

better_swap_names <- function(old, key){
  for(x in unique(key[,1])){
    old[old == x] <- key[key[,1] == x, 2]
  }
  return(old)
}


get_top_time <- function(df){
  out_list <- list()
  for(i in 1:length(unique(df$type))){
    out_list[[i]] <- df[df$time == max(df$time[df$type == unique(df$type)[i]]) & df$type == unique(df$type)[i],]
  }
  return(do.call("rbind", out_list))
}

custom_rep <- function(x, y){
  z <- list()
  for(i in 1:length(x)){
    z[[i]] <- rep(x[i], sum(y == unique(y)[i] ))
  }
  return(unlist(z))
}

custom_tri_rep <- function(x, y){
  z <- list()
  for(i in 1:length(x)){
    z[[i]] <- rep(x[i], sum(y == unique(y)[i] ))
  }
  return(unlist(z))
}

scale_by_intermediate <- function(plot_df){
  for(utype in unique(plot_df$type)){
    for(utime in seq(365, 365*14, 365)){
      for(urisk in c("high", "low")) {
        x <- plot_df$val[plot_df$time <= utime & plot_df$time > utime-300 & plot_df$type == utype & plot_df$risk_group == urisk]
        y <- plot_df$val[plot_df$time <= utime & plot_df$time > utime-300 & plot_df$type == utype & plot_df$risk_group == "inter"]

        plot_df$val[plot_df$time <= utime & plot_df$time > utime-300 & plot_df$type == utype & plot_df$risk_group == urisk] <- (x-y)/abs(y)

      }
    }
  }
  return(plot_df)
}


join_save <- function(abs_df, imp_df, diff_abs_df, diff_imp_df, stat_name){
  abs_df <- abs_df[abs_df[,1] != "Plain",]
  imp_df <- imp_df[imp_df[,1] != "Plain",]

  colnames(abs_df) <- c("type", "disease", "abs_val")
  colnames(imp_df) <- c("type", "disease", "imp_val")
  colnames(diff_abs_df) <- c("type", "disease", "diff_abs_val")
  colnames(diff_imp_df) <- c("type", "disease", "diff_imp_val")

  abs_df$combo <- paste(as.character(abs_df[,1]) , "-", as.character(abs_df[,2]))
  imp_df$combo <- paste(as.character(imp_df[,1]) , "-", as.character(imp_df[,2]))
  diff_abs_df$combo <- paste(as.character(diff_abs_df[,1]) , "-", as.character(diff_abs_df[,2]))
  diff_imp_df$combo <- paste(as.character(diff_imp_df[,1]) , "-", as.character(diff_imp_df[,2]))


  togo <- merge(abs_df, imp_df, by = "combo", all = T)
  print(colnames(togo))
  togo <- togo[,-c(which(grepl(".x", colnames(togo), fixed = T)), which(grepl(".y", colnames(togo), fixed = T)))]
  print(colnames(togo))
  togo <- merge(togo, diff_abs_df, by = "combo", all = T)
  print(colnames(togo))
  #togo <- togo[,-c(which(grepl(".x", colnames(togo), fixed = T)), which(grepl(".y", colnames(togo), fixed = T)))]
  print(colnames(togo))
  togo <- merge(togo, diff_imp_df, by = "combo", all = T)
  print(colnames(togo))
  togo <- togo[,c("combo", "abs_val", "imp_val", "diff_abs_val", "diff_imp_val")]

  togo <- data.frame(str_split(togo$combo, " - ", simplify = T), togo[,2:5])
  colnames(togo) <- c("Analysis", "Disease", paste("Abs.", stat_name), paste("Imp.", stat_name),
		      paste("Diff. Abs.", stat_name), paste("Diff. Imp.", stat_name))
  togo[,3] <- signif(togo[,3], 3)
  togo[,4] <- signif(togo[,4], 3)
  togo[,5] <- signif(togo[,5], 3)
  togo[,6] <- signif(togo[,6], 3)

  return(togo)

}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
