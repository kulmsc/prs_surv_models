

disease_names <- read.table("helper_files/disease_names", stringsAsFactors = F, sep = "\t")
if(analysis_type == "age_covar"){
  model_names <- data.frame("old" = c("unchanged", "square", "interaction", "dob", "age_match"),
                            "new" = c("Plain", "Square", "Interaction", "Date of Birth", "Age Match Cohort"),
                            stringsAsFactors = F)
}else if(analysis_type == "censoring"){
  model_names <- data.frame("old" = c("no_censor", "unchanged", "left_unchanged_censor", "ehr", "left_ehr_censor"),
                            "new" = c("No Censoring", "Plain", "Left Censor", "Censor at EHR", "Left Censor at EHR"),
                            stringsAsFactors = F)
}else if(analysis_type == "competing_risks"){
  model_names <- data.frame("old" = c("fg", "ms", "rr", "icare_gbd", "icare_holt", "icare_combo", "gbd_icare", "holt_icare", "combo_icare"),
                            "new" = c("Fine-Gray", "Multi-State", "Risk Regression", "iCARE: GBD", "iCARE: Holt", "iCARE: Combo", "iCARE: GBD", "iCARE: Holt", "iCARE: Combo"),
                            stringsAsFactors = F)
}else if(analysis_type == "disease_labels"){
  model_names <- data.frame("old" = c("reclassify_950", "reclassify_990", "reclassify_995",
                                      "remove_950", "remove_990", "remove_995"),
                            "new" = c("Reclassify: 95%", "Reclassify: 99%", "Reclassify: 99.5%",
                                      "Remove: 95%", "Remove: 99%", "Remove: 99.5%"),
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

