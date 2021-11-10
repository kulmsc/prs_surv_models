library(epitools)

get_or <- function(vec, iden, hi_cut, low_cut){
  if(length(vec) != length(iden)){
    print("MISMATCH")
    exit()
  }

  group <- rep(1, length(vec))
  group[vec < quantile(vec, low_cut)] <- 0
  group[vec > quantile(vec, hi_cut)] <- 2

  odds_table <- matrix(c(sum(iden == 1 & group == 2), sum(iden == 0 & group == 2),
                         sum(iden == 1 & group == 0), sum(iden == 0 & group == 0)), nrow = 2)

  odds_ratio <- oddsratio.wald(odds_table)

  return(c(odds_ratio$measure[2,c(2,1,3)], hi_cut, low_cut))
}

get_prev <- function(vec, iden, hi_cut){
  if(length(vec) != length(iden)){
    print("MISMATCH")
    exit()
  }

  prev_val <- sum(iden[vec > quantile(vec, hi_cut)] == 1)/sum(vec > quantile(vec, hi_cut))
  return(c(prev_val, hi_cut))
}


get_diff <- function(base_vec, score_vec, iden, hi_cut){
  if(length(score_vec) != length(iden) | length(base_vec) != length(iden)){
    print("MISMATCH")
    exit()
  }

  base_ind <- which(base_vec > quantile(base_vec, hi_cut))
  score_ind <- which(score_vec > quantile(score_vec, hi_cut))

  u_base_ind <- base_ind[!(base_ind %in% intersect(base_ind, score_ind))]
  u_score_ind <- score_ind[!(score_ind %in% intersect(base_ind, score_ind))]

  base_prev <- sum(iden[u_base_ind])/length(u_base_ind)
  score_prev <- sum(iden[u_score_ind])/length(u_score_ind)

  return(c(score_prev, base_prev, hi_cut))
}


get_surv <- function(author){
  #read in the UKBB data
  surv_df <- readRDS(paste0("../init_data/surv_data/survdf.", author, ".RDS"))
  surv_df$time <- as.numeric(as.Date(surv_df$end_date) - as.Date(surv_df$date_attend)) #time is counted from the date of assessment

  surv_df$age <- as.numeric(surv_df$date_attend - surv_df$dob)
  surv_df$age_start_time <- as.numeric(surv_df$date_attend - surv_df$dob)
  surv_df$age_end_time <-  as.numeric(as.Date(as.character(surv_df$end_date)) - surv_df$dob)


  #split the dataset
  train_frac <- 0.6
  test_frac <- 1 - train_frac
  train_eid <- read.table(paste0("~/athena/doc_score/qc/cv_files/train_eid.", train_frac, ".txt"), stringsAsFactors=F)
  test_eid <- read.table(paste0("~/athena/doc_score/qc/cv_files/test_eid.", test_frac, ".txt"), stringsAsFactors=F)
  train_surv_df <- surv_df[surv_df$eid %in% train_eid[,1],]
  test_surv_df <- surv_df[surv_df$eid %in% test_eid[,1],]

  return(list(train_surv_df, test_surv_df))
}





#####################################################################
#####################################################################
#####################################################################
#age
print("age")

all_files <- list.files("../age_covar/res")
all_author <- rep(NA, length(all_files))
all_prev <- list()
all_or <- list()
all_diff <- list()

for(i in 1:length(all_files)){
  res <- readRDS(paste0("../age_covar/res/", all_files[i]))
  author_name <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  all_author[i] <- author_name
  surv_dfs <- get_surv(author_name)

  k <- 1
  temp_prev_base <- list()
  temp_prev_score <- list()
  temp_or_base <- list()
  temp_or_score <- list()
  temp_diff <- list()


  for(j in 1:length(res)){


#can calculate the typical odds ratio
#prevalence of high risk subset
#then also look at the difference in prevalence between those moved out of the score risk group compared to the base
    for(hicut in c(0.9, 0.95, 0.99)){

      temp_prev_base[[k]] <-  get_prev(res[[j]][["base"]][["risk"]][,14], surv_dfs[[2]]$pheno, hicut)
      temp_prev_score[[k]] <-  get_prev(res[[j]][["score"]][["risk"]][,14], surv_dfs[[2]]$pheno, hicut)

      temp_or_base[[k]] <- get_or(res[[j]][["base"]][["risk"]][,14], surv_dfs[[2]]$pheno, hicut, 0.5)
      temp_or_score[[k]] <- get_or(res[[j]][["score"]][["risk"]][,14], surv_dfs[[2]]$pheno, hicut, 0.5)

      temp_diff[[k]] <- get_diff(res[[j]][["base"]][["risk"]][,14], res[[j]][["score"]][["risk"]][,14], surv_dfs[[2]]$pheno, hicut)

      k <- k + 1

    }
  }

  all_prev[[i]] <- as.data.frame(rbind(do.call("rbind", temp_prev_base), do.call("rbind", temp_prev_score)))
  colnames(all_prev[[i]]) <- c("prev", "cutoff")
  all_prev[[i]]$type <- rep(c("base", "score"), each = length(temp_prev_base))
  all_prev[[i]]$model <- rep(names(res), each = 3)

  all_or[[i]] <- as.data.frame(rbind(do.call("rbind", temp_or_base), do.call("rbind", temp_or_score)))
  colnames(all_or[[i]]) <- c("low_or", "or", "hi_or", "hi_cutoff", "lo_cutoff")
  all_or[[i]]$type <- rep(c("base", "score"), each = length(temp_or_base))
  all_or[[i]]$model <- rep(names(res), each = 3)


  all_diff[[i]] <- as.data.frame(do.call("rbind", temp_diff))
  colnames(all_diff[[i]]) <- c("score_diff", "base_diff", "cutoff")
  all_diff[[i]]$model <- rep(names(res), each = 3)


}


names(all_prev) <- all_author
names(all_or) <- all_author
names(all_diff) <- all_author

saveRDS(list("prev" = all_prev, "or" = all_or, "diff" = all_diff), "ready_to_plot/relative_risk.age_covar.RDS")







#####################################################################
#####################################################################
#####################################################################
#censoring
print("censoring")


all_files <- list.files("../censoring/res")
all_author <- rep(NA, length(all_files))
all_prev <- list()
all_or <- list()
all_diff <- list()

for(i in 1:length(all_files)){
  res <- readRDS(paste0("../censoring/res/", all_files[i]))
  author_name <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  all_author[i] <- author_name
  surv_dfs <- get_surv(author_name)

  k <- 1
  temp_prev_base <- list()
  temp_prev_score <- list()
  temp_or_base <- list()
  temp_or_score <- list()
  temp_diff <- list()


  for(j in 1:5){

    if(j < 4){
      use_phen <- res[["phen"]][[1]]
    } else if(j == 4){
      use_phen <- res[["phen"]][[2]]
    } else {
      use_phen <- res[["phen"]][[3]]
    }

    for(hicut in c(0.9, 0.95, 0.99)){

      max_col <- ncol(res[[j]][["base"]][["risk"]])
      temp_prev_base[[k]] <-  get_prev(res[[j]][["base"]][["risk"]][,max_col], use_phen, hicut)
      temp_prev_score[[k]] <-  get_prev(res[[j]][["score"]][["risk"]][,max_col], use_phen, hicut)

      temp_or_base[[k]] <- get_or(res[[j]][["base"]][["risk"]][,max_col], use_phen, hicut, 0.5)
      temp_or_score[[k]] <- get_or(res[[j]][["score"]][["risk"]][,max_col], use_phen, hicut, 0.5)

      temp_diff[[k]] <- get_diff(res[[j]][["base"]][["risk"]][,max_col], res[[j]][["score"]][["risk"]][,14], use_phen, hicut)

      k <- k + 1

    }
  }

  all_prev[[i]] <- as.data.frame(rbind(do.call("rbind", temp_prev_base), do.call("rbind", temp_prev_score)))
  colnames(all_prev[[i]]) <- c("prev", "cutoff")
  all_prev[[i]]$type <- rep(c("base", "score"), each = length(temp_prev_base))
  all_prev[[i]]$model <- rep(names(res)[1:5], each = 3)

  all_or[[i]] <- as.data.frame(rbind(do.call("rbind", temp_or_base), do.call("rbind", temp_or_score)))
  colnames(all_or[[i]]) <- c("low_or", "or", "hi_or", "hi_cutoff", "lo_cutoff")
  all_or[[i]]$type <- rep(c("base", "score"), each = length(temp_or_base))
  all_or[[i]]$model <- rep(names(res)[1:5], each = 3)


  all_diff[[i]] <- as.data.frame(do.call("rbind", temp_diff))
  colnames(all_diff[[i]]) <- c("score_diff", "base_diff", "cutoff")
  all_diff[[i]]$model <- rep(names(res)[1:5], each = 3)


}


names(all_prev) <- all_author
names(all_or) <- all_author
names(all_diff) <- all_author

saveRDS(list("prev" = all_prev, "or" = all_or, "diff" = all_diff), "ready_to_plot/relative_risk.censoring.RDS")




###############################################################################
###############################################################################
###############################################################################
###############################################################################
#competing risks
print("competing risks")


all_files <- list.files("../competing_risks/res")
all_author <- rep(NA, length(all_files))
all_prev <- list()
all_or <- list()
all_diff <- list()

for(i in 1:length(all_files)){
  res <- readRDS(paste0("../competing_risks/res/", all_files[i]))
  author_name <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  all_author[i] <- author_name
  surv_dfs <- get_surv(author_name)

  k <- 1
  temp_prev_base <- list()
  temp_prev_score <- list()
  temp_or_base <- list()
  temp_or_score <- list()
  temp_diff <- list()


  for(j in 1:length(res)){

    for(hicut in c(0.9, 0.95, 0.99)){

      #note that I can look at sub_base and sub_score, but Im not to keep things simple and consistent
      if(j < 4){
      max_col <- ncol(res[[j]][["base"]][["risk"]])
      temp_prev_base[[k]] <-  get_prev(res[[j]][["base"]][["risk"]][,max_col], surv_dfs[[2]]$pheno, hicut)
      temp_prev_score[[k]] <-  get_prev(res[[j]][["score"]][["risk"]][,max_col], surv_dfs[[2]]$pheno, hicut)

      temp_or_base[[k]] <- get_or(res[[j]][["base"]][["risk"]][,max_col], surv_dfs[[2]]$pheno, hicut, 0.5)
      temp_or_score[[k]] <- get_or(res[[j]][["score"]][["risk"]][,max_col], surv_dfs[[2]]$pheno, hicut, 0.5)

      temp_diff[[k]] <- get_diff(res[[j]][["base"]][["risk"]][,max_col], res[[j]][["score"]][["risk"]][,max_col], surv_dfs[[2]]$pheno, hicut)

      k <- k + 1
      } else {
        if(!is.null(res[[j]])){
          for(l in 1:3){
            max_col <- ncol(res[[j]][["base"]][[l]])
            temp_prev_base[[k]] <-  get_prev(res[[j]][["base"]][[l]][,max_col], surv_dfs[[2]]$pheno, hicut)
            temp_prev_score[[k]] <-  get_prev(res[[j]][["score"]][[l]][,max_col], surv_dfs[[2]]$pheno, hicut)

            temp_or_base[[k]] <- get_or(res[[j]][["base"]][[l]][,max_col], surv_dfs[[2]]$pheno, hicut, 0.5)
            temp_or_score[[k]] <- get_or(res[[j]][["score"]][[l]][,max_col], surv_dfs[[2]]$pheno, hicut, 0.5)

            temp_diff[[k]] <- get_diff(res[[j]][["base"]][[l]][,max_col], res[[j]][["score"]][[l]][,max_col], surv_dfs[[2]]$pheno, hicut)

            k <- k + 1
          }
        }
      }

    }
  }

  if(!is.null(res[[4]])){
    model_names <- c(names(res)[1:3], paste0(names(res)[4], "_", names(res[[4]][[1]])))
  } else {
    model_names <- names(res)[1:3]
  }

  all_prev[[i]] <- as.data.frame(rbind(do.call("rbind", temp_prev_base), do.call("rbind", temp_prev_score)))
  colnames(all_prev[[i]]) <- c("prev", "cutoff")
  all_prev[[i]]$type <- rep(c("base", "score"), each = length(temp_prev_base))
  all_prev[[i]]$model <- rep(model_names, each = 3)

  all_or[[i]] <- as.data.frame(rbind(do.call("rbind", temp_or_base), do.call("rbind", temp_or_score)))
  colnames(all_or[[i]]) <- c("low_or", "or", "hi_or", "hi_cutoff", "lo_cutoff")
  all_or[[i]]$type <- rep(c("base", "score"), each = length(temp_or_base))
  all_or[[i]]$model <- rep(model_names, each = 3)


  all_diff[[i]] <- as.data.frame(do.call("rbind", temp_diff))
  colnames(all_diff[[i]]) <- c("score_diff", "base_diff", "cutoff")
  all_diff[[i]]$model <- rep(model_names, each = 3)


}


names(all_prev) <- all_author
names(all_or) <- all_author
names(all_diff) <- all_author

saveRDS(list("prev" = all_prev, "or" = all_or, "diff" = all_diff), "ready_to_plot/relative_risk.competing_risks.RDS")






###############################################################################
###############################################################################
###############################################################################
###############################################################################
#disease labels
print("disease labels")


all_files <- list.files("../disease_labels/res", "RDS")
all_author <- rep(NA, length(all_files))
all_prev <- list()
all_or <- list()
all_diff <- list()
all_python_conc <- list()

for(i in 1:length(all_files)){
  res <- readRDS(paste0("../disease_labels/res/", all_files[i]))
  author_name <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  all_author[i] <- author_name

  all_python_conc[[i]] <- c("conc" = res[[7]][["conc"]][6], "std" = res[[7]][["conc"]][7])

  k <- 1
  temp_prev_base <- list()
  temp_prev_score <- list()
  temp_or_base <- list()
  temp_or_score <- list()
  temp_diff <- list()

  for(j in 1:6){
    for(hicut in c(0.9, 0.95, 0.99)){
      max_col <- ncol(res[[j]][["base"]][["risk"]])
      temp_prev_base[[k]] <-  get_prev(res[[j]][["base"]][["risk"]][,max_col], res[[j]][["pheno"]], hicut)
      temp_prev_score[[k]] <-  get_prev(res[[j]][["score"]][["risk"]][,max_col], res[[j]][["pheno"]], hicut)

      temp_or_base[[k]] <- get_or(res[[j]][["base"]][["risk"]][,max_col], res[[j]][["pheno"]], hicut, 0.5)
      temp_or_score[[k]] <- get_or(res[[j]][["score"]][["risk"]][,max_col], res[[j]][["pheno"]], hicut, 0.5)

      temp_diff[[k]] <- get_diff(res[[j]][["base"]][["risk"]][,max_col], res[[j]][["score"]][["risk"]][,14], res[[j]][["pheno"]], hicut)

      k <- k + 1
    }
  }

  all_prev[[i]] <- as.data.frame(rbind(do.call("rbind", temp_prev_base), do.call("rbind", temp_prev_score)))
  colnames(all_prev[[i]]) <- c("prev", "cutoff")
  all_prev[[i]]$type <- rep(c("base", "score"), each = length(temp_prev_base))
  all_prev[[i]]$model <- rep(names(res)[1:6], each = 3)

  all_or[[i]] <- as.data.frame(rbind(do.call("rbind", temp_or_base), do.call("rbind", temp_or_score)))
  colnames(all_or[[i]]) <- c("low_or", "or", "hi_or", "hi_cutoff", "lo_cutoff")
  all_or[[i]]$type <- rep(c("base", "score"), each = length(temp_or_base))
  all_or[[i]]$model <- rep(names(res)[1:6], each = 3)


  all_diff[[i]] <- as.data.frame(do.call("rbind", temp_diff))
  colnames(all_diff[[i]]) <- c("score_diff", "base_diff", "cutoff")
  all_diff[[i]]$model <- rep(names(res)[1:6], each = 3)
}

all_python_conc <- as.data.frame(do.call("rbind", all_python_conc))
all_python_conc$author <- all_author

names(all_prev) <- all_author
names(all_or) <- all_author
names(all_diff) <- all_author

saveRDS(list("prev" = all_prev, "or" = all_or, "diff" = all_diff, "python" = all_python_conc), "ready_to_plot/relative_risk.disease_labels.RDS")






###############################################################################
###############################################################################
###############################################################################
###############################################################################
#adjustment
print("adjustment")


all_files <- list.files("../adjustment/res", "RDS")
all_author <- rep(NA, length(all_files))
all_prev <- list()
all_or <- list()
all_diff <- list()
all_python_conc <- list()

for(i in 1:length(all_files)){
  res <- readRDS(paste0("../adjustment/res/", all_files[i]))
  author_name <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  all_author[i] <- author_name
  #surv_dfs <- get_surv(author_name)

  all_python_conc[[i]] <- c("conc" = res[[8]][["conc"]][6], "std" = res[[8]][["conc"]][7])

  k <- 1
  temp_prev_base <- list()
  temp_prev_score <- list()
  temp_or_base <- list()
  temp_or_score <- list()
  temp_diff <- list()

  non_na_names <- rep("", 7)
 for(j in 1:7){
    if(class(res[[j]]) == "list" ){
      if(!is.null(res[[j]][["score"]])){

    non_na_names[j] <- names(res)[j]
    for(hicut in c(0.9, 0.95, 0.99)){
      max_col <- ncol(res[[j]][["base"]][["risk"]])
      temp_prev_base[[k]] <-  get_prev(res[[j]][["base"]][["risk"]][,max_col], res[["pheno"]], hicut)
      temp_prev_score[[k]] <-  get_prev(res[[j]][["score"]][["risk"]][,max_col], res[["pheno"]], hicut)

      temp_or_base[[k]] <- get_or(res[[j]][["base"]][["risk"]][,max_col], res[["pheno"]], hicut, 0.5)
      temp_or_score[[k]] <- get_or(res[[j]][["score"]][["risk"]][,max_col], res[["pheno"]], hicut, 0.5)

      temp_diff[[k]] <- get_diff(res[[j]][["base"]][["risk"]][,max_col], res[[j]][["score"]][["risk"]][,14], res[["pheno"]], hicut)

      k <- k + 1
    }
    }
    }
  }

  non_na_names <- non_na_names[non_na_names != ""]
  all_prev[[i]] <- as.data.frame(rbind(do.call("rbind", temp_prev_base), do.call("rbind", temp_prev_score)))
  colnames(all_prev[[i]]) <- c("prev", "cutoff")
  all_prev[[i]]$type <- rep(c("base", "score"), each = length(temp_prev_base))
  all_prev[[i]]$model <- rep(non_na_names, each = 3)

  all_or[[i]] <- as.data.frame(rbind(do.call("rbind", temp_or_base), do.call("rbind", temp_or_score)))
  colnames(all_or[[i]]) <- c("low_or", "or", "hi_or", "hi_cutoff", "lo_cutoff")
  all_or[[i]]$type <- rep(c("base", "score"), each = length(temp_or_base))
  all_or[[i]]$model <- rep(non_na_names, each = 3)

  all_diff[[i]] <- as.data.frame(do.call("rbind", temp_diff))
  colnames(all_diff[[i]]) <- c("score_diff", "base_diff", "cutoff")
  all_diff[[i]]$model <- rep(non_na_names, each = 3)
}

all_python_conc <- as.data.frame(do.call("rbind", all_python_conc))
all_python_conc$author <- all_author

names(all_prev) <- all_author
names(all_or) <- all_author
names(all_diff) <- all_author

saveRDS(list("prev" = all_prev, "or" = all_or, "diff" = all_diff, "python" = all_python_conc), "ready_to_plot/relative_risk.adjustment.RDS")

