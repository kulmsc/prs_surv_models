

get_km_df <- function(all_risk, time_vec, author, type){
  high_risk_ind <- which(all_risk[,ncol(all_risk)] >= quantile(all_risk[,ncol(all_risk)], 0.8))
  inter_risk_ind <- which(all_risk[,ncol(all_risk)] < quantile(all_risk[,ncol(all_risk)], 0.8) & all_risk[,ncol(all_risk)] > quantile(all_risk[,ncol(all_risk)], 0.2))
  low_risk_ind <- which(all_risk[,ncol(all_risk)] <= quantile(all_risk[,ncol(all_risk)], 0.2))

  temp_list <- list(rep(time_vec, 3),
                    rep(c("high", "inter", "low"), each = length(time_vec)),
                    c(apply(all_risk[high_risk_ind,], 2, mean),
                    apply(all_risk[inter_risk_ind,], 2, mean),
                    apply(all_risk[low_risk_ind,], 2, mean)),
                    c(apply(all_risk[high_risk_ind,], 2, sd),
                    apply(all_risk[inter_risk_ind,], 2, sd),
                    apply(all_risk[low_risk_ind,], 2, sd)))
  comp_risk <- as.data.frame(do.call("cbind", temp_list))
  colnames(comp_risk) <- c("time", "risk_group", "val", "sd")
  comp_risk$time <- as.numeric(as.character(comp_risk$time))
  comp_risk$val <- as.numeric(as.character(comp_risk$val))
  comp_risk$sd <- as.numeric(as.character(comp_risk$sd))
  comp_risk$risk_group <- as.character(comp_risk$risk_group)
  comp_risk$author <- author
  comp_risk$model <- type
  return(comp_risk)
}




#################################################################
#################################################################
#################################################################
#age
print("age")

all_files <- list.files("../age_covar/res")
all_author <- rep(NA, length(all_files))
keep_coef <- list()


k <- 1
for(i in 1:length(all_files)){
  all_author[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../age_covar/res/", all_files[i]))
  temp_list <- list()
  for(j in 1:length(res)){

    score_risk <- get_km_df(res[[j]][["score"]][["risk"]], res[[j]][["score"]][["time"]], all_author[i], "score")
    base_risk <- get_km_df(res[[j]][["base"]][["risk"]], res[[j]][["base"]][["time"]], all_author[i], "base")

    temp_list[[j]] <- rbind(score_risk, base_risk)
    temp_list[[j]]$type <- names(res)[j]
  }
  keep_coef[[i]] <- do.call("rbind", temp_list)
}

names(keep_coef) <- all_author
saveRDS(keep_coef, "ready_to_plot/km.age_covar.RDS")


#########################################################################
#########################################################################
#########################################################################
#censoring
print("censoring")

all_files <- list.files("../censoring/res")
all_author <- rep(NA, length(all_files))
keep_coef <- list()

k <- 1
for(i in 1:length(all_files)){
  all_author[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../censoring/res/", all_files[i]))
  temp_list <- list()
  for(j in 1:length(res)){

    score_risk <- get_km_df(res[[j]][["score"]][["risk"]], res[[j]][["score"]][["time"]], all_author[i], "score")
    base_risk <- get_km_df(res[[j]][["base"]][["risk"]], res[[j]][["base"]][["time"]], all_author[i], "base")    

    temp_list[[j]] <- rbind(score_risk, base_risk)
    temp_list[[j]]$type <- names(res)[j]
  }
  keep_coef[[i]] <- do.call("rbind", temp_list)
}

names(keep_coef) <- all_author
saveRDS(keep_coef, "ready_to_plot/km.censoring.RDS")




#########################################################################
#########################################################################
#########################################################################
#competing risks
print("competing risks")

all_files <- list.files("../competing_risks/res")
all_author <- rep(NA, length(all_files))
keep_coef <- list()
extra_names <- c("icare_gbd", "icare_holt", "icare_combo")

k <- 1
for(i in 1:length(all_files)){
  all_author[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../competing_risks/res/", all_files[i]))
  temp_list <- list()
  for(j in 1:length(res)){
    if(!is.null(res[[j]])){

      if(names(res)[j] == "icare"){
        for(sj in 1:3){
          score_risk <- get_km_df(res[[j]][["score"]][[sj]], res[["ms"]][["score"]][["time"]], all_author[i], "score")
          base_risk <- get_km_df(res[[j]][["base"]][[sj]], res[["ms"]][["score"]][["time"]], all_author[i], "base")

          temp_list[[j]] <- rbind(score_risk, base_risk)
          temp_list[[j]]$type <- paste0(names(res)[j], "_", names(res[[j]][["score"]])[sj])
          k <- k + 1
        }
      } else {
        score_risk <- get_km_df(res[[j]][["score"]][["risk"]], res[[j]][["score"]][["time"]], all_author[i], "score")
        base_risk <- get_km_df(res[[j]][["base"]][["risk"]], res[[j]][["base"]][["time"]], all_author[i], "base")

        temp_list[[j]] <- rbind(score_risk, base_risk)
        temp_list[[j]]$type <- names(res)[j]
        k <- k + 1
      }
    } 
  }
  keep_coef[[i]] <- do.call("rbind", temp_list)
}

names(keep_coef) <- all_author
saveRDS(keep_coef, "ready_to_plot/km.competing_risks.RDS")





#########################################################################
#########################################################################
#########################################################################
#disease_labels
print("disease labels")

all_files <- list.files("../disease_labels/res", "RDS")
all_author <- rep(NA, length(all_files))
keep_coef <- list()

k <- 1
for(i in 1:length(all_files)){
  all_author[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../disease_labels/res/", all_files[i]))
  temp_list <- list()
  for(j in 1:6){

    score_risk <- get_km_df(res[[j]][["score"]][["risk"]], res[[j]][["score"]][["time"]], all_author[i], "score")
    base_risk <- get_km_df(res[[j]][["base"]][["risk"]], res[[j]][["base"]][["time"]], all_author[i], "base")

    temp_list[[j]] <- rbind(score_risk, base_risk)
    temp_list[[j]]$type <- names(res)[j]
  }
  keep_coef[[i]] <- do.call("rbind", temp_list)
}

names(keep_coef) <- all_author
saveRDS(keep_coef, "ready_to_plot/km.disease_labels.RDS")





#########################################################################
#########################################################################
#########################################################################
#adjustment
print("adjustment")

all_files <- list.files("../adjustment/res", "RDS")
all_author <- rep(NA, length(all_files))
keep_coef <- list()

k <- 1
for(i in 1:length(all_files)){
  all_author[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../adjustment/res/", all_files[i]))
  temp_list <- list()
  for(j in 1:7){
    if(class(res[[j]]) == "list"){
      if(!is.null(res[[j]][["score"]])){
        score_risk <- get_km_df(res[[j]][["score"]][["risk"]], res[[j]][["score"]][["time"]], all_author[i], "score")
        base_risk <- get_km_df(res[[j]][["base"]][["risk"]], res[[j]][["base"]][["time"]], all_author[i], "base")

        temp_list[[j]] <- rbind(score_risk, base_risk)
        temp_list[[j]]$type <- names(res)[j]
      } else {
        temp_list[[j]] <- NULL
      }
    } else {
      temp_list[[j]] <- NULL
    }
  }
  keep_coef[[i]] <- do.call("rbind", temp_list)
}

names(keep_coef) <- all_author
saveRDS(keep_coef, "ready_to_plot/km.adjustment.RDS")


