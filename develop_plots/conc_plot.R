
####################################################################3
#age
print("age")

all_files <- list.files("../age_covar/res")
keep_conc <- list()
all_authors <- rep(NA, length(all_files))

for(i in 1:length(all_files)){
  all_authors[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../age_covar/res/", all_files[i]))

  vals <- data.frame(matrix(NA, nrow = length(res), ncol = 6))
  rownames(vals) <- names(res)
  colnames(vals) <- c("base_val", "score_val", "diff_val", "base_var", "score_var", "diff_var")  

  for(j in 1:length(res)){
    vals[j,] <- c(res[[j]][["base"]][["conc"]][[1]],
                  res[[j]][["score"]][["conc"]][[1]],
                  res[[j]][["score"]][["conc"]][[1]] - res[[j]][["base"]][["conc"]][[1]],
                  res[[j]][["base"]][["conc"]][[4]],
                  res[[j]][["score"]][["conc"]][[4]],
                  mean(c(res[[j]][["base"]][["conc"]][[4]], res[[j]][["score"]][["conc"]][[4]])))

  }
  keep_conc[[i]] <- vals
}

names(keep_conc) <- all_authors
saveRDS(keep_conc, "ready_to_plot/conc.age_covar.RDS")




####################################################################3
#censoring
print("censoring")

all_files <- list.files("../censoring/res")
keep_conc <- list()
all_authors <- rep(NA, length(all_files))

for(i in 1:length(all_files)){
  all_authors[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../censoring/res/", all_files[i]))

  vals <- data.frame(matrix(NA, nrow = length(res), ncol = 6))
  rownames(vals) <- names(res)
  colnames(vals) <- c("base_val", "score_val", "diff_val", "base_var", "score_var", "diff_var")

  for(j in 1:length(res)){
    vals[j,] <- c(res[[j]][["base"]][["conc"]][[1]],
                  res[[j]][["score"]][["conc"]][[1]],
                  res[[j]][["score"]][["conc"]][[1]] - res[[j]][["base"]][["conc"]][[1]],
                  res[[j]][["base"]][["conc"]][[4]],
                  res[[j]][["score"]][["conc"]][[4]],
                  mean(c(res[[j]][["base"]][["conc"]][[4]], res[[j]][["score"]][["conc"]][[4]])))

  }
  keep_conc[[i]] <- vals
}


names(keep_conc) <- all_authors
saveRDS(keep_conc, "ready_to_plot/conc.censoring.RDS")




####################################################################3
#competing risks
print("competing risks")

all_files <- list.files("../competing_risks/res")
keep_conc <- list()
all_authors <- rep(NA, length(all_files))

for(i in 1:length(all_files)){
  all_authors[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../competing_risks/res/", all_files[i]))


  if(is.null(res[[4]])){
    realj <- length(res)-1
    vals <- data.frame(matrix(NA, nrow = realj, ncol = 6))
    rownames(vals) <- names(res)[1:3]
    colnames(vals) <- c("base_val", "score_val", "diff_val", "base_var", "score_var", "diff_var")
  } else {
    realj <- 6
    vals <- data.frame(matrix(NA, nrow = realj, ncol = 6))
    rownames(vals) <- c(names(res)[1:3], paste0("icare_",names(res[[4]][[1]])))
    colnames(vals) <- c("base_val", "score_val", "diff_val", "base_var", "score_var", "diff_var")
  }

  for(j in 1:realj){
    if(j < 4){
    vals[j,] <- c(res[[j]][["base"]][["conc"]][6],
                  res[[j]][["score"]][["conc"]][6],
                  res[[j]][["score"]][["conc"]][6] - res[[j]][["base"]][["conc"]][6],
                  res[[j]][["base"]][["conc"]][7],
                  res[[j]][["score"]][["conc"]][7],
                  mean(c(res[[j]][["base"]][["conc"]][7], res[[j]][["score"]][["conc"]][7])))
    } else {
    vals[j,] <- c(res[[4]][["conc_base"]][[j-3]][6],
                  res[[4]][["conc_score"]][[j-3]][6],
                  res[[4]][["conc_score"]][[j-3]][6] - res[[4]][["conc_base"]][[j-3]][6],
                  res[[4]][["conc_base"]][[j-3]][7],
                  res[[4]][["conc_score"]][[j-3]][7],
                  mean(c(res[[4]][["conc_base"]][[j-3]][7], res[[4]][["conc_score"]][[j-3]][7])))
    }
  }
  keep_conc[[i]] <- vals
}


names(keep_conc) <- all_authors
saveRDS(keep_conc, "ready_to_plot/conc.competing_risks.RDS")






###################################################################
#disease labels
print("disease labels")


all_files <- list.files("../disease_labels/res", "RDS")
keep_conc <- list()
all_authors <- rep(NA, length(all_files))

for(i in 1:length(all_files)){
  all_authors[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../disease_labels/res/", all_files[i]))

  vals <- data.frame(matrix(NA, nrow = length(res)-1, ncol = 6))
  rownames(vals) <- names(res)[1:6]
  colnames(vals) <- c("base_val", "score_val", "diff_val", "base_var", "score_var", "diff_var")

  for(j in 1:6){
    vals[j,] <- c(res[[j]][["base"]][["conc"]][[1]],
                  res[[j]][["score"]][["conc"]][[1]],
                  res[[j]][["score"]][["conc"]][[1]] - res[[j]][["base"]][["conc"]][[1]],
                  res[[j]][["base"]][["conc"]][[4]],
                  res[[j]][["score"]][["conc"]][[4]],
                  mean(c(res[[j]][["base"]][["conc"]][[4]], res[[j]][["score"]][["conc"]][[4]])))

  }
  keep_conc[[i]] <- vals
}


names(keep_conc) <- all_authors
saveRDS(keep_conc, "ready_to_plot/conc.disease_labels.RDS")




###################################################################
#adjustment
print("adjustment")


all_files <- list.files("../adjustment/res", "RDS")
keep_conc <- list()
all_authors <- rep(NA, length(all_files))

for(i in 1:length(all_files)){
  all_authors[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../adjustment/res/", all_files[i]))

  vals <- data.frame(matrix(NA, nrow = length(res)-1, ncol = 6))
  rownames(vals) <- names(res)[1:7]
  colnames(vals) <- c("base_val", "score_val", "diff_val", "base_var", "score_var", "diff_var")

  for(j in 1:7){
    if(class(res[[j]]) == "list" ){
      if(!is.null(res[[j]][["score"]])){

    vals[j,] <- c(res[[j]][["base"]][["conc"]][[1]],
                  res[[j]][["score"]][["conc"]][[1]],
                  res[[j]][["score"]][["conc"]][[1]] - res[[j]][["base"]][["conc"]][[1]],
                  res[[j]][["base"]][["conc"]][[4]],
                  res[[j]][["score"]][["conc"]][[4]],
                  mean(c(res[[j]][["base"]][["conc"]][[4]], res[[j]][["score"]][["conc"]][[4]])))
      }
    }
  }
  keep_conc[[i]] <- vals
}


names(keep_conc) <- all_authors
saveRDS(keep_conc, "ready_to_plot/conc.adjustment.RDS")

