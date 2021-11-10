


####################################################################3
#age
print("age")

all_files <- list.files("../age_covar/res")
keep_coef <- list()

k <- 1
for(i in 1:length(all_files)){
  res <- readRDS(paste0("../age_covar/res/", all_files[i]))
  for(j in 1:length(res)){
    all_coef <- as.data.frame(res[[j]][["score"]][["coef"]])
    all_coef$author <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
    all_coef$model <- names(res)[j]
    keep_coef[[k]] <- all_coef[rownames(all_coef) == "score",]
    k <- k + 1
  }
}

keep_coef <- as.data.frame(do.call("rbind", keep_coef))
saveRDS(keep_coef, "ready_to_plot/forest.age_covar.RDS")




######################################################################
######################################################################
######################################################################
#censoring
print("censoring")


all_files <- list.files("../censoring/res")
keep_coef <- list()

k <- 1
for(i in 1:length(all_files)){
  res <- readRDS(paste0("../censoring/res/", all_files[i]))
  for(j in 1:length(res)){
    all_coef <- as.data.frame(res[[j]][["score"]][["coef"]])
    all_coef$author <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
    all_coef$model <- names(res)[j]
    keep_coef[[k]] <- all_coef[rownames(all_coef) == "score",]
    k <- k + 1
  }
}

keep_coef <- as.data.frame(do.call("rbind", keep_coef))
saveRDS(keep_coef, "ready_to_plot/forest.censoring.RDS")


######################################################################
######################################################################
######################################################################
#competing risks
print("competing risks")


all_files <- list.files("../competing_risks/res")
keep_coef <- list()

k <- 1
for(i in 1:length(all_files)){
  res <- readRDS(paste0("../competing_risks/res/", all_files[i]))
  for(j in 1:3){
      all_coef <- as.data.frame(t(res[[j]][["score"]][["coef"]]))
      all_coef$author <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
      all_coef$model <- names(res)[j]
      all_coef$type <- "by_time"

      #if("sub_score" %in% names(res[[j]])){
      #  sub_coef <- as.data.frame(t(res[[j]][["sub_score"]][["coef"]]))
      #  sub_coef$author <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
      #  sub_coef$model <- names(res)[j]
      #  sub_coef$type <- "by_age"
      #  all_coef <- rbind(all_coef, sub_coef)
      #}

      #keep_coef[[k]] <- all_coef[rownames(all_coef) == "score",]
      if("robust se" %in% colnames(all_coef)){
        all_coef <- all_coef[,colnames(all_coef) != "robust se"]
      }
      colnames(all_coef) <- c("coef", "exp_coef", "se", "z", "p", "author", "model", "type")

      keep_coef[[k]] <- all_coef
      k <- k + 1
  }
}


keep_coef <- as.data.frame(do.call("rbind", keep_coef))
saveRDS(keep_coef, "ready_to_plot/forest.competing_risks.RDS")



######################################################################
######################################################################
######################################################################
#disease labels
print("disease labels")

all_files <- list.files("../disease_labels/res", "RDS")
keep_coef <- list()

k <- 1
for(i in 1:length(all_files)){
  res <- readRDS(paste0("../disease_labels/res/", all_files[i]))
  for(j in 1:6){
    all_coef <- as.data.frame(res[[j]][["score"]][["coef"]])
    all_coef$author <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
    all_coef$model <- names(res)[j]
    keep_coef[[k]] <- all_coef[rownames(all_coef) == "score",]
    k <- k + 1
  }
}

keep_coef <- as.data.frame(do.call("rbind", keep_coef))
saveRDS(keep_coef, "ready_to_plot/forest.disease_labels.RDS")




######################################################################
######################################################################
######################################################################
#adjustment
print("adjustment")


all_files <- list.files("../adjustment/res", "RDS")
keep_coef <- list()

k <- 1
for(i in 1:length(all_files)){
  res <- readRDS(paste0("../adjustment/res/", all_files[i]))
  for(j in 1:7){
    if(class(res[[j]]) == "list" ){
      if(!is.null(res[[j]][["score"]])){
        all_coef <- as.data.frame(res[[j]][["score"]][["coef"]])
        all_coef$author <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
        all_coef$model <- names(res)[j]
        keep_coef[[k]] <- all_coef[rownames(all_coef) == "score",]
      }
    }
    k <- k + 1
  }
}

keep_coef <- as.data.frame(do.call("rbind", keep_coef))
saveRDS(keep_coef, "ready_to_plot/forest.adjustment.RDS")


