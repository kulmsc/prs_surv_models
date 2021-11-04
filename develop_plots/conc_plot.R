
####################################################################3
#age

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

#conc can only be calculated with a normal coxph model
#could try to use the linear predictors of the compete risk models to form a dummy coxph
#that then goes into the concordance function, but would have to redo everything
#still unclear if correct
