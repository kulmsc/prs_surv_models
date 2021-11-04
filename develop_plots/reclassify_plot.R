

get_reclass_rate <- function(orig, new, add_label, time_label, cut_offs = c(0.9, 0.95, 0.99)){
  ans <- rep(NA, length(cut_offs))
  for(i in 1:length(ans)){
    orig_inds <- which(orig[,ncol(orig)] > quantile(orig[,ncol(orig)], cut_offs[i]))
    new_inds <- which(new[,ncol(new)] > quantile(new[,ncol(new)], cut_offs[i]))
    ans[i] <- length(intersect(orig_inds, new_inds))/length(new_inds) #this is actually the stay the same rate, but will fix in plot script
  }
 
  names(ans) <- paste0("cut_", cut_offs)
  ans <- as.data.frame(t(ans))
  ans$"model" <- add_label  

  return(ans)
}



#################################################################
#################################################################
#################################################################
#age
all_files <- list.files("../age_covar/res")
over_author_list <- list()
all_authors <- rep(NA, length(all_files))

k <- 1
for(i in 1:length(all_files)){

  all_authors[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../age_covar/res/", all_files[i]))
  unchanged_score <- res[["unchanged"]][["score"]][["risk"]]
  unchanged_base <- res[["unchanged"]][["base"]][["risk"]]

  reclass_list <- list()
  jj <-1
  for(j in which(names(res) != "unchanged")){

    #want change in people classified at a few cutoffs and a few times
    temp_list <- list(get_reclass_rate(unchanged_base, res[[j]][["base"]][["risk"]], "base", res[[j]][["base"]][["time"]]),
                      get_reclass_rate(unchanged_score, res[[j]][["score"]][["risk"]], "score", res[[j]][["score"]][["time"]]),
                      get_reclass_rate(unchanged_score - unchanged_base, res[[j]][["score"]][["risk"]] - res[[j]][["base"]][["risk"]], "diff", res[[j]][["score"]][["time"]]))
    reclass_list[[jj]] <- do.call("rbind", temp_list)
    reclass_list[[jj]]$type <- names(res)[j]
    jj <- jj + 1
  }

  over_author_list[[i]] <- do.call("rbind", reclass_list)
}

names(over_author_list) <- all_authors
saveRDS(over_author_list, "ready_to_plot/reclass.age_covar.RDS")





#################################################################
#################################################################
#################################################################
#censoring
all_files <- list.files("../censoring/res")
over_author_list <- list()
all_authors <- rep(NA, length(all_files))

k <- 1
for(i in 1:length(all_files)){

  all_authors[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../censoring/res/", all_files[i]))
  unchanged_score <- res[["unchanged"]][["score"]][["risk"]]
  unchanged_base <- res[["unchanged"]][["base"]][["risk"]]

  reclass_list <- list()
  jj <- 1
  for(j in  which(names(res) != "unchanged")){

    #want change in people classified at a few cutoffs and a few times
    temp_list <- list(get_reclass_rate(unchanged_base, res[[j]][["base"]][["risk"]], "base", res[[j]][["base"]][["time"]]),
                      get_reclass_rate(unchanged_score, res[[j]][["score"]][["risk"]], "score", res[[j]][["score"]][["time"]]),
                      get_reclass_rate(unchanged_score - unchanged_base, res[[j]][["score"]][["risk"]] - res[[j]][["base"]][["risk"]], "diff", res[[j]][["score"]][["time"]]))
    reclass_list[[jj]] <- do.call("rbind", temp_list)
    reclass_list[[jj]]$type <- names(res)[j]
    jj <- jj + 1
  }

  over_author_list[[i]] <- do.call("rbind", reclass_list)
}

names(over_author_list) <- all_authors
saveRDS(over_author_list, "ready_to_plot/reclass.censoring.RDS")





#################################################################
#################################################################
#################################################################
#competing risks
all_files <- list.files("../competing_risks/res")
over_author_list <- list()
all_authors <- rep(NA, length(all_files))

k <- 1
for(i in 1:length(all_files)){

  all_authors[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
  res <- readRDS(paste0("../competing_risks/res/", all_files[i]))

  unchange_res <- readRDS(paste0("../age_covar/res/", all_authors[i], ".RDS"))
  unchanged_score <- unchange_res[["unchanged"]][["score"]][["risk"]]
  unchanged_base <- unchange_res[["unchanged"]][["base"]][["risk"]]

  reclass_list <- list()
  jj <- 1
  for(j in which(names(res) != "unchanged")){

    if(j < 4){
      temp_list <- list(get_reclass_rate(unchanged_base, res[[j]][["base"]][["risk"]], "base", res[[j]][["base"]][["time"]]),
                      get_reclass_rate(unchanged_score, res[[j]][["score"]][["risk"]], "score", res[[j]][["score"]][["time"]]),
                      get_reclass_rate(unchanged_score - unchanged_base, res[[j]][["score"]][["risk"]] - res[[j]][["base"]][["risk"]], "diff", res[[j]][["score"]][["time"]]))
      reclass_list[[jj]] <- do.call("rbind", temp_list)
      reclass_list[[jj]]$type <- names(res)[j]
      jj <- jj + 1
    } else {
      if(!is.null(res[[j]])){
        for(sj in 1:3){
          temp_list <- list(get_reclass_rate(unchanged_base, res[[j]][["base"]][[sj]], "base", 61:80),
                            get_reclass_rate(unchanged_score, res[[j]][["score"]][[sj]], "score", 61:80),
                            get_reclass_rate(unchanged_score - unchanged_base, res[[j]][["score"]][[sj]] - res[[j]][["base"]][[sj]], "diff", 61:80))
          reclass_list[[jj]] <- do.call("rbind", temp_list)
          reclass_list[[jj]]$type <- paste0(names(res[[4]][[2]])[sj], "_", names(res)[j])

          jj <- jj + 1
        }
      }
    }
  }

  over_author_list[[i]] <- do.call("rbind", reclass_list)
}

names(over_author_list) <- all_authors
saveRDS(over_author_list, "ready_to_plot/reclass.competing_risks.RDS")






#################################################################
#################################################################
#################################################################
#disease labels
#all_files <- list.files("../disease_labels/res", "RDS")
#over_author_list <- list()
#all_authors <- rep(NA, length(all_files))

#k <- 1
#for(i in 1:length(all_files)){

#  all_authors[i] <- strsplit(all_files[i], split = ".", fixed = T)[[1]][1]
#  res <- readRDS(paste0("../disease_labels/res/", all_files[i]))

#  unchange_res <- readRDS(paste0("../age_covar/res/", all_authors[i], ".RDS"))
#  unchanged_score <- unchange_res[["unchanged"]][["score"]][["risk"]]
#  unchanged_base <- unchange_res[["unchanged"]][["base"]][["risk"]]

#  reclass_list <- list()
#  jj <- 1
#  for(j in 1:6){

#    #want change in people classified at a few cutoffs and a few times
#    temp_list <- list(get_reclass_rate(unchanged_base, res[[j]][["base"]][["risk"]], "base", res[[j]][["base"]][["time"]]),
#                      get_reclass_rate(unchanged_score, res[[j]][["score"]][["risk"]], "score", res[[j]][["score"]][["time"]]),
#                      get_reclass_rate(unchanged_score - unchanged_base, res[[j]][["score"]][["risk"]] - res[[j]][["base"]][["risk"]], "diff", res[[j]][["score"]][["time"]]))
#    reclass_list[[jj]] <- do.call("rbind", temp_list)
#    reclass_list[[jj]]$type <- names(res)[j]
#    jj <- jj + 1
#  }

#  over_author_list[[i]] <- do.call("rbind", reclass_list)
#}

#names(over_author_list) <- all_authors
#saveRDS(over_author_list, "ready_to_plot/reclass.disease_labels.RDS")


