#python_acc <- readRDS(paste0("fuel_for_plot/python.", analysis_type, ".RDS"))
python_acc <- readRDS(paste0("fuel_for_plot/python.", "adjustment", ".RDS"))
#if(analysis_type == "adjustment"){
#  conc_limit <- 0
#}

if(length(res) == 18){
  new_res <- list()
  new_names <- c()
  
  j <- 1
  for(i in 1:length(res)){
    if(python_acc$concordance[python_acc$author == names(res)[i]] > conc_limit){
      new_res[[j]] <- res[[i]]
      new_names <- c(new_names, names(res)[i])
      j <- j + 1
    }
  }
  
  res <- new_res
  names(res) <- new_names

} else {
  #relative risk
  for(k in 1:3){
    new_res <- list()
    new_names <- c()
    
    j <- 1
    for(i in 1:length(res[[k]])){
      if(python_acc$concordance[python_acc$author == names(res[[k]])[i]] > conc_limit){
        new_res[[j]] <- res[[k]][[i]]
        new_names <- c(new_names, names(res[[k]])[i])
        j <- j + 1
      }
    }
    
    res[[k]] <- new_res
    names(res[[k]]) <- new_names 
    
  }
}
