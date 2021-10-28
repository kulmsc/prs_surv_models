

check_model_inputs <- function(apply.cov.profile , model.log.RR, model.ref.dataset, model.cov.info, model.formula){
  return_and_print <- print
  if(!is.data.frame(apply.cov.profile ) && !is.matrix(apply.cov.profile )){
    return_and_print("ERROR: apply.cov.profile  must be a dataframe.")
    stop()
  }
  if(!is.data.frame(model.ref.dataset)&& !is.matrix(model.ref.dataset)){
    return_and_print("ERROR: model.ref.dataset must be a dataframe (with appropriate column names).")
    stop()
  }
  if(sum(rowSums(is.na(model.ref.dataset)))>0){
    return_and_print("ERROR: model.ref.dataset cannot contain NAs.")
    stop()
  }
  check_triple_check(model.cov.info)
  model.cov.info = process_triple_check(model.cov.info)

  ### Will verify all names and proper order of things against model.cov.info - so check against model first
  if(is.null(colnames(model.cov.info))){
    return_and_print("ERROR: model.cov.info must have same names and order as predictors in model.formula.")
    stop()
  }
  if(length(colnames(model.cov.info))!=length(all.vars(model.formula)[2:length(all.vars(model.formula))])){
    return_and_print("ERROR: model.cov.info must have same names and order as predictors in model.formula.")
    stop()
  }
  if( sum(colnames(model.cov.info)!=all.vars(model.formula)[2:length(all.vars(model.formula))])>0 ){
    return_and_print("ERROR: model.cov.info must have same names and order as predictors in model.formula.")
    stop()
  }
  if(prod(colnames(model.cov.info)==colnames(apply.cov.profile ))!=1){
    return_and_print("ERROR: Variable names of model.cov.info should match column names of apply.cov.profile .")
    stop()
  }
  for(k in 1:ncol(model.cov.info)){
    if(is.factor(model.cov.info[,k])){
      if(sum(is.element(apply.cov.profile [,k], model.cov.info[,k])+(is.na(c(apply.cov.profile [,k])))==0)>0){
        return_and_print(paste("ERROR: apply.cov.profile  categorical Variable (", colnames(model.cov.info)[k] ,") has levels outside levels modeled."))
        stop()
      }
    }
  }
  if(prod(colnames(model.cov.info)==colnames(model.ref.dataset))!=1){
    return_and_print("ERROR: Variable names of model.cov.info should match column names of model.ref.dataset.")
    stop()
  }
  if(nrow(model.ref.dataset) < 200){
    print(paste("WARNING: Samples in referent distribution model.ref.dataset should be large.  Currently only size ", nrow(model.ref.dataset),".", sep=""))
  }
  if(is.null(model.ref.dataset.weights)){
    model.ref.dataset.weights = rep(1/nrow(model.ref.dataset), nrow(model.ref.dataset))
  }else{
    if(length(model.ref.dataset.weights)!=nrow(model.ref.dataset)){
      return_and_print("ERROR: If model.ref.dataset.weights is provided it must be same length as the number of rows in model.ref.dataset.")
      stop()
    }
    if( sum((model.ref.dataset.weights<0), na.rm=TRUE)!=0){
      return_and_print("ERROR: model.ref.dataset.weights must not contain negative values.")
      stop()
    }
    model.ref.dataset.weights = model.ref.dataset.weights / sum(model.ref.dataset.weights, na.rm=TRUE)
  }
  ## for now don't allow NA's here
  for(k in 1:ncol(model.cov.info)){
    if(is.factor(model.cov.info[,k])){
      if(sum(is.element(c(model.ref.dataset[,k]), model.cov.info[,k])==0)>0){
        return_and_print(paste("ERROR: model.ref.dataset categorical Variable (", colnames(model.cov.info)[k] ,") has levels outside levels modeled."))
        stop()
      }
    }
  }
  variables <- unique(all.vars(model.formula))[-1]
  data_use  <- subset(model.cov.info, select = variables)

  if(is.vector(model.log.RR)==FALSE & is.matrix(model.log.RR)==FALSE){
    return_and_print("ERROR: model.log.RR must be either a vector, or a matrix with only one column.")
    stop()
  }
  if(is.matrix(model.log.RR)==TRUE){
    if(dim(model.log.RR)[2]>1){
      return_and_print("ERROR: model.log.RR must be either a vector, or a matrix with only one column.")
      stop()
    }
  }
  if(is.vector(model.log.RR)==TRUE){
    temp = names(model.log.RR)
    model.log.RR = cbind(model.log.RR)
    rownames(model.log.RR) = temp
  }
  res = list(); res[[1]] = data_use; res[[2]] = model.log.RR; res[[3]] = model.ref.dataset.weights; res
}

check_design_matrix <- function(model.log.RR, design_covs){

    if(is.null(rownames(model.log.RR)) ){
      return_and_print("ERROR: row names of model.log.RR must match names and order in design matrix.")
      print(paste("Row names of model.log.RR = NULL"))
      print(paste("Names of design_matrix = ", paste0(colnames(design_covs), collapse=", ") ))
      stop()
    }
    if( sum( rownames(model.log.RR)!= colnames(design_covs) )>0){
      return_and_print("ERROR: row names of model.log.RR must match names and order in design matrix.")
      print(cbind(names(model.log.RR), colnames(design_covs)))
      stop()
    }
}




check_triple_check <- function(model.cov.info){

  if( is.list(model.cov.info)!=TRUE){
    return_and_print("ERROR: model.cov.info must be a list.")
    stop()
  }
  for(i in 1:length(model.cov.info)){
    if( is.list(model.cov.info[[i]])!=TRUE ){
      return_and_print("ERROR: Each element of model.cov.info must be a list (that contains the information for a variable). ")
      stop()
    }
    if( is.null(model.cov.info[[i]]$name) ){
      return_and_print("ERROR: Each list within model.cov.info describing a variable must have a 'name' field. ")
      stop()
    }
    if( is.null(model.cov.info[[i]]$type) ){
      return_and_print("ERROR: Each list within model.cov.info describing a variable must have a 'type' field. ")
      stop()
    }}
  for(i in 1:length(model.cov.info)){

    if( model.cov.info[[i]]$type == "continuous"  & is.null(model.cov.info[[i]]$levels)!=TRUE){
      print(paste("WARNING: You have specified that the variable '", model.cov.info[[i]]$name, "' is continuous, so the 'levels' input for that variable is not needed and will not be used.", sep="") )
    }# allow to continue
    if( model.cov.info[[i]]$type == "factor"  & is.null(model.cov.info[[i]]$levels)==TRUE){
      return_and_print(paste("ERROR: You must specify the 'levels' field of the variable '", model.cov.info[[i]]$name, "', because you gave its type to be 'factor'.", sep="") )
      stop()
    }
    if( model.cov.info[[i]]$type == "factor"  & is.null(model.cov.info[[i]]$levels)==FALSE){
      if(is.vector(model.cov.info[[i]]$levels)!=TRUE ){
        return_and_print(paste("ERROR: You must specify the 'levels' field of the variable '", model.cov.info[[i]]$name, "', as a vector.", sep="") )
        stop()
      }
      if(is.null(model.cov.info[[i]]$ref)==FALSE){
        if(is.element(model.cov.info[[i]]$ref, model.cov.info[[i]]$levels)==FALSE){
          return_and_print(paste("ERROR: The 'ref' field of the variable '", model.cov.info[[i]]$name, "' must specify a referent level that is contained in the 'levels' field for that variable.", sep="") )
          stop()
        }
      }
    }
  }
}

check_disease_rates <- function(filename){

  lambda = read.table(filename, sep=",", header=TRUE)
  lambda = check_flexible_rate_inputs(lambda, "model.disease.incidence.rates")
  lambda
}

check_competing_rates <- function(filename, lambda){

  if(is.null(filename) || filename==""){
    model.competing.incidence.rates= data.frame( cbind(lambda[,1], rep(0, length(lambda[,1]))) )
  }else{
    model.competing.incidence.rates = read.table(filename, sep=",", header=TRUE)
  }
  model.competing.incidence.rates = check_flexible_rate_inputs(model.competing.incidence.rates, "model.competing.incidence.rates")
  model.competing.incidence.rates
}

## function to convert (start, end, rate) matrix into (integer_ages, rate) matrix
format_flexible_rate_inputs<-function(mat){

  if(ncol(mat)==3){
    start <- mat[,1]
    end   <- mat[,2]
    rate  <- mat[,3]

    integer_ages <- seq(min(start), max(end))
    formatted    <- cbind( integer_ages, rep(NA, length(integer_ages))  )
    for(i in 1:nrow(mat)){
      these              <- which(formatted[,1]>=start[i] & formatted[,1]<=end[i])
      formatted[these,2] <- rate[i]/length(these)
    }
    colnames(formatted) <- c("ages", "rates")
    formatted
    }else{
      formatted = mat
    }
  formatted
}

check_flexible_rate_inputs<-function(mat, name){
  if(!is.data.frame(mat) && !is.matrix(mat)){
    return_and_print(paste("ERROR: ", name, " must be provided as a matrix.", sep=""))
    stop()
  }
  if(ncol(mat)!=3 & ncol(mat)!=2){
    return_and_print(paste("ERROR: ", name, " must have either 2 columns: [Ages,Rates] or 3 columns:[Start_Ages, End_Ages, Rates].", sep=""))
    stop()
  }
  ll = nrow(mat)
  if(ll>1 & ncol(mat)==3){
    if(sum(mat[2:ll,1] - mat[1:(ll-1),2])!=0){
      return_and_print(paste("ERROR: The rates provided in ", name, " must cover sequential age intervals (i.e. if an interval ends at age 30, the next interval must start at age 31).", sep=""))
      stop()
    }
  }
  if(ncol(mat)==2 && sum(mat[,1]%%1)!=0){
    return_and_print(paste("ERROR: The first column of ", name, " should be integer ages.", sep=""))
    stop()
  }
  if(sum(mat[,ncol(mat)]<0, na.rm=TRUE) + sum(mat[,ncol(mat)]>1, na.rm=TRUE)!=0){
    return_and_print("ERROR: The rates should be probabilities between 0 and 1.")
    stop()
  }
  format_flexible_rate_inputs(mat)
}

check_rates <- function(model.competing.incidence.rates, lambda, apply.age.start, apply.age.interval.length){

  lambda = check_flexible_rate_inputs(lambda, "model.disease.incidence.rates")

  if(is.null(model.competing.incidence.rates)){  model.competing.incidence.rates= cbind(lambda[,1], rep(0, length(lambda[,1])))   }

  model.competing.incidence.rates = check_flexible_rate_inputs(model.competing.incidence.rates, "model.competing.incidence.rates")

  if( prod(is.element(seq(range(apply.age.start)[1], range(apply.age.start + apply.age.interval.length)[2]), lambda[,1])) == 0){
    return_and_print("ERROR: The 'model.disease.incidence.rates' input must have age-specific rates for each integer age covered by the prediction intervals defined by 'apply.age.start' and 'apply.age.interval.length.'  You must make these inputs consistent with one another to proceed.")
    stop()
  }
  if( prod(is.element(seq(range(apply.age.start)[1], range(apply.age.start + apply.age.interval.length)[2]), model.competing.incidence.rates[,1])) == 0){
    return_and_print("ERROR: The 'model.competing.incidence.rates' input must have age-specific rates for each integer age covered by the prediction intervals defined by 'apply.age.start' and 'apply.age.interval.length.'  You must make these inputs consistent with one another to proceed.")
    stop()
  }
  res=list(); res[[1]] = lambda; res[[2]] = model.competing.incidence.rates; res
}




# transform model.cov.info into a form easily usable in the function
process_triple_check <- function(model.cov.info){
  maximum_dim = 0
  for(i in 1:length(model.cov.info)){
    maximum_dim = max(maximum_dim, length(model.cov.info[[i]]$levels) )
  }
  matt = data.frame(starter_column = rep(0, maximum_dim))

  for(i in 1:length(model.cov.info)){

    if(model.cov.info[[i]]$type == "continuous"){
        matt[,model.cov.info[[i]]$name] = rep(0,maximum_dim)
    }
    if(model.cov.info[[i]]$type == "factor"){
        matt[,model.cov.info[[i]]$name] = factor( c(model.cov.info[[i]]$levels, rep(model.cov.info[[i]]$levels[1], maximum_dim - length(model.cov.info[[i]]$levels))), levels = unique(c(model.cov.info[[i]]$levels, rep(model.cov.info[[i]]$levels[1], maximum_dim - length(model.cov.info[[i]]$levels)))))
        ### if there is a referent level specified, recode that way
        if(is.null(model.cov.info[[i]]$ref)==FALSE){
            matt[,model.cov.info[[i]]$name] = relevel(matt[,model.cov.info[[i]]$name], ref = as.character(model.cov.info[[i]]$ref))
        }
    }
  }
  matt[,1] <- NULL
  matt
}
