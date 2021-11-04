library(survival)

#author <- "michailidou"
author = commandArgs(trailingOnly=TRUE)

#ideally want abs risk val for each year between assessment and now (or from 60 to 80)

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

train_surv_df$end_date <- as.Date(as.character(train_surv_df$end_date))
test_surv_df$end_date <- as.Date(as.character(test_surv_df$end_date))

train_surv_df$sqage <- train_surv_df$age^2
test_surv_df$sqage <- test_surv_df$age^2

train_surv_df$real_age <- train_surv_df$age
test_surv_df$real_age <- test_surv_df$age

if("sex" %in% colnames(surv_df)){
  use_vars <- c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "score", "sqage", "dob")
} else {
  use_vars <- c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "score", "sqage", "dob")
}
for(i in which(colnames(train_surv_df) %in% use_vars)){
  train_surv_df[,i] <- as.numeric(train_surv_df[,i])
  train_surv_df[,i] <- (train_surv_df[,i] - min(train_surv_df[,i]))/(max(train_surv_df[,i]) - min(train_surv_df[,i]))
}
for(i in which(colnames(test_surv_df) %in% use_vars)){
  test_surv_df[,i] <- as.numeric(test_surv_df[,i])
  test_surv_df[,i] <- (test_surv_df[,i] - min(test_surv_df[,i]))/(max(test_surv_df[,i]) - min(test_surv_df[,i]))
}


if("sex" %in% colnames(surv_df)){
  base_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")

  base_sq_covars <- c("age + sex + sqage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_sq_covars <- c("age + sex + sqage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")

  base_int_covars <- c("age + sex + sqage + age:sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_int_covars <- c("age + sex + sqage + age:sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")

  base_dob_covars <- c("age + dob + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_dob_covars <- c("age + dob + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")

} else {
  base_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")

  base_sq_covars <- c("age + sqage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_sq_covars <- c("age + sqage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")

  base_int_covars <- c("age + sqage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_int_covars <- c("age + sqage +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")

  base_dob_covars <- c("dob + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_dob_covars <- c("dob +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")
}

all_base_covars <- list(base_covars, base_sq_covars, base_int_covars, base_dob_covars)
all_score_covars <- list(score_covars, score_sq_covars, score_int_covars, score_dob_covars)
all_name_covars <- c("unchanged", "square", "interaction", "dob")

#######################################33



get_inds <- function(real_times, want_times){
  inds <- rep(NA, length(want_times))
  for(i in 1:length(want_times)){
    inds[i] <- which.min(abs(real_times - want_times[i]))
  }
  return(inds)
}

get_survfit <- function(fit_model, new_data, desire_time, input_model = "ms"){
  all_res <- list()
  counter <- 1
  interval <- 5000
  for(ind in seq(1, nrow(new_data), interval)){
    print(ind)
    fit_prod <- survfit(fit_model, newdata = new_data[ind:(ind+interval-1),])
    fit_prod <- cbind(fit_prod$time, fit_prod$cumhaz)
    all_res[[counter]] <- fit_prod[get_inds(fit_prod[,1], desire_time),]
    counter <- counter + 1
  }
  return(all_res) #rows are the time, cols are the people
}


#age refers to age at time of assessment
#can do no age, normal age, +age^2, + age^2 + age:sex
#year of birth instead of age
#age match to same proportion as cases

all_res <- list()
################################## UNCHANGED #############################################
for(i in 1:length(all_base_covars)){

  cox_model_base <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", all_base_covars[[i]])), data = train_surv_df)

  base_conc_obj <- concordance(cox_model_base, newdata=test_surv_df)

  cox_fit_base <- get_survfit(cox_model_base, test_surv_df, seq(365, max(test_surv_df$time), 365))
  cox_fit_base <- t(do.call("cbind", cox_fit_base)) #rows are the individuals, cols are the time
  cox_fit_base <- list("coef" = summary(cox_model_base)$coef, "time" = cox_fit_base[1,], "risk" = cox_fit_base[cox_fit_base[,1] != cox_fit_base[1,1],], "conc" = base_conc_obj)



  cox_model_score <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", all_score_covars[[i]])), data = train_surv_df)

  score_conc_obj <- concordance(cox_model_score, newdata=test_surv_df)

  cox_fit_score <- get_survfit(cox_model_score, test_surv_df, seq(365, max(test_surv_df$time), 365))
  cox_fit_score <- t(do.call("cbind", cox_fit_score)) #rows are the individuals, cols are the time
  cox_fit_score <- list("coef" = summary(cox_model_score)$coef, "time" = cox_fit_score[1,], "risk" = cox_fit_score[cox_fit_score[,1] != cox_fit_score[1,1],], "conc" = score_conc_obj)


  all_res[[i]] <- list("base" = cox_fit_base, "score" = cox_fit_score)

}

names(all_res) <- all_name_covars



######################################################3

train_surv_df$age <- train_surv_df$real_age/365
test_surv_df$age <- test_surv_df$real_age/365

x <- table(round(train_surv_df$age)[train_surv_df$pheno == 1])
y <- table(round(train_surv_df$age)[train_surv_df$pheno == 0])
y <- y[names(y) %in% names(x)]
x <- x[names(x) %in% names(y)]
borrow_ratio <- floor(min(y/x))

parts_train_df <- list()
for(i in 1:length(x)){
  control_pick <- sample(which(round(train_surv_df$age) == names(x)[i] & train_surv_df$pheno == 0), x[i]*borrow_ratio)
  case_pick <- which(round(train_surv_df$age) == names(x)[i] & train_surv_df$pheno == 1)
  if(length(case_pick) > 1){
    case_pick <- sample(case_pick, x[i])
  }
  picked_inds <- sort(c(control_pick, case_pick))
  parts_train_df[[i]] <- train_surv_df[picked_inds,]
}
use_train_surv_df <- do.call("rbind", parts_train_df)



cox_model_base <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", all_base_covars[[1]])), data = use_train_surv_df)

base_conc_obj <- concordance(cox_model_base, newdata=test_surv_df)

cox_fit_base <- get_survfit(cox_model_base, test_surv_df, seq(365, max(test_surv_df$time), 365))
cox_fit_base <- t(do.call("cbind", cox_fit_base)) #rows are the individuals, cols are the time
cox_fit_base <- list("coef" = summary(cox_model_base)$coef, "time" = cox_fit_base[1,], "risk" = cox_fit_base[cox_fit_base[,1] != cox_fit_base[1,1],], "conc" = base_conc_obj)



cox_model_score <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", all_score_covars[[1]])), data = use_train_surv_df)

score_conc_obj <- concordance(cox_model_score, newdata=test_surv_df)

cox_fit_score <- get_survfit(cox_model_score, test_surv_df, seq(365, max(test_surv_df$time), 365))
cox_fit_score <- t(do.call("cbind", cox_fit_score)) #rows are the individuals, cols are the time
cox_fit_score <- list("coef" = summary(cox_model_score)$coef, "time" = cox_fit_score[1,], "risk" = cox_fit_score[cox_fit_score[,1] != cox_fit_score[1,1],], "conc" = score_conc_obj)


all_res[[5]] <- list("base" = cox_fit_base, "score" = cox_fit_score)

names(all_res) <- c(all_name_covars, "age_match")

saveRDS(all_res, paste0("res/", author, ".RDS"))
