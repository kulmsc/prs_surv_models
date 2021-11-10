library(survival)

#want to remove age as a covariate from left censor
#can try to write a function to pull risk from the same time as in right and left censor




#author <- "christophersen"
author = commandArgs(trailingOnly=TRUE)

#ideally want abs risk val for each year between assessment and now (or from 60 to 80)

#read in the UKBB data
surv_df <- readRDS(paste0("../init_data/surv_data/survdf.", author, ".RDS"))
surv_df$time <- as.numeric(as.Date(surv_df$end_date) - as.Date(surv_df$date_attend)) #time is counted from the date of assessment

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


if("sex" %in% colnames(surv_df)){
  base_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")
  left_base_covars <- c("sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  left_score_covars <- c("sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")

} else {
  base_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")
  left_base_covars <- c("PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  left_score_covars <- c("PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")

}



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


faster_fit <- function(fit_model, new_data, desire_time, model_type = "cox"){
  if(model_type == "ms"){
    coefs <- fit_model$coef[grep("1:3", names(fit_model$coef))]
    coef_names <- unlist(lapply(strsplit(names(coefs), "_"), function(x) x[1]))
  } else {
    coefs <- fit_model$coef
    coef_names <- names(fit_model$coef)
  }
  new_data <- new_data[,colnames(new_data) %in% coef_names]
  new_data <- new_data[,order(colnames(new_data))]
  coefs <- coefs[order(coef_names)]

  simple_preds <- colSums(t(new_data) * coefs)

  fit_prod <- survfit(fit_model, newdata = new_data[1:1000,], se.fit = FALSE)
  if(model_type == "ms"){
    fit_prod <- cbind(fit_prod$time, fit_prod$cumhaz[,,2])
    fit_prod <- fit_prod[get_inds(fit_prod[,1], desire_time),]
  } else if(model_type == "left"){
    #fit_prod <- cbind(fit_prod$time, fit_prod$cumhaz)
    print("LEFT")
  } else {
    fit_prod <- cbind(fit_prod$time, fit_prod$cumhaz)
    fit_prod <- fit_prod[get_inds(fit_prod[,1], desire_time),]
  }


  if(model_type == "left"){

    full_preds <- rbind(do.call("rbind", lapply(1:1000, function(x) fit_prod$cumhaz[get_inds(fit_prod$time, desire_time[x,]),x])),
                        matrix(0, length(simple_preds)-1000, ncol = ncol(desire_time)))
    #get inds for everyone
    mat_inds <- rbind(matrix(-10, nrow = 1000, ncol = ncol(desire_time)),
                       do.call("rbind", lapply(1001:nrow(full_preds), function(x) get_inds(fit_prod$time, desire_time[x,]))))
    for(i in 1:length(fit_prod$time)){
      mod <- lm(fit_prod$cumhaz[i,] ~ exp(simple_preds[1:1000]))
      temp_ans <- mod$coef[1] + exp(simple_preds)*mod$coef[2]
      for(j in 1:ncol(full_preds)){
        full_preds[mat_inds[,j] == i, j] <- temp_ans[mat_inds[,j] == i]
      }
    }
    desire_time <- 1:ncol(full_preds)

  } else {
    full_preds <- matrix(0, nrow = length(simple_preds), ncol = nrow(fit_prod))
    for(i in 1:ncol(full_preds)){
      mod <- lm(fit_prod[i,-1] ~ exp(simple_preds[1:1000]))
      full_preds[,i] <- mod$coef[1] + exp(simple_preds)*mod$coef[2]
    }
  }




  return(rbind(desire_time, full_preds))
}





################################## UNCHANGED #############################################

cox_model_base <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", base_covars)), data = train_surv_df)

base_conc_obj <- concordance(cox_model_base, newdata=test_surv_df)

#cox_fit_base <- get_survfit(cox_model_base, test_surv_df, seq(365, max(test_surv_df$time), 365))
cox_fit_base <- faster_fit(cox_model_base, test_surv_df, seq(365, max(test_surv_df$time), 365))
#cox_fit_base <- t(do.call("cbind", cox_fit_base)) #rows are the individuals, cols are the time
cox_fit_base <- list("coef" = summary(cox_model_base)$coef, "time" = cox_fit_base[1,], "risk" = cox_fit_base[cox_fit_base[,1] != cox_fit_base[1,1],], "conc" = base_conc_obj)




cox_model_score <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", score_covars)), data = train_surv_df)

score_conc_obj <- concordance(cox_model_score, newdata=test_surv_df)

#cox_fit_score <- get_survfit(cox_model_score, test_surv_df, seq(365, max(test_surv_df$time), 365))
cox_fit_score <- faster_fit(cox_model_score, test_surv_df, seq(365, max(test_surv_df$time), 365))
#cox_fit_score <- t(do.call("cbind", cox_fit_score)) #rows are the individuals, cols are the time
cox_fit_score <- list("coef" = summary(cox_model_score)$coef, "time" = cox_fit_score[1,], "risk" = cox_fit_score[cox_fit_score[,1] != cox_fit_score[1,1],], "conc" = score_conc_obj)


unchanged_res <- list("base" = cox_fit_base, "score" = cox_fit_score)



#################################### LEFT CENSOR ####################################
get_left_time <- do.call("rbind", lapply(test_surv_df$age_start_time + 365, function(x) seq(x, x+(13*365), 365)))

cox_model_base <- coxph(as.formula(paste0("Surv(age_start_time, age_end_time, pheno) ~ ", left_base_covars)), data = train_surv_df)

base_conc_obj <- concordance(cox_model_base, newdata=test_surv_df)

cox_fit_base <- faster_fit(cox_model_base, test_surv_df, get_left_time, "left")
#cox_fit_base <- t(do.call("cbind", cox_fit_base)) #rows are the individuals, cols are the time
cox_fit_base <- list("coef" = summary(cox_model_base)$coef, "time" = cox_fit_base[1,], "risk" = cox_fit_base[cox_fit_base[,1] != cox_fit_base[1,1],], "conc" = base_conc_obj)



cox_model_score <- coxph(as.formula(paste0("Surv(age_start_time, age_end_time, pheno) ~ ", left_score_covars)), data = train_surv_df)

score_conc_obj <- concordance(cox_model_score, newdata=test_surv_df)

cox_fit_score <- faster_fit(cox_model_score, test_surv_df, get_left_time, "left")
#cox_fit_score <- t(do.call("cbind", cox_fit_score)) #rows are the individuals, cols are the time
cox_fit_score <- list("coef" = summary(cox_model_score)$coef, "time" = cox_fit_score[1,], "risk" = cox_fit_score[cox_fit_score[,1] != cox_fit_score[1,1],], "conc" = score_conc_obj)


left_censor_res <- list("base" = cox_fit_base, "score" = cox_fit_score)





saveRDS(train_surv_df, "train_surv_df.RDS")
saveRDS(test_surv_df, "test_surv_df.RDS")


################################## CENSOR AT LAST EHR #############################################

time_res <- readRDS("time_res.RDS")
time_eid <- read.table("time_res_eid.txt", stringsAsFactors=F)
time_res$last_diag <- as.Date(time_res$last_diag, "1970-01-01")
time_res$first_diag <- as.Date(time_res$first_diag, "1970-01-01")
time_res$delay_last_diag <- time_res$last_diag + 365
time_res$delay_start_diag <- time_res$first_diag - 365

train_surv_df <- train_surv_df[train_surv_df$eid %in% time_eid[,1],]
train_time_res <- time_res[time_eid[,1] %in% train_surv_df$eid,]
train_time_eid <- time_eid[time_eid[,1] %in% train_surv_df$eid,,drop=F]
train_time_res <- train_time_res[order(train_time_eid[,1])[rank(train_surv_df$eid)],]

test_surv_df <- test_surv_df[test_surv_df$eid %in% time_eid[,1],]
test_time_res <- time_res[test_eid[,1] %in% test_surv_df$eid,]
test_time_eid <- test_eid[test_eid[,1] %in% test_surv_df$eid,,drop=F]
test_time_res <- test_time_res[order(test_time_eid[,1])[rank(test_surv_df$eid)],]

#train_surv_df$end_date <- as.Date(as.character(train_surv_df$end_date))
#test_surv_df$end_date <- as.Date(as.character(test_surv_df$end_date))

train_time_res$delay_last_diag[train_time_res$delay_last_diag >= train_surv_df$end_date] <- train_surv_df$end_date[train_time_res$delay_last_diag >= train_surv_df$end_date]
test_time_res$delay_last_diag[test_time_res$delay_last_diag >= test_surv_df$end_date] <- test_surv_df$end_date[test_time_res$delay_last_diag >= test_surv_df$end_date]

train_time_res$delay_start_diag[train_time_res$delay_start_diag <= train_surv_df$date_attend] <- train_surv_df$date_attende[train_time_res$delay_first_diag <= train_surv_df$date_attend]
test_time_res$delay_start_diag[test_time_res$delay_start_diag <= test_surv_df$date_attend] <- test_surv_df$date_attend[test_time_res$delay_first_diag <= test_surv_df$date_attend]

train_surv_df$ehr_time <- as.numeric(train_time_res$delay_last_diag - as.Date(train_surv_df$date_attend))
test_surv_df$ehr_time <- as.numeric(test_time_res$delay_last_diag - as.Date(test_surv_df$date_attend))

train_time_res <- train_time_res[!is.na(train_surv_df$ehr_time),]
test_time_res <- test_time_res[!is.na(test_surv_df$ehr_time),]

train_surv_df <- train_surv_df[!is.na(train_surv_df$ehr_time),]
test_surv_df <- test_surv_df[!is.na(test_surv_df$ehr_time),]

train_time_res <- train_time_res[is.finite(train_surv_df$ehr_time),]
test_time_res <- test_time_res[is.finite(test_surv_df$ehr_time),]

train_surv_df <- train_surv_df[is.finite(train_surv_df$ehr_time),]
test_surv_df <- test_surv_df[is.finite(test_surv_df$ehr_time),]





cox_model_base <- coxph(as.formula(paste0("Surv(ehr_time, pheno) ~ ", base_covars)), data = train_surv_df)

base_conc_obj <- concordance(cox_model_base, newdata=test_surv_df)

cox_fit_base <- faster_fit(cox_model_base, test_surv_df, seq(365, max(test_surv_df$ehr_time), 365))
#cox_fit_base <- t(do.call("cbind", cox_fit_base)) #rows are the individuals, cols are the time
cox_fit_base <- list("coef" = summary(cox_model_base)$coef, "time" = cox_fit_base[1,], "risk" = cox_fit_base[cox_fit_base[,1] != cox_fit_base[1,1],], "conc" = base_conc_obj)


cox_model_score <- coxph(as.formula(paste0("Surv(ehr_time, pheno) ~ ", score_covars)), data = train_surv_df)

score_conc_obj <- concordance(cox_model_score, newdata=test_surv_df)

cox_fit_score <- faster_fit(cox_model_score, test_surv_df, seq(365, max(test_surv_df$ehr_time), 365))
#cox_fit_score <- t(do.call("cbind", cox_fit_score)) #rows are the individuals, cols are the time
cox_fit_score <- list("coef" = summary(cox_model_score)$coef, "time" = cox_fit_score[1,], "risk" = cox_fit_score[cox_fit_score[,1] != cox_fit_score[1,1],], "conc" = score_conc_obj)


last_ehr_res <- list("base" = cox_fit_base, "score" = cox_fit_score)




##################################### CENSOR BY BOTH EHR ######################################



train_surv_df$ehr_start_time <- as.numeric(train_time_res$delay_start_diag - train_surv_df$dob)
train_surv_df$ehr_end_time <-  as.numeric(train_time_res$delay_last_diag - train_surv_df$dob)

test_surv_df$ehr_start_time <- as.numeric(test_time_res$delay_start_diag - test_surv_df$dob)
test_surv_df$ehr_end_time <-  as.numeric(test_time_res$delay_last_diag - test_surv_df$dob)

train_surv_df <- train_surv_df[!is.na(train_surv_df$ehr_start_time),]
test_surv_df <- test_surv_df[!is.na(test_surv_df$ehr_start_time),]
train_surv_df <- train_surv_df[is.finite(train_surv_df$ehr_start_time),]
test_surv_df <- test_surv_df[is.finite(test_surv_df$ehr_start_time),]

train_surv_df <- train_surv_df[train_surv_df$ehr_start_time < train_surv_df$ehr_end_time,]
test_surv_df <- test_surv_df[test_surv_df$ehr_start_time < test_surv_df$ehr_end_time,]

get_left_time <- do.call("rbind", lapply(test_surv_df$ehr_start_time + 365, function(x) seq(x, x+(13*365), 365)))



cox_model_base <- coxph(as.formula(paste0("Surv(ehr_start_time, ehr_end_time, pheno) ~ ", left_base_covars)), data = train_surv_df)

base_conc_obj <- concordance(cox_model_base, newdata=test_surv_df)

cox_fit_base <- faster_fit(cox_model_base, test_surv_df, get_left_time, "left")
#cox_fit_base <- t(do.call("cbind", cox_fit_base)) #rows are the individuals, cols are the time
cox_fit_base <- list("coef" = summary(cox_model_base)$coef, "time" = cox_fit_base[1,], "risk" = cox_fit_base[cox_fit_base[,1] != cox_fit_base[1,1],], "conc" = base_conc_obj)




cox_model_score <- coxph(as.formula(paste0("Surv(ehr_start_time, ehr_end_time, pheno) ~ ", left_score_covars)), data = train_surv_df)

score_conc_obj <- concordance(cox_model_score, newdata=test_surv_df)

cox_fit_score <- faster_fit(cox_model_score, test_surv_df, get_left_time, "left")
#cox_fit_score <- t(do.call("cbind", cox_fit_score)) #rows are the individuals, cols are the time
cox_fit_score <- list("coef" = summary(cox_model_score)$coef, "time" = cox_fit_score[1,], "risk" = cox_fit_score[cox_fit_score[,1] != cox_fit_score[1,1],], "conc" = score_conc_obj)


both_ehr_res <- list("base" = cox_fit_base, "score" = cox_fit_score)



#################################### NO CENSOR #####################################

train_surv_df <- readRDS("train_surv_df.RDS")
test_surv_df <- readRDS("test_surv_df.RDS")

train_surv_df$end_date <- as.Date("2020-05-31")
test_surv_df$end_date <- as.Date("2020-05-31")

train_surv_df$time <- as.numeric(train_surv_df$end_date - as.Date(train_surv_df$date_attend))
test_surv_df$time <- as.numeric(test_surv_df$end_date - as.Date(test_surv_df$date_attend))


cox_model_base <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", base_covars)), data = train_surv_df)

base_conc_obj <- concordance(cox_model_base, newdata=test_surv_df)

cox_fit_base <- faster_fit(cox_model_base, test_surv_df, seq(365, max(test_surv_df$time), 365))
#cox_fit_base <- t(do.call("cbind", cox_fit_base)) #rows are the individuals, cols are the time
cox_fit_base <- list("coef" = summary(cox_model_base)$coef, "time" = cox_fit_base[1,], "risk" = cox_fit_base[cox_fit_base[,1] != cox_fit_base[1,1],], "conc" = base_conc_obj)




cox_model_score <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", score_covars)), data = train_surv_df)

score_conc_obj <- concordance(cox_model_score, newdata=test_surv_df)

cox_fit_score <- faster_fit(cox_model_score, test_surv_df, seq(365, max(test_surv_df$time), 365))
#cox_fit_score <- t(do.call("cbind", cox_fit_score)) #rows are the individuals, cols are the time
cox_fit_score <- list("coef" = summary(cox_model_score)$coef, "time" = cox_fit_score[1,], "risk" = cox_fit_score[cox_fit_score[,1] != cox_fit_score[1,1],], "conc" = score_conc_obj)


no_censor_res <- list("base" = cox_fit_base, "score" = cox_fit_score)



##########################################################################################

all_res <- list("no_censor" = no_censor_res, "unchanged" = unchanged_res, "left_unchanged_censor" = left_censor_res, "ehr" = last_ehr_res, "left_ehr_censor" = both_ehr_res)

saveRDS(all_res, paste0("res/", author, ".res.RDS"))
