library(survival)
library(iCARE)
library(riskRegression)

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



if("sex" %in% colnames(surv_df)){
  base_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")
} else {
  base_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")
}


#train_surv_df <- train_surv_df[1:50000,]
#test_surv_df <- test_surv_df[1:3200,]
train_surv_df$id <- 1:nrow(train_surv_df)
test_surv_df$id <- 1:nrow(test_surv_df)


print("RR")
#####################################################################################
#####################################################################################
#riskregression #########################################
#cannote handle left censoring, so I can either have risk one the scale of time to assessment, or analyze risk in small age-match group

rr_model_base <- CSC(formula = as.formula(paste0("Hist(time, event_type) ~ ", base_covars)), data = train_surv_df)
rr_fit_base <- predict(rr_model_base, newdata = test_surv_df, cause = "diagnosis", times = seq(365, max(test_surv_df$time), 365))
cox_model_conc <- coxph(Surv(time, pheno) ~ rr_fit_base$absRisk[,14], data = test_surv_df)
rr_fit_base <- list("time" = rr_fit_base$times, "risk" = rr_fit_base$absRisk, "conc" = cox_model_conc$conc)

rr_model_score <- CSC(formula = as.formula(paste0("Hist(time, event_type) ~ ", score_covars)), data = train_surv_df)
rr_fit_score <- predict(rr_model_score, newdata = test_surv_df, cause = "diagnosis", times = seq(365, max(test_surv_df$time), 365))
cox_model_conc <- coxph(Surv(time, pheno) ~ rr_fit_score$absRisk[,14], data = test_surv_df)
rr_fit_coef <- summary(rr_model_score$models[[3]])$coefficients[nrow(summary(rr_model_score$models[[3]])$coefficients),]
rr_fit_score <- list("time" = rr_fit_score$times, "risk" = rr_fit_score$absRisk, "coef" = rr_fit_coef, "conc" = cox_model_conc$conc)



#rr_fit_subset_base <- list()
#rr_fit_subset_score <- list()
#rr_subset_eids <- list()
#rr_subset_coefs <- list()
#u_age <- sort(unique(round(train_surv_df$age)))
#u_age <- u_age[3:(length(u_age)-2)]
#ii <- 1
#for(i in 1:length(u_age)){
#  sub_train_df <- train_surv_df[round(train_surv_df$age) - 1 < u_age[i] & round(train_surv_df$age) + 1 > u_age[i],]
#  sub_test_df <- test_surv_df[round(test_surv_df$age) - 1 < u_age[i] & round(test_surv_df$age) + 1 > u_age[i],]
#  if(sum(sub_train_df$event_type == "diagnosis") > 10 & sum(sub_test_df$event_type == "diagnosis") > 10){
#
#    rr_model_base <- CSC(formula =  as.formula(paste0("Hist(time, event_type) ~ ", base_covars)), data = sub_train_df)
#    rr_fit_subset_base[[ii]] <- predict(rr_model_base, newdata = sub_test_df, cause = "diagnosis", times = seq(365, max(test_surv_df$time), 365))
#
#    rr_model_score <- CSC(formula = as.formula(paste0("Hist(time, event_type) ~ ", score_covars)), data = sub_train_df)
#    rr_fit_subset_score[[ii]] <- predict(rr_model_score, newdata = sub_test_df[i,], cause = "diagnosis", times = seq(365, max(test_surv_df$time), 365))
#    rr_subset_coefs[[ii]] <- summary(rr_model_score$models[[3]])$coefficients[nrow(summary(rr_model_score$models[[3]])$coefficients),]
#
#    rr_subset_eids[[ii]] <- sub_test_df$eid
#
#    ii <- ii + 1
#  } 
#}

#rr_subset_coefs <- Reduce("+", rr_subset_coefs)/length(rr_subset_coefs)

#rr_fit_subset_base <- list("time" = seq(365, max(test_surv_df$time), 365), "risk" = do.call("rbind", lapply(rr_fit_subset_base, function(x) x$absRisk)), "eids" = unlist(rr_subset_eids))
#rr_fit_subset_score <- list("time" = seq(365, max(test_surv_df$time), 365), "risk" = do.call("rbind", lapply(rr_fit_subset_score, function(x) x$absRisk)), "eids" = unlist(rr_subset_eids), "coef" = rr_subset_coefs)



all_rr_res <- list("base" = rr_fit_base, "score" = rr_fit_score)





print("MS")
#####################################################################################
#####################################################################################
#multi-state model ######################################3
train_surv_df$id <- 1:nrow(train_surv_df)
test_surv_df$id <- 1:nrow(test_surv_df)

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
  interval <- 1000
  for(ind in seq(1, nrow(new_data), interval)){
    print(ind)
    fit_prod <- survfit(fit_model, newdata = new_data[ind:(ind+interval-1),], se.fit = FALSE)
    if(input_model == "ms"){
      fit_prod <- cbind(fit_prod$time, fit_prod$cumhaz[,,2])
    } else {
      fit_prod <- cbind(fit_prod$time, fit_prod$cumhaz)
    }
    all_res[[counter]] <- fit_prod[get_inds(fit_prod[,1], desire_time),]
    counter <- counter + 1
  }
  return(all_res) #rows are the time, cols are the people
}

faster_fit <- function(fit_model, new_data, desire_time, model_type = "ms"){
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
  } else {
    fit_prod <- cbind(fit_prod$time, fit_prod$cumhaz)
  }
  fit_prod <- fit_prod[get_inds(fit_prod[,1], desire_time),]

  full_preds <- matrix(0, nrow = length(simple_preds), ncol = nrow(fit_prod))
  for(i in 1:ncol(full_preds)){
    mod <- lm(fit_prod[i,-1] ~ exp(simple_preds[1:1000]))
    full_preds[,i] <- mod$coef[1] + exp(simple_preds)*mod$coef[2]
  }
  return(rbind(desire_time, full_preds))
}

test_surv_df <- as.data.frame(test_surv_df)

#the third dimension of cumhaz follows the same order of states as the states element
#cannot get cumhaz at a specific time
ms_model_base <- coxph(as.formula(paste0("Surv(time, as.factor(event_type)) ~ ", base_covars)), id = id, data = train_surv_df)
ms_fit_base <- faster_fit(ms_model_base, test_surv_df, seq(365, max(test_surv_df$time), 365))
cox_model_conc <- coxph(Surv(time, pheno) ~ ms_fit_base[-1,14], data = test_surv_df)
ms_fit_base <- list("time" = ms_fit_base[1,], "risk" = ms_fit_base[-1,], "conc" = cox_model_conc$conc)


ms_model_score <- coxph(as.formula(paste0("Surv(time, as.factor(event_type)) ~ ", score_covars)), id = id, data = train_surv_df)
ms_model_coef <- summary(ms_model_score)$coef[nrow(summary(ms_model_score)$coef),]
ms_fit_score <- faster_fit(ms_model_score, test_surv_df, seq(365, max(test_surv_df$time), 365))
cox_model_conc <- coxph(Surv(time, pheno) ~ ms_fit_score[-1,14], data = test_surv_df)
ms_fit_score <- list("time" = ms_fit_score[1,], "risk" = ms_fit_score[-1,], "coef" = ms_model_coef, "conc" = cox_model_conc$conc)



#left censoring
#ms_model_base <- coxph(as.formula(paste0("Surv(time = age_start_time, time2 = age_end_time, event =  as.factor(event_type)) ~ ", base_covars)), id = id, data = train_surv_df)
#ms_fit_subset_base <- faster_fit(ms_model_base, test_surv_df, (60:80)*365)
#ms_fit_subset_base <- list("time" = ms_fit_subset_base[1,], "risk" = ms_fit_subset_base[-1,])


#ms_model_score <- coxph(as.formula(paste0("Surv(time = age_start_time, time2 = age_end_time, event =  as.factor(event_type)) ~ ", score_covars)), id = id, data = train_surv_df)
#ms_model_coef <- summary(ms_model_score)$coef[nrow(summary(ms_model_score)$coef),]
#ms_fit_subset_score <- faster_fit(ms_model_score, test_surv_df, (60:80)*365)
#ms_fit_subset_score <- list("time" = ms_fit_subset_score[1,], "risk" = ms_fit_subset_score[-1,], "coef" = ms_model_coef)


all_ms_res <- list("base" = ms_fit_base, "score" = ms_fit_score)





print("FG")
##############################################################################################3
###############################################################################################
#Fine and Gray ##########################################
fg_diag <- finegray(Surv(time, event_type) ~ ., data = train_surv_df, etype="diagnosis")

fg_model_base <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", base_covars)), data = fg_diag, weight=fgwt)
fg_fit_base <- faster_fit(fg_model_base, test_surv_df, seq(365, max(test_surv_df$time), 365), "fg")
cox_model_conc <- coxph(Surv(time, pheno) ~ fg_fit_base[-1,14], data = test_surv_df)
fg_fit_base <- list("time" = fg_fit_base[1,], "risk" = fg_fit_base[-1,], "conc" = cox_model_conc$conc)


fg_model_score <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", score_covars)), data = fg_diag, weight=fgwt)
fg_fit_coef <- summary(fg_model_score)$coef[nrow(summary(fg_model_score)$coef),]
fg_fit_score <- faster_fit(fg_model_score, test_surv_df, seq(365, max(test_surv_df$time), 365), "fg")
cox_model_conc <- coxph(Surv(time, pheno) ~ fg_fit_score[-1,14], data = test_surv_df)
fg_fit_score <- list("time" = fg_fit_score[1,], "risk" = fg_fit_score[-1,], "coef" = fg_fit_coef, "conc" = cox_model_conc$conc)



#left censoring
#fg_diag <- finegray(Surv(time = age_start_time, time2 = age_end_time, event_type) ~ ., id = id, data = train_surv_df, etype="diagnosis")

#fg_sub_model_base <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", base_covars)), data = fg_diag, weight=fgwt)
#fg_sub_fit_base <- faster_fit(fg_sub_model_base, test_surv_df, (60:80)*365, "fg")
#fg_sub_fit_base <- list("time" = fg_sub_fit_base[1,], "risk" = fg_sub_fit_base[-1,])


#fg_sub_model_score <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", score_covars)), data = fg_diag, weight=fgwt)
#fg_sub_fit_coef <- summary(fg_sub_model_score)$coef[nrow(summary(fg_sub_model_score)$coef),]
#fg_sub_fit_score <- faster_fit(fg_sub_model_score, test_surv_df, (60:80)*365, "fg")
#fg_sub_fit_score <- list("time" = fg_sub_fit_score[1,], "risk" = fg_sub_fit_score[-1,], "coef" = fg_sub_fit_coef)


all_fg_res <- list("base" = fg_fit_base, "score" = fg_fit_score)




print("IC")
######################################################################################
###################################################################################
#iCare ################################################
#note that this approach can also take factors, but requires much more work, so might as well stick with 
#converting everything to 1-hot encoding up front
#data("bc_data", package="iCARE")

if(file.exists(paste0("gbd_data/", author, ".inci_df.RDS"))){
death_rate <- readRDS("gbd_data/death.inci_df.RDS")
disease_prev <- readRDS(paste0("gbd_data/", author, ".inci_df.RDS"))

train_surv_df$diagnosis_yes <- (train_surv_df$event_type == "diagnosis")*1
test_surv_df$diagnosis_yes <- (test_surv_df$event_type == "diagnosis")*1


if("sex" %in% colnames(surv_df)){
base_model_formula <- as.formula(paste0("diagnosis_yes ~ ", "sex + age +", paste(paste0("PC", 1:10), collapse="+")))
base_cov_info <- list(list("name" = "sex", "type" = "continuous"),
                      list("name" = "age", "type" = "continuous"),
                      list("name" = "PC1", "type" = "continuous"),
                      list("name" = "PC2", "type" = "continuous"),
                      list("name" = "PC3", "type" = "continuous"),
                      list("name" = "PC4", "type" = "continuous"),
                      list("name" = "PC5", "type" = "continuous"),
                      list("name" = "PC6", "type" = "continuous"),
                      list("name" = "PC7", "type" = "continuous"),
                      list("name" = "PC8", "type" = "continuous"),
                      list("name" = "PC9", "type" = "continuous"),
                      list("name" = "PC10", "type" = "continuous"))

score_model_formula <- as.formula(paste0("diagnosis_yes ~ ", "sex + age +", paste(paste0("PC", 1:10), collapse="+"), "+score"))
score_cov_info <- list(list("name" = "sex", "type" = "continuous"),
		      list("name" = "age", "type" = "continuous"),
                      list("name" = "PC1", "type" = "continuous"),
                      list("name" = "PC2", "type" = "continuous"),
                      list("name" = "PC3", "type" = "continuous"),
                      list("name" = "PC4", "type" = "continuous"),
                      list("name" = "PC5", "type" = "continuous"),
                      list("name" = "PC6", "type" = "continuous"),
                      list("name" = "PC7", "type" = "continuous"),
                      list("name" = "PC8", "type" = "continuous"),
                      list("name" = "PC9", "type" = "continuous"),
                      list("name" = "PC10", "type" = "continuous"),
                      list("name" = "score", "type" = "continuous"))



ms_model <- coxph(as.formula(paste0("Surv(time, as.factor(event_type)) ~ ", "sex + age +", paste(paste0("PC", 1:10), collapse="+"))), id = id, data = train_surv_df)
base_model_log_or_1 <- ms_model$coef[grep("1:3", names(ms_model$coef))]
names(base_model_log_or_1) <- unlist(lapply(base_cov_info, function(x) x[[1]]))

ms_model <- coxph(as.formula(paste0("Surv(time, as.factor(event_type)) ~ ", "sex + age +", paste(paste0("PC", 1:10), collapse="+"), "+score")), id = id, data = train_surv_df)
score_model_log_or_1 <- ms_model$coef[grep("1:3", names(ms_model$coef))]
names(score_model_log_or_1) <- unlist(lapply(score_cov_info, function(x) x[[1]]))

} else {

base_model_formula <- as.formula(paste0("diagnosis_yes ~ ", "age +", paste(paste0("PC", 1:10), collapse="+")))
base_cov_info <- list(list("name" = "age", "type" = "continuous"),
                      list("name" = "PC1", "type" = "continuous"),
                      list("name" = "PC2", "type" = "continuous"),
                      list("name" = "PC3", "type" = "continuous"),
                      list("name" = "PC4", "type" = "continuous"),
                      list("name" = "PC5", "type" = "continuous"),
                      list("name" = "PC6", "type" = "continuous"),
                      list("name" = "PC7", "type" = "continuous"),
                      list("name" = "PC8", "type" = "continuous"),
                      list("name" = "PC9", "type" = "continuous"),
                      list("name" = "PC10", "type" = "continuous"))

score_model_formula <- as.formula(paste0("diagnosis_yes ~ ", "age +", paste(paste0("PC", 1:10), collapse="+"), "+score"))
score_cov_info <- list(list("name" = "age", "type" = "continuous"),
                      list("name" = "PC1", "type" = "continuous"),
                      list("name" = "PC2", "type" = "continuous"),
                      list("name" = "PC3", "type" = "continuous"),
                      list("name" = "PC4", "type" = "continuous"),
                      list("name" = "PC5", "type" = "continuous"),
                      list("name" = "PC6", "type" = "continuous"),
                      list("name" = "PC7", "type" = "continuous"),
                      list("name" = "PC8", "type" = "continuous"),
                      list("name" = "PC9", "type" = "continuous"),
                      list("name" = "PC10", "type" = "continuous"),
                      list("name" = "score", "type" = "continuous"))



ms_model <- coxph(as.formula(paste0("Surv(time, as.factor(event_type)) ~ ", "age +", paste(paste0("PC", 1:10), collapse="+"))), id = id, data = train_surv_df)
base_model_log_or_1 <- ms_model$coef[grep("1:3", names(ms_model$coef))]
names(base_model_log_or_1) <- unlist(lapply(base_cov_info, function(x) x[[1]]))

ms_model <- coxph(as.formula(paste0("Surv(time, as.factor(event_type)) ~ ", "age +", paste(paste0("PC", 1:10), collapse="+"), "+score")), id = id, data = train_surv_df)
score_model_log_or_1 <- ms_model$coef[grep("1:3", names(ms_model$coef))]
names(score_model_log_or_1) <- unlist(lapply(score_cov_info, function(x) x[[1]]))

}



curr_death_rate <- death_rate[,c(1,2)]
colnames(curr_death_rate) <- c("age", "val")
curr_disease_prev <- disease_prev[,c(1,2)]
colnames(curr_disease_prev) <- c("age", "val")


if("sex" %in% colnames(surv_df)){
  use_train_surv_df_base <- train_surv_df[,colnames(train_surv_df) %in% c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
  use_test_surv_df_base <- test_surv_df[,colnames(test_surv_df) %in% c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]

  use_train_surv_df_score <- train_surv_df[,colnames(train_surv_df) %in% c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "score")]
  use_test_surv_df_score <- test_surv_df[,colnames(test_surv_df) %in% c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "score")]
} else {

  use_train_surv_df_base <- train_surv_df[,colnames(train_surv_df) %in% c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
  use_test_surv_df_base <- test_surv_df[,colnames(test_surv_df) %in% c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]

  use_train_surv_df_score <- train_surv_df[,colnames(train_surv_df) %in% c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "score")]
  use_test_surv_df_score <- test_surv_df[,colnames(test_surv_df) %in% c("age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "score")]
}

#check_model_inputs(apply.cov.profile = use_test_surv_df , model.log.RR = base_mod_log_or_1, model.ref.dataset = use_train_surv_df, model.cov.info = base_cov_info, model.formula = base_model_formula)
get_left_time <- do.call("rbind", lapply(test_surv_df$age_start_time + 365, function(x) seq(x, x+(13*365), 365)))
current_time <- (51:90)*365

fix_func <- function(p, current_risk){
  funcy <- splinefun(current_time, current_risk[p,])
  new_vec <- funcy(get_left_time[p,])
  new_vec[get_left_time[p,] < current_time[1]] <- NA
  return(new_vec)
}


final_icare_base_res <- list()
final_icare_score_res <- list()
final_base_conc <- list()
final_score_conc <- list()

for(kk in 1:3){
curr_death_rate <- death_rate[,c(1,kk+1)]
colnames(curr_death_rate) <- c("age", "val")
curr_disease_prev <- disease_prev[,c(1,kk+1)]
colnames(curr_disease_prev) <- c("age", "val")


icare_full_res <- list()
icare_score_res <- list()
for(i in 1:40){
  icare_full_res[[i]] <-  computeAbsoluteRisk(model.formula = base_model_formula,
                                         model.cov.info = base_cov_info,
                                         model.log.RR = base_model_log_or_1,
                                         model.ref.dataset = use_train_surv_df_base,
                                         model.disease.incidence.rates = curr_disease_prev,
                                         model.competing.incidence.rates = curr_death_rate,
                                         apply.age.start = 51,
                                         apply.age.interval.length = i,
                                         apply.cov.profile = use_test_surv_df_base,
                                         return.refs.risk = TRUE)
  icare_full_res[[i]] <- icare_full_res[[i]]$details
  icare_full_res[[i]] <- icare_full_res[[i]]$Risk_Estimate



  icare_score_res[[i]] <-  computeAbsoluteRisk(model.formula = score_model_formula,
                                         model.cov.info = score_cov_info,
                                         model.log.RR = score_model_log_or_1,
                                         model.ref.dataset = use_train_surv_df_score,
                                         model.disease.incidence.rates = curr_disease_prev,
                                         model.competing.incidence.rates = curr_death_rate,
                                         apply.age.start = 51,
                                         apply.age.interval.length = i,
                                         apply.cov.profile = use_test_surv_df_score,
                                         return.refs.risk = TRUE)
  icare_score_res[[i]] <- icare_score_res[[i]]$details
  icare_score_res[[i]] <- icare_score_res[[i]]$Risk_Estimate

}

final_icare_score_res[[kk]] <- do.call("cbind", icare_score_res)
final_icare_base_res[[kk]] <- do.call("cbind", icare_full_res)

final_icare_score_res[[kk]] <- do.call("rbind", lapply(1:nrow(final_icare_score_res[[kk]]), fix_func, final_icare_score_res[[kk]]))
final_icare_base_res[[kk]] <- do.call("rbind", lapply(1:nrow(final_icare_base_res[[kk]]), fix_func, final_icare_base_res[[kk]]))

cox_model_conc <- coxph(Surv(time, pheno) ~ final_icare_base_res[[kk]][,14], data = test_surv_df)
final_base_conc[[kk]] <- cox_model_conc$conc
cox_model_conc <- coxph(Surv(time, pheno) ~ final_icare_score_res[[kk]][,14], data = test_surv_df)
final_score_conc[[kk]] <- cox_model_conc$conc
}

names(final_icare_score_res) <- c("gbd", "holt", "combo")
names(final_icare_base_res) <- c("gbd", "holt", "combo")
names(final_base_conc) <- c("gbd", "holt", "combo")
names(final_score_conc) <- c("gbd", "holt", "combo")


all_icare_res <- list("base" = final_icare_base_res, "score" = final_icare_score_res, "conc_base" = final_base_conc, "conc_score" = final_score_conc)




} else {

all_icare_res <- NULL
}


###########################################

all_res <- list("rr" = all_rr_res, "ms" = all_ms_res, "fg" = all_fg_res, "icare" = all_icare_res)
saveRDS(all_res, paste0("res/", author, ".RDS"))

