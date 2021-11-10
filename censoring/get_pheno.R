library(survival)

#author = "bentham"
author = commandArgs(trailingOnly=TRUE)
res <- readRDS(paste0("res/", author, ".res.RDS"))

###########################################################################

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

first_phen <- test_surv_df$pheno

#############################################################################


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

second_phen <- test_surv_df$pheno


##################################################################################


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

third_phen <- test_surv_df$pheno

###################################################################################

res[["phen"]] <- list(first_phen, second_phen, third_phen)
saveRDS(res, paste0("res/", author, ".res.RDS"))
