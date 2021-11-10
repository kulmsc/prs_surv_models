
author = commandArgs(trailingOnly=TRUE)
#author = "shah"

system("rm train_df.gz")
system("rm test_df.gz")

#read in the UKBB data
surv_df <- readRDS(paste0("../init_data/surv_data/survdf.", author, ".RDS"))
surv_df$time <- as.numeric(as.Date(surv_df$end_date) - as.Date(surv_df$date_attend)) #time is counted from the date of assessment



#split the dataset
train_frac <- 0.6
test_frac <- 1 - train_frac
train_eid <- read.table(paste0("~/athena/doc_score/qc/cv_files/train_eid.", train_frac, ".txt"), stringsAsFactors=F)
test_eid <- read.table(paste0("~/athena/doc_score/qc/cv_files/test_eid.", test_frac, ".txt"), stringsAsFactors=F)
train_surv_df <- surv_df[surv_df$eid %in% train_eid[,1],]
test_surv_df <- surv_df[surv_df$eid %in% test_eid[,1],]


#get big data
big_data <- readRDS(paste0("../init_data/clean_big_data/", author, ".big_data.RDS"))
bad_val <- apply(big_data[[1]], 2, function(x) sum(is.na(x)))
big_data[[1]] <- big_data[[1]][,bad_val == 0]
bad_val <- apply(big_data[[2]], 2, function(x) sum(is.na(x)))
big_data[[2]] <- big_data[[2]][,bad_val == 0]

big_data[[1]] <- big_data[[1]][,colnames(big_data[[1]]) %in% colnames(big_data[[2]])]
big_data[[2]] <- big_data[[2]][,colnames(big_data[[2]]) %in% colnames(big_data[[1]])]


train_eid <- read.table("../init_data/clean_big_data/train_full_eids_big_data_order.txt", stringsAsFactors=F)
if(!all(big_data[[1]]$eid[1:10] == train_eid[1:10,1])){
mod <- lm(train_eid[1:3,1] ~ big_data[[1]]$eid[1:3])
if(summary(mod)$r.squared != 1){
  print("STOP")
  mod <- lm(train_eid[1:2,1] ~ big_data[[1]]$eid[1:2])
  if(all(round(big_data[[1]]$eid*mod$coef[2] + mod$coef[1])[1:10] %in% train_eid[,1])){
    print("GOOD")
  } else {
    print("BAD!!!!!!!!!!!!!!!!!!!!!!!!!")
  }
}

big_data[[1]]$eid <- round(big_data[[1]]$eid*mod$coef[2] + mod$coef[1])
print("R SQUARED")
print(summary(mod)$r.squared)
}

test_eid <- read.table("../init_data/clean_big_data/test_full_eids_big_data_order.txt", stringsAsFactors=F)
if(!all(big_data[[2]]$eid[1:10] == test_eid[1:10,1])){
mod <- lm(test_eid[1:3,1] ~ big_data[[2]]$eid[1:3])
if(summary(mod)$r.squared != 1){
  print("STOP")
  mod <- lm(test_eid[1:2,1] ~ big_data[[2]]$eid[1:2])
  if(all(round(big_data[[2]]$eid*mod$coef[2] + mod$coef[1])[1:10] %in% test_eid[,1])){
    print("GOOD")
  } else {
    print("BAD!!!!!!!!!!!!!!!!!!!!!!!!!")
  }
}
big_data[[2]]$eid <- round(big_data[[2]]$eid*mod$coef[2] + mod$coef[1])
print("R SQUARED")
print(summary(mod)$r.squared)
}

train_surv_df <- train_surv_df[train_surv_df$eid %in% big_data[[1]]$eid,]
big_data[[1]] <- big_data[[1]][big_data[[1]]$eid %in% train_surv_df$eid,]
big_data[[1]] <- big_data[[1]][order(big_data[[1]]$eid)[rank(train_surv_df$eid)],]
gz1 <- gzfile("train_df.gz", "w")
write.csv(big_data[[1]][,2:ncol(big_data[[1]])], gz1, row.names=F,quote=F)
close(gz1)

test_surv_df <- test_surv_df[test_surv_df$eid %in% big_data[[2]]$eid,]
big_data[[2]] <- big_data[[2]][big_data[[2]]$eid %in% test_surv_df$eid,]
big_data[[2]] <- big_data[[2]][order(big_data[[2]]$eid)[rank(test_surv_df$eid)],]
gz1 <- gzfile("test_df.gz", "w")
write.csv(big_data[[2]][,2:ncol(big_data[[2]])], gz1, row.names=F,quote=F)
close(gz1)

train_surv_df <- data.frame("time" = train_surv_df$time, "status" = train_surv_df$pheno)
write.table(train_surv_df, "y_train.csv", row.names = F, col.names = T, quote = F, sep = ",")

test_surv_df <- data.frame("time" = test_surv_df$time, "status" = test_surv_df$pheno)
write.table(test_surv_df, "y_test.csv", row.names = F, col.names = T, quote = F, sep = ",")


