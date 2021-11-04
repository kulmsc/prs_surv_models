library(survival)

author = commandArgs(trailingOnly=TRUE)
#author = "wray"


system("zcat test_df.gz | cut -f1 -d',' > test_eid")
system("zcat train_df.gz | cut -f1 -d',' > train_eid")

#system("rm train_df.gz")
#system("rm test_df.gz")


#read in the UKBB data
surv_df <- readRDS(paste0("../init_data/surv_data/survdf.", author, ".RDS"))
surv_df$time <- as.numeric(as.Date(surv_df$end_date) - as.Date(surv_df$date_attend)) #time is counted from the date of assessment


#split the dataset
train_eid <- read.table("train_eid", stringsAsFactors = F, header=T)
test_eid <- read.table("test_eid", stringsAsFactors = F, header=T)

train_surv_df <- surv_df[surv_df$eid %in% train_eid[,1],]
test_surv_df <- surv_df[surv_df$eid %in% test_eid[,1],]



if("sex" %in% colnames(surv_df)){
  base_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")
} else {
  base_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  score_covars <- c("age +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score")
}

####################################################################

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



##################################################################3###

train_cox_preds <- read.table("train_pred.txt.gz", stringsAsFactors=F, sep=",")
cox_preds <- read.table("static_pred.txt.gz", stringsAsFactors=F, sep = ",")
test_surv_df$python_preds <- cox_preds[,2]

cox_model_python <- coxph(Surv(time, pheno) ~ python_preds, data = test_surv_df)
python_res <- list("coef" = summary(cox_model_python)$coef, "conc" = cox_model_python$conc, "pred" = cox_preds[,2],
                   "df" = test_surv_df[,colnames(test_surv_df) %in% c("time", "pheno", "score")])





########################################################################33

coef <- read.table("coef.txt.gz", stringsAsFactors=F, sep = ",")

good_inds <- which(abs(coef[,2]) >= sort(abs(coef[,2]), decreasing=T)[10])
system(paste0("zcat test_df.gz | cut -d',' -f", paste(good_inds, collapse=","), " > use_test_data"))
system(paste0("zcat train_df.gz | cut -d',' -f", paste(good_inds, collapse=","), " > use_train_data"))
use_test_data <- read.table("use_test_data", stringsAsFactors=F, sep = ",", header=T)
use_train_data <- read.table("use_train_data", stringsAsFactors=F, sep = ",", header=T)
#pred, top 5, 10, each of top 3 alone

write.table(data.frame("coef" = colnames(use_test_data), "val" = coef[good_inds,2]), paste0("res/", author, ".coef.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
colnames(use_test_data) <- paste("col", good_inds, sep="_")
colnames(use_train_data) <- paste("col", good_inds, sep="_")

train_surv_df <- cbind(train_surv_df, use_train_data)
test_surv_df <- cbind(test_surv_df, use_test_data)

train_surv_df$big_pred <- train_cox_preds[,2]
test_surv_df$big_pred <- cox_preds[,2]


#####################################################################


all_base_covars <- list()
all_score_covars <- list()

top_takes <- c(3, 5, 10)

for(i in 1:3){
  curr_inds <- which(abs(coef[,2]) >= sort(abs(coef[,2]), decreasing=T)[top_takes[i]])
  curr_names <- paste("col", curr_inds, sep="_")

  all_base_covars[[i]] <- paste(base_covars, "+", paste(curr_names, collapse="+"))
  all_score_covars[[i]] <- paste(score_covars, "+", paste(curr_names, collapse="+"))
}

top_takes <- 1:3

for(i in 1:3){
  curr_inds <- which(abs(coef[,2]) == sort(abs(coef[,2]), decreasing=T)[top_takes[i]])
  curr_names <- paste("col", curr_inds, sep="_")

  all_base_covars[[i+3]] <- paste(base_covars, "+", paste(curr_names, collapse="+"))
  all_score_covars[[i+3]] <- paste(score_covars, "+", paste(curr_names, collapse="+"))
}



all_base_covars[[7]] <- paste0(base_covars, " + big_pred")
all_score_covars[[7]] <- paste0(score_covars, " + big_pred")


################################################################

all_res <- list()



for(i in 1:7){
  cox_model_base <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", all_base_covars[[i]])), data = train_surv_df)

  base_conc_obj <- concordance(cox_model_base, newdata=test_surv_df)

  #cox_fit_base <- get_survfit(cox_model_base, test_surv_df, seq(365, max(test_surv_df$time), 365))
  cox_fit_base <- faster_fit(cox_model_base, test_surv_df, seq(365, max(test_surv_df$time), 365))
  #cox_fit_base <- t(do.call("cbind", cox_fit_base)) #rows are the individuals, cols are the time
  cox_fit_base <- list("coef" = summary(cox_model_base)$coef, "time" = cox_fit_base[1,], "risk" = cox_fit_base[-1,], "conc" = base_conc_obj)


  cox_model_score <- coxph(as.formula(paste0("Surv(time, pheno) ~ ", all_score_covars[[i]])), data = train_surv_df)

  score_conc_obj <- concordance(cox_model_score, newdata=test_surv_df)

  #cox_fit_score <- get_survfit(cox_model_score, test_surv_df, seq(365, max(test_surv_df$time), 365))
  cox_fit_score <- faster_fit(cox_model_score, test_surv_df, seq(365, max(test_surv_df$time), 365))
  #cox_fit_score <- t(do.call("cbind", cox_fit_score)) #rows are the individuals, cols are the time
  cox_fit_score <- list("coef" = summary(cox_model_score)$coef, "time" = cox_fit_score[1,], "risk" = cox_fit_score[-1,], "conc" = score_conc_obj)

  all_res[[i]] <- list("base" = cox_fit_base, "score" = cox_fit_score)

}



all_res[[8]] <- python_res
names(all_res) <- c("top3_coef", "top5_coef", "top10_coef", "just1_coef", "just2_coef", "just_3_coef", "big_pred", "python")

saveRDS(all_res, paste0("res/", author, ".RDS"))
