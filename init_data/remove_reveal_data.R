library(Hmisc)

#author <- "shah"
author = commandArgs(trailingOnly=TRUE)

big_data <- read.table("clean_big_data/init_clean_big_data.txt.gz", stringsAsFactors=F, header=T, sep="\t")
#manually add back in some survey data
survey_data <- read.table("raw_big_data/big_data.just_survey.stop_at_baseline.txt.gz", stringsAsFactors=F, header=T)
na_count <- apply(survey_data, 2, function(x) sum(is.na(x)))
survey_data <- survey_data[,na_count/nrow(survey_data) < 0.5]
survey_data <- survey_data[,!(colnames(survey_data) %in% big_data)]
big_data <- cbind(big_data, survey_data)

#big_data <- readRDS("../get_data/big_data.ehr_stop_at_baseline.RDS")
#colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))] <- paste0("X", colnames(big_data)[which(!is.na(as.numeric(colnames(big_data))))])

pheno_decoder <- read.table("~/athena/doc_score/analyze_score/descript_defs/author_defs", stringsAsFactors=F, header=T)
useful_decoder <- as.character(pheno_decoder[tolower(pheno_decoder[,1]) == author,4:9][1,])
useful_decoder[useful_decoder == "NA"] <- NA
real_names <- c("CANCER", "NONCANCER", "ICD9", "ICD10", "OPCS", "MEDS")
bad_names <- c()
bad_people <- list()
j <- 1


#remove individuals who were already diagnosed
for(i in 1:length(useful_decoder)){
  if(!is.na(useful_decoder[i])){
    if(grepl("|", useful_decoder[i])){
      for(subcode in strsplit(useful_decoder[i], "|", fixed = T)[[1]]){
        bad_names <- c(bad_names, paste0(real_names[i], "_", subcode))
        if(real_names[i] %in% c("ICD9", "ICD10", "OPCS")){
          bad_people[[j]] <- big_data$eid[big_data[,which(colnames(big_data) == paste0(real_names[i], "_", subcode))] == 1]
          j <- j + 1
        }
      }

    } else {
      bad_names <- c(bad_names, paste0(real_names[i], "_", useful_decoder[i]))
      bad_people[[j]] <- big_data$eid[big_data[,which(colnames(big_data) == paste0(real_names[i], "_", useful_decoder[i]))] == 1]
      j <- j + 1
    }

  }
}


bad_people <- unique(unlist(bad_people))
big_data <- big_data[!(big_data$eid %in% bad_people),]

big_data <- big_data[,-which(colnames(big_data) %in% bad_names)]




###########################################################################################
###########################################################################################

train_frac <- 0.6
test_frac <- 1 - train_frac
train_eid <- read.table(paste0("~/athena/doc_score/qc/cv_files/train_eid.", train_frac, ".txt"), stringsAsFactors=F)
test_eid <- read.table(paste0("~/athena/doc_score/qc/cv_files/test_eid.", test_frac, ".txt"), stringsAsFactors=F)

inc_vars <- c("SURVEY_X31.0.0", #sex
              "SURVEY_X34.0.0", #year of birth
              "SURVEY_X50.0.0", #standing height
              "SURVEY_X135.0.0", #number illnesses
              "SURVEY_X137.0.0", #number of medication
              "SURVEY_X1647.0.0", #born in UK
              "SURVEY_X2178.0.0", #health ratio
              "SURVEY_X20458.0.0", #general happiness
              "SURVEY_X23104.0.0", #BMI
              "SURVEY_X845.0.0", #years education
              "SURVEY_X738.0.0", #income
              "SURVEY_X20401.0.0", #ever addicted
              "SURVEY_X30020.0.0") #hemoglobin

if(length(unique(big_data$SURVEY_X31.0.0)) == 1){
  big_data <- big_data[,-which(colnames(big_data) == "SURVEY_X31.0.0")]
  inc_vars <- inc_vars[inc_vars != "SURVEY_X31.0.0"]
}

inc_vars <- inc_vars[inc_vars %in% colnames(big_data)]

train_data <- big_data[big_data$eid %in% train_eid[,1],]
test_data <- big_data[big_data$eid %in% test_eid[,1],]

#Train Data
the_formula <- as.formula(paste0("~ ", paste(inc_vars, collapse = " + ")))
f <- aregImpute(the_formula, nk=0, tlinear=FALSE, data = train_data, B = 50, n.impute = 1)
for(i in 1:length(inc_vars)){
  if(!is.null(f$imputed[[i]]) & any(is.na(train_data[inc_vars[i]]))){
    train_data[inc_vars[i]][is.na(train_data[inc_vars[i]])] <- unname(f$imputed[[inc_vars[i]]][,1])
  }
}


na_count <- apply(train_data, 2, function(x) sum(is.na(x)))

bad_names <- names(na_count)[na_count > 0]
for(i in 1:length(bad_names)){
  the_formula <- as.formula(paste0("~ ", paste(c(inc_vars, bad_names[i]), collapse = " + ")))
  f <- aregImpute(the_formula, nk=0, tlinear=FALSE, data = train_data, B = 50, n.impute = 1)
  train_data[bad_names[i]][is.na(train_data[bad_names[i]])] <- unname(f$imputed[[bad_names[i]]][,1])
}




#Test Data
the_formula <- as.formula(paste0("~ ", paste(inc_vars, collapse = " + ")))
f <- aregImpute(the_formula, nk=0, tlinear=FALSE, data = test_data, B = 50, n.impute = 1)
for(i in 1:length(inc_vars)){
  if(!is.null(f$imputed[[i]]) & any(is.na(test_data[inc_vars[i]]))){
    test_data[inc_vars[i]][is.na(test_data[inc_vars[i]])] <- unname(f$imputed[[inc_vars[i]]][,1])
  }
}


na_count <- apply(test_data, 2, function(x) sum(is.na(x)))

bad_names <- names(na_count)[na_count > 0]
for(i in 1:length(bad_names)){
  the_formula <- as.formula(paste0("~ ", paste(c(inc_vars, bad_names[i]), collapse = " + ")))
  f <- aregImpute(the_formula, nk=0, tlinear=FALSE, data = test_data, B = 50, n.impute = 1)
  test_data[bad_names[i]][is.na(test_data[bad_names[i]])] <- unname(f$imputed[[bad_names[i]]][,1])
}








for(i in 2:ncol(train_data)){
  train_data[,i] <- (train_data[,i] - min(train_data[,i]))/(max(train_data[,i]) - min(train_data[,i]))
  test_data[,i] <- (test_data[,i] - min(test_data[,i]))/(max(test_data[,i]) - min(test_data[,i]))
}


saveRDS(list("train" = train_data, "test" = test_data),  paste0("clean_big_data/", author, ".big_data.RDS"))
