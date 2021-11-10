
#author <- "bentham"
author <- commandArgs(trailingOnly=TRUE)
print(author)

system(paste0("Rscript get_y_data.R ", author))

res <- readRDS(paste0("res/", author, ".RDS"))

system("zcat test_df.gz | cut -f1 -d',' > test_eid")
system("zcat train_df.gz | cut -f1 -d',' > train_eid")

surv_df <- readRDS(paste0("../init_data/surv_data/survdf.", author, ".RDS"))
surv_df$time <- as.numeric(as.Date(surv_df$end_date) - as.Date(surv_df$date_attend)) #time is counted from the date of assessment

train_eid <- read.table("train_eid", stringsAsFactors = F, header=T)
test_eid <- read.table("test_eid", stringsAsFactors = F, header=T)

train_surv_df <- surv_df[surv_df$eid %in% train_eid[,1],]
test_surv_df <- surv_df[surv_df$eid %in% test_eid[,1],]

res[["eid"]] <- test_surv_df$eid
res[["pheno"]] <- test_surv_df$pheno

print(author)
saveRDS(res, paste0("res/", author, ".RDS"))

system("rm test_df.gz")
system("rm train_df.gz")
