library(stringr)
library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(rmda)
library(forecast)
theme_set(theme_cowplot())


author="death"
#author = commandArgs(trailingOnly=TRUE)



if(author == "death"){
  surv_data <- readRDS(paste0("surv_data/survdf.jin.RDS"))
  surv_data <- surv_data[surv_data$pheno == 0,]
  surv_data$pheno[surv_data$is_death_date == 1] <- 1
} else {
  surv_data <- readRDS(paste0("surv_data/survdf.", author, ".RDS"))
}

surv_data$age <- (as.Date(surv_data$final_rec_date) - as.Date(surv_data$dob))/365
surv_data$time <- as.numeric(as.Date(surv_data$end_date) - as.Date(surv_data$date_attend)) #time is counted from the date of assessment


############### GBD #########################

translator <- read.table("author_translator.txt", stringsAsFactors=F, sep="\t")

if(author == "death"){
  gen_inci <- read.table("disease_prev/combo_mort.txt", stringsAsFactors = F, header = T)
  gen_inci <- gen_inci[gen_inci$age > 45,]
  gbd_inci <- list("2019" = gen_inci)
  
} else {
  gen_inci <- read.csv("gbd_inci.csv", stringsAsFactors = F, header=T)
  gen_inci <- gen_inci[gen_inci$cause == translator[translator[,2] == author,1] & gen_inci$sex == "Both",]  

  u_year <- sort(unique(gen_inci$year))
  gbd_inci <- list()
  
  for(i in 1:length(u_year)){
    #inci
    sub_inci <- gen_inci[gen_inci$year == u_year[i],]
    sub_inci$age <- seq(2, 92, by = 5)
    sub_inci <- sub_inci[sub_inci$age > 37,]
    
    new_inci <- data.frame("age" = min(sub_inci$age):max(sub_inci$age))
    
    funcy <- splinefun(sub_inci$age, sub_inci$val, method = "fmm")
    new_inci$val <- funcy(new_inci$age)
    
    funcy <- splinefun(sub_inci$age, sub_inci$upper, method = "fmm")
    new_inci$upper <- funcy(new_inci$age)
    
    funcy <- splinefun(sub_inci$age, sub_inci$lower, method = "fmm")
    new_inci$lower <- funcy(new_inci$age)
    
    new_inci[,2:4] <- new_inci[,2:4]/100000
    gbd_inci[[i]] <- new_inci
  }
  
  names(gbd_inci) <- u_year

}


################# UKBB #######################


all_years <- 2010:2019
ukb_prev <- list()


for(k in 1:length(all_years)){
  
  comp_ukb_prev <- data.frame(matrix(0, nrow = nrow(gbd_inci[[1]]), ncol = 5))
  colnames(comp_ukb_prev) <-  c(colnames(gbd_inci[[1]]), "n")
  comp_ukb_prev[,1] <- gbd_inci[[1]][,1]

  #current age is measured from July 1 - with floor function
  surv_data$current_age <- floor(as.numeric(as.Date(paste0(all_years[k], "-07-01")) - surv_data$dob)/365)
  for(j in 1:nrow(gbd_inci[[1]])){
    age_surv_data <- surv_data[surv_data$current_age == gbd_inci[[1]][j,1],]
    
    if(nrow(age_surv_data) > 0){

      comp_ukb_prev[j,2] <- sum(substr(age_surv_data$end_date[age_surv_data$pheno == 1], 1, 4) == all_years[k])/nrow(age_surv_data)
      comp_ukb_prev[j,5] <- nrow(age_surv_data)
      
      if(comp_ukb_prev[j,2] > 0){
        pro_se <- sqrt((comp_ukb_prev[j,2] * (1 - comp_ukb_prev[j,2]))/nrow(age_surv_data))
        comp_ukb_prev[j,3] <- (comp_ukb_prev[j,2] + pro_se)
        comp_ukb_prev[j,4] <- (comp_ukb_prev[j,2] - pro_se)
        comp_ukb_prev[j,2] <- comp_ukb_prev[j,2]
      }
    }
  }
  
  ukb_prev[[k]] <- comp_ukb_prev
}

names(ukb_prev) <- all_years




#################################################3##########################


curr_inci <- gbd_inci[["2019"]]
curr_inci$val <- cumsum(curr_inci$val)
curr_inci$upper <- cumsum(curr_inci$upper)
curr_inci$lower <- cumsum(curr_inci$lower)
curr_inci$n <- 500000

curr_uk <- ukb_prev[["2019"]]
curr_uk$val <- cumsum(curr_uk$val)
curr_uk$upper <- cumsum(curr_uk$upper)
curr_uk$lower <- cumsum(curr_uk$lower)

df <- rbind(curr_inci, curr_uk)

df$type <- "UK"
df$type[df$n == 500000] <- "GBD"
df$val[df$n < 800] <- 0

#Now making the time series objects
ts_df <- data.frame("gbd" = df$val[df$type == "GBD"], "uk" = df$val[df$type == "UK"], "age" = df$age[df$type == "UK"])
ts_df <- ts_df[ts_df$age > 50,]
ts_df$uk[ts_df$uk == 0] <- NA

simple_obj <- ts(ts_df$uk, start = ts_df$age[1], end = max(ts_df$age[!is.na(ts_df$uk)]) )
ts_df$nholt <- c(ts_df$uk[!is.na(ts_df$uk)], as.numeric(holt(simple_obj, h=sum(is.na(ts_df$uk)), damped = T)$mean))
ts_df$gholt <- c(ts_df$uk[1:10], rep(NA, nrow(ts_df)-10))
for(i in 11:sum(!is.na(ts_df$uk))){
  ts_df$gholt[i] <- as.numeric(holt(ts(ts_df$uk[1:i], start = ts_df$age[1], end =ts_df$age[1]+i ), 1, damped = T)$mean)
}

ts_obj <- ts(ts_df[,1:5], start = ts_df$age[1], end = max(ts_df$age[!is.na(ts_df$uk)]) )
ts_mod <- tslm(uk ~ gholt + gbd + gbd:age, data = ts_obj)

cast_data <- data.frame("gbd" = ts_df$gbd[is.na(ts_df$uk)],
                        "age" = ts_df$age[is.na(ts_df$uk)],
                        "gholt" = ts_df$nholt[is.na(ts_df$uk)])
cast_res <- forecast(ts_mod, newdata =cast_data, h = sum(is.na(tf_df$uk)))

ts_df$combo_pred <- ts_df$uk
ts_df$combo_pred[is.na(ts_df$uk)] <- as.numeric(cast_res$mean)




plot_df <- data.frame("val" = c(ts_df$gbd, ts_df$uk[!is.na(ts_df$uk)], ts_df$nholt[is.na(ts_df$uk)], ts_df$combo_pred[is.na(ts_df$uk)]),
                      "age" = c(ts_df$age, ts_df$age[!is.na(ts_df$uk)], ts_df$age[is.na(ts_df$uk)], ts_df$age[is.na(ts_df$uk)]),
                      "type" = c(rep("GBD", nrow(ts_df)), rep("UK", sum(!is.na(ts_df$uk))), rep("HOLT", sum(is.na(ts_df$uk))), rep("COMBO", sum(is.na(ts_df$uk)))))

the_plot <- ggplot(plot_df, aes(age, val, color = type)) + geom_point() +
  labs(x = "Age", y = "Cumulative Incidence", color = "Source")
plot(the_plot)
ggsave(paste0("correct_inci_dfs/", author, ".inci.png"), the_plot, width = 5.5, height=4)

#more complicated calculation needed
inci_df <- data.frame("age" = ts_df$age,
                      "gbd" = ts_df$gbd - c(0, ts_df$gbd[-nrow(ts_df)]),
                      "holt" = ts_df$nholt - c(0, ts_df$nholt[-nrow(ts_df)]),
                      "combo" = ts_df$nholt - c(0, ts_df$nholt[-nrow(ts_df)]))
saveRDS(inci_df, paste0("correct_inci_dfs/", author, ".inci_df.RDS"))
