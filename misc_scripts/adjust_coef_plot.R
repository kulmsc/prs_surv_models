library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
library(forestplot)
library(data.table)
theme_set(theme_cowplot())


analysis_type <- "adjustment" 
source("common.R")

all_files <- list.files("fuel_for_plot/adjustment_coefs/")
all_coefs <- list()
ten_coefs <- list()
five_coefs <- list()
three_coefs <- list()

for(i in 1:length(all_files)){
  all_coefs[[i]] <- read.table(paste0("fuel_for_plot/adjustment_coefs/", all_files[i]), stringsAsFactors = F, header = T)
  all_coefs[[i]]$disease <- swap_names(strsplit(all_files[i], ".", fixed = T)[[1]][[1]], disease_names)

  
  ten_coefs[[i]] <- all_coefs[[i]][abs(all_coefs[[i]]$val) >= sort(abs(all_coefs[[i]]$val), decreasing = T)[ifelse(nrow(all_coefs[[i]]) >= 10, 10, nrow(all_coefs[[i]]))] & all_coefs[[i]]$val != 0, ]
  five_coefs[[i]] <- all_coefs[[i]][abs(all_coefs[[i]]$val) >= sort(abs(all_coefs[[i]]$val), decreasing = T)[ifelse(nrow(all_coefs[[i]]) >= 5, 5, nrow(all_coefs[[i]]))] & all_coefs[[i]]$val != 0, ]
  three_coefs[[i]] <- all_coefs[[i]][abs(all_coefs[[i]]$val) >= sort(abs(all_coefs[[i]]$val), decreasing = T)[ifelse(nrow(all_coefs[[i]]) >= 3, 3, nrow(all_coefs[[i]]))] & all_coefs[[i]]$val != 0, ]
}

all_coefs <- do.call("rbind", all_coefs)
ten_coefs <- do.call("rbind", ten_coefs)
five_coefs <- do.call("rbind", five_coefs)
three_coefs <- do.call("rbind", three_coefs)

###################################################
#"SURVEY_X20110.0.0.1" = Severe Depression
#"SURVEY_X20110.0.0.1.1" = Severe Depression
#"SURVEY_X20111.0.0.6" = Unknown
#SURVEY_X20086.0.0 = Vegan
#SURVEY_X20086.0.0.6 = Unknown
#SURVEY_X20110.0.0.1 = Severe Depression
#SURVEY_X20110.0.0.1.1 = Unknown



meds <- unique(grep("MEDS", ten_coefs$coef, value = T))
dict <- as.data.frame(fread("~/Downloads/coding4.tsv"))
new_meds <- better_swap_names(stringr::str_split(meds, "_", simplify = T)[,2], dict)
#write.table(cbind(meds, stringr::str_to_title(new_meds)), "helper_files/med_dict.txt", row.names = F, col.names = F, sep = "\t", quote = F)
meds_dict <- read.table("helper_files/med_dict.txt", stringsAsFactors = F, sep = "\t")
meds_dict[,2] <- paste0("Meds: ", meds_dict[,2])
ten_coefs$coef <- better_swap_names(ten_coefs$coef, meds_dict)

opcs <- unique(grep("OPCS", ten_coefs$coef, value = T))
dict <- as.data.frame(fread("~/Downloads/coding240.tsv"))
new_opcs <- better_swap_names(stringr::str_split(opcs, "_", simplify = T)[,2], dict)
#write.table(cbind(opcs, stringr::str_to_title(new_opcs)), "helper_files/opcs_dict.txt", row.names = F, col.names = F, sep = "\t", quote = F)
opcs_dict <- read.table("helper_files/opcs_dict.txt", stringsAsFactors = F, sep = "\t")
opcs_dict[,2] <- paste0("OPCS: ", opcs_dict[,2])
ten_coefs$coef <- better_swap_names(ten_coefs$coef, opcs_dict)

icd <- unique(grep("ICD", ten_coefs$coef, value = T))
dict <- as.data.frame(fread("~/Downloads/coding19.tsv"))
new_icd <- better_swap_names(stringr::str_split(icd, "_", simplify = T)[,2], dict)
#write.table(cbind(icd, stringr::str_to_title(new_icd)), "helper_files/icd_dict.txt", row.names = F, col.names = F, sep = "\t", quote = F)
icd_dict <- read.table("helper_files/icd_dict.txt", stringsAsFactors = F, sep = "\t")
icd_dict[,2] <- paste0("ICD: ", icd_dict[,2])
ten_coefs$coef <- better_swap_names(ten_coefs$coef, icd_dict)

survey <- unique(grep("SURVEY", ten_coefs$coef, value = T))
survey <- unlist(lapply(strsplit(survey, ".", fixed = T), function(x) x[1]))
ten_coefs$coef[grep("SURVEY", ten_coefs$coef)] <- unlist(lapply(strsplit(grep("SURVEY", ten_coefs$coef, value = T), ".", fixed = T), function(x) x[1]))
survey <- unique(survey)
#write.table(survey, "helper_files/survey_dict.txt", row.names = F, col.names = F, quote = F)
survey_dict <- read.table("helper_files/survey_dict.txt", stringsAsFactors = F, sep = "\t")
ten_coefs$coef <- better_swap_names(ten_coefs$coef, survey_dict)

census <- unique(grep("CENSUS", ten_coefs$coef, value = T))
#write.table(census, "helper_files/census_dict.txt", row.names = F, quote = F, col.names = F)
census_dict <- read.table("helper_files/census_dict.txt", stringsAsFactors = F, sep = "\t")
census_dict[,2] <- paste0("Census: ", census_dict[,2])
ten_coefs$coef <- better_swap_names(ten_coefs$coef, census_dict)




####################################################



ten_coefs$type <- trimws(stringr::str_split(ten_coefs$coef, ":", simplify = T)[,1])
#ten_coefs <- ten_coefs[ten_coefs$type != "eid",]
# ten_coefs$good <- ""
ten_coefs$for_box <- ""



# ten_coefs$good[ten_coefs$coef == "Survey: Other Bread Type"] <- "Survey: Other Bread Type"
# ten_coefs$good[ten_coefs$coef == "Meds: Salbutamol" & ten_coefs$disease == "Asthma"] <- "Meds: Salbutamol"
# ten_coefs$good[ten_coefs$coef == "Meds: Warfarin"] <- "Meds: Warfarin"
# ten_coefs$good[ten_coefs$coef == "Survey: Number Days Walked 10+ Min"] <- "Survey: Number Days Walked 10+ Min"
# ten_coefs$good[ten_coefs$coef == "Biomarker: Urate" & ten_coefs$disease == "Gout"] <- "Biomarker: Urate"
# ten_coefs$good[ten_coefs$coef == "Meds: Insulin" & ten_coefs$disease == "Type 1 Diabetes"] <- "Meds: Insulin"
# ten_coefs$good[ten_coefs$coef == "Survey: Ease of Getting Up In Morning"] <- "Survey: Ease of Getting Up In Morning"



ten_coefs$for_box[ten_coefs$coef == "Survey: Other Bread Type" & ten_coefs$disease == "Celiac Disease"] <- "Celiac Disease:\nOther Bread Type"
ten_coefs$for_box[ten_coefs$coef == "Survey: Number Days Walked 10+ Min" & ten_coefs$disease == "Multiple Sclerosis"] <- "MS: Num. Days\nWalked 10+ Min."
ten_coefs$for_box[ten_coefs$coef == "Metric: Year of Birth" & ten_coefs$disease == "Prostate Cancer"] <- "Prostate Cancer: Age"
#ten_coefs$for_box[ten_coefs$coef == "Meds: Warfarin" & ten_coefs$disease == "Atrial Fibrillation"] <- "Atrial Fibrillation:\nWarfarin"
ten_coefs$for_box[ten_coefs$coef == "Survey: Ease of Getting Up In Morning" & ten_coefs$disease == "Depression"] <- "Depression:Ease of\nGettin Up In Morning"
ten_coefs$for_box[ten_coefs$coef == "Biomarker: Urate" & ten_coefs$disease == "Gout"] <- "Gout: Urate"
ten_coefs$for_box[ten_coefs$coef == "Survey: Mother Diag Breast Cancer" & ten_coefs$disease == "Ovarian Cancer"] <- "Ovarian Cancer: Mother\nDiag. Breast Cancer"
ten_coefs$for_box[ten_coefs$coef == "Metric: Systolic Blood Pressure" & ten_coefs$disease == "Stroke"] <- "Stroke: Blood\nPressure"
ten_coefs$for_box[ten_coefs$coef == "OPCS: Diagnostic Spinal Puncture" & ten_coefs$disease == "Multiple Sclerosis"] <- "Multiple Sclerosis:\nSpinal Puncture"
ten_coefs$for_box[ten_coefs$coef == "Census: Other Qualifications" & ten_coefs$disease == "Hypertension"] <- "Hypertension: Other\nQualifications\n"
ten_coefs$for_box[ten_coefs$coef == "Meds: Carbimazole" & ten_coefs$disease == "Hyperthyroidism"] <- "Hyperthyroidism:\nCarbimazole"


ten_coefs$dup <- paste0(ten_coefs$disease, "_", ten_coefs$coef)
ten_coefs <- ten_coefs[!duplicated(ten_coefs$dup),]
ten_coefs$type <- unlist(lapply(strsplit(ten_coefs$coef, ":"), function(x) x[1]))


the_plot <- ggplot(ten_coefs, aes(val, disease)) + geom_vline(xintercept = 0) + geom_point(aes(color = type)) + 
  geom_hline(yintercept = 1:length(unique(ten_coefs$disease)) + 0.5, color = "grey80") +
  labs(x = "Model Coefficient", y = "", color = "Coef.\nType") 
plot(the_plot)
ggsave("out_plots/meta/python/coef_point.png", the_plot, width = 5, height = 5)


ten_coefs$type <- factor(as.character(ten_coefs$type), levels = c("Census", "Survey", "ICD", "Biomarker", "Meds", "OPCS", "Metric"))
the_plot <- ggplot(ten_coefs, aes(type, val, label = for_box)) + geom_boxplot() + geom_point() + 
  geom_text_repel(box.padding = 0.8, min.segment.length = 0, max.overlaps = Inf, nudge_y = 0.1, size = 3) +
  labs(x = "Coefficient Type", y = "Coefficient Value")
plot(the_plot)
ggsave("out_plots/meta/python/coef_box.png", the_plot, width = 8, height = 5)



#############################################################

conc <- readRDS("python.adjustment.RDS")
conc$author <- swap_names(conc$author, disease_names)
conc$author <- factor(conc$author, levels = conc$author[order(conc$concordance)])

the_plot <- ggplot(conc, aes(concordance, author)) + geom_point() +
  geom_errorbarh(aes(xmin = concordance - std, xmax = concordance + std), height = 0) +
  geom_hline(yintercept = 1:length(unique(conc$author)) + 0.5, color = "grey80") +
  labs(x = "Concordance", y = "")
plot(the_plot)
ggsave("out_plots/meta/python/conc.png", the_plot, width = 5, height = 5.5)

#coef_translator <- read.table("helper_files/coef_translate.txt", stringsAsFactors = F, sep = "\t")

# ten_coefs$coef <- better_swap_names(ten_coefs$coef, coef_translator)
# all_coefs$coef <- better_swap_names(all_coefs$coef, coef_translator)
# five_coefs$coef <- better_swap_names(five_coefs$coef, coef_translator)
# three_coefs$coef <- better_swap_names(three_coefs$coef, coef_translator)

# ten_coefs <- ten_coefs[ten_coefs$coef %in% coef_translator[,2],]
# 
# ggplot(ten_coefs, aes(val, coef, color = disease)) + geom_point()
