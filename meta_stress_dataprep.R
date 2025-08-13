################### Data preparation ##################################################
#---- 
# Clear global environment and load required packages
rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, knitr, MAd, metafor, reactable, readr, readxl, rmarkdown, 
               R.rsp, stringr, tidyr) 

#setwd("Y:/Paper/Drafts/pulsedtaVNS_pupil_meta/code")
#setwd("D:/SynologyDrive/Paper/Drafts/pulsedtaVNS_pupil_meta/code")
#setwd("/Users/lisa/Desktop/code_tVNS_pupil_meta")

# Import Master xlsx file, in which articles were coded
Master <- read_excel("/mnt/big_data/SynologyDrive/Paper/drafts/Review_Stress_metabolism/Analysis/master_table.xlsx",sheet = "Data", col_types = c("text", 
                                                                                                                                                  "text", "text", "text", "text", "text", 
                                                                                                                                                  "text", "text", "text", "text", "text", 
                                                                                                                                                  "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                                                                                                                  "numeric", "numeric", "numeric", "text", "text",
                                                                                                                                                  "text", "text", "numeric", "text", "text", 
                                                                                                                                                  "text", "text", "numeric", "text", "numeric", 
                                                                                                                                                  "text", "text", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                                                                                                                  "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                                                                                                                  "numeric", "numeric", "text", "text", "numeric", "numeric","numeric",
                                                                                                                                                  "numeric", "numeric", "numeric", "numeric", "text", 
                                                                                                                                                  "numeric", "numeric", "numeric", "numeric", "numeric", "text",
                                                                                                                                                  "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                                                                                                                  "text", "numeric", "numeric", "numeric", "numeric", 
                                                                                                                                                  "numeric", "text", "numeric", "numeric", 
                                                                                                                                                  "numeric", "text"), na = "NA")

# Generate column "study" serving as publication identifier (last name of 1st author + publication year)
firstauthor <- word(string = Master$Author, sep = ",") #Grab last name of the 1st author
study<-as.data.frame(cbind(firstauthor,Master$Publication_Year)) #Dataframe with 1st authors last name and publication year
study <- study %>% unite(study, sep = ", ") #Unite it into 1 column
Master <- cbind(study,Master)

# # Table showing reasons for exclusions
 excluded <- Master %>% filter(Master$Exclusion == 1)
 excluded <- unique(setDT(excluded) [sort.list(study)], by = "DOI")
 reasons_exc <- table(excluded$Reason.for.exclusion)

# Generate data frame "df" consisting of all studies 
Master$Exclusion[is.na(Master$Exclusion)] <- 0
Master$`Sex diff`[is.na(Master$`Sex diff`)] <- 0
df <- Master%>% filter(Master$Exclusion!=1)

# Calculate cort increase for groups
df$cort_change_low <- if_else(!is.na(df$`AUC Cortisol_low_hormone`),df$`AUC Cortisol_low_hormone`,df$`Peak mean_cortisol_low_hormone` - df$`Stress baseline mean_cortisol_low_hormone`)
df$cort_change_high <- if_else(!is.na(df$`AUC Cortisol_high_hormone`),df$`AUC Cortisol_high_hormone`,df$`Peak mean_cortisol_high_hormone`-df$`Stress onset mean_cortisol_high_hormone`)
df$cort_change_sd_low <- if_else(!is.na(df$`AUC_SD Cortisol_low_hormone`),df$`AUC_SD Cortisol_low_hormone`,df$`Peak SD_cortisol_low_hormone`)
df$cort_change_sd_high <- if_else(!is.na(df$`AUC_SD Cortisol_high_hormone`),df$`AUC_SD Cortisol_high_hormone`,df$`Peak SD_cortisol_high_hormone`)



# # Calculate effect sizes (Hedges' g) and store them in "es"
es <- as.data.frame(es<-escalc(m1i = df$cort_change_high,
                               sd1i =df$cort_change_sd_high,
                               n1i =df$N_low,
                               m2i = df$cort_change_low,
                               sd2i = df$cort_change_sd_low,
                               n2i = df$N_low,
                               measure = "SMD",
                               yi, sei, vi, correct=FALSE))

### Manual change!
# Convert Correlation Coefficient r to Effect Size d

# Estimation of Standard Error and Variance
se <- (1 - df$correlation_cortisol**2) / sqrt(df$N_high - 2)
se
var <- (1 - (df$correlation_cortisol ** 2)) **2 / (df$N_high - 2)
var
es_d <- res(df$correlation_cortisol, var.r = var, df$N_high)
#result_es <- hedges_g(d = es_d['d'], totaln = df$N_high)
#result_es

es <- cbind(es,es_d[,2:3])

es[is.na(es[,1]),1] <- es[is.na(es[,1]),3]
es[is.na(es[,2]),2] <- es[is.na(es[,2]),4]



df<-cbind(es[,1:2],df)



########################################

### Manual change!
df$Publication_Year <- as.numeric(df$Publication_Year)

# #Adds a ' in front of numbers to prevent transformation into date format in Excel. Make sure to remove before using output!!!
#df$Mean.Age<-paste0("'",df$Mean.Age)


df <-  df[complete.cases(df$yi), ] 
df <- df[complete.cases(df$vi), ] 


# Add "https://doi.org/" infront of dois
Master$DOI <- paste0("https://doi.org/", Master$DOI)
Master$DOI <- paste0("<a href='", Master$DOI,"'target='_blank'>",Master$DOI,"</a>")

# Write "df" as .csv
write.csv2(df, "df.csv")

# Save "df" as .RDa
save(df, file = "df.RDa")

# Save overview of all screened studies
screened <- Master[,c(1,7,9,14:15)]
screened$Included <- factor(screened$Included,levels = c(0,1), labels = c("no","yes"))
screened <- unique(setDT(screened) [sort.list(study)], by = "DOI")
colnames(screened) <- c("Study", "Title", "DOI", "Design", "Sample")
save(screened, file = "screened.RDa")

