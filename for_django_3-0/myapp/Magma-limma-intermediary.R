library(limma)
library(data.table)

# Read the input filename from the command line argument
args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]
input_filename <- paste(args[1],'.txt',sep = "")
df_final = read.delim(input_filename)
head(df_final)

# setwd("/Users/user/Downloads/TMT_Analysis/Paper_analysis/SILAC/")
# df_final = read.delim("METTL3_SILAC_Will.txt")

# df_final = read.delim("/Users/user/Downloads/TMT_IPMS_Server/src/for_django_3-0/media/documents/tmtipms_files/four-one_prelimma.txt")

columns_check = grepl('rep[0-9]',colnames(df_final))
# Code for FC and pvalue calculation


# # Code if unmoderated one sample t-test
# myfit <- limma::lmFit(df_final[,columns_check], method="robust", maxit=1000)
# logFC <- rowMeans(df_final[,columns_check],na.rm=TRUE)
# unmod.t <- myfit$coefficients/myfit$stdev.unscaled/myfit$sigma
# pvalue <- 2*pt(-abs(unmod.t), myfit$df.residual)
# FDR <- p.adjust(pvalue, method="BH")
# result <- data.frame(cbind(df_final, logFC, pvalue, FDR))

# library(broom)
# #Code if running student t-test one sample
# logFC <- rowMeans(df_final[,columns_check], na.rm=TRUE)
# sd_all <- apply(df_final[,columns_check],1,function(i)sd(i, na.rm = TRUE))
# sample_size <- apply(df_final[,columns_check],1,function(i)sum(!is.na(i)))
# pvalue <- rep(NA, nrow(df_final))
# for (i in 1:nrow(df_final)){
#   # pvalue[i]<-NA
#   if ((sd_all[i]>0)&!(is.na(sd_all[i]))){
#     tstat <- (logFC[i] - 0)/(sd_all[i]/sqrt(sample_size[i]))
#     pvalue[i] <- 2*pt(-abs(tstat), sample_size[i]-1)
#   }
# }
# FDR <- p.adjust(pvalue, method="BH")
# sd_all <- NA
# sample_size <- NA
# original_dt <- as.data.table(df_final)
# original_dt[, logFC := logFC]
# original_dt[, pvalue := pvalue]
# original_dt[, FDR := FDR]
# # result <- data.frame(cbind(df_final, logFC, pvalue, FDR))
# result <- as.data.frame(original_dt)

#Code if moderated one sample t-test
myfit <- limma::lmFit(df_final[,columns_check], method="robust", maxit=1000)
efit <- limma::eBayes(myfit)
# Code to do get FDR (Benjamini Hochberg)
modtest <- limma::topTable(efit, number=nrow(myfit), sort.by='none')
colnames(modtest)[4:5] <- c("pvalue","FDR")
result <- data.frame(cbind(df_final, modtest[,-c(2,3,6)]))

head(result)

write.table(result, paste(filename,"_out.txt",sep = ""), sep="\t", row.names = FALSE)
