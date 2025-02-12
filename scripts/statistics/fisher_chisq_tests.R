library(data.table)
library(dplyr)

path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/statistics/"
file <- "hnsc_freeze7_for-fisher_10042024.txt"
outname <- "hnsc_freeze7_locations_statistics_10042024.txt"

df <- fread(paste0(path, file))
df$FE_raqtlvscontrol <- NA
df$Fisher_OR_raqtlvscontrol <- NA
df$Fisher_p_raqtlvscontrol <- NA
#df$chisq_p_raqtlvscontrol <- NA

df$FE_raqtlvsall <- NA
df$Fisher_OR_raqtlvsall <- NA
df$Fisher_p_raqtlvsall <- NA
#df$chisq_p_raqtlvsall <- NA

for (row in 1:nrow(df)) {
  raqtl <- df$n_raqtl[row]
  no_raqtl <- df$raqtl_total[row] - df$n_raqtl[row]
  control <- df$n_control[row]
  no_control <- df$control_total[row] - df$n_control[row]
  
  df$FE_raqtlvscontrol[row] <- (df$n_raqtl[row]/df$raqtl_total[row])/(df$n_control[row]/df$control_total[row])
  
  a <- data.frame(
    "hit" = c(raqtl, control),
    "no_hit" = c(no_raqtl, no_control),
    row.names = c("raQTL", "control"),
    stringsAsFactors = FALSE
  )
  colnames(a) <- c("hit", "no hit")
  
  a
  
  mosaicplot(a, main = "Mosaic plot", color = TRUE)
  
  test <- fisher.test(a)
  test
  test$p.value
  df$Fisher_OR_raqtlvscontrol[row] <- test$estimate
  df$Fisher_p_raqtlvscontrol[row] <- test$p.value
 
  
  #test2 <- chisq.test(a)
  #test2
  #test2$p.value
  #df$chisq_p_raqtlvscontrol[row] <- test2$p.value
  
}

for (row in 1:nrow(df)) {
  raqtl <- df$n_raqtl[row]
  no_raqtl <- df$raqtl_total[row] - df$n_raqtl[row]
  control <- df$n_all[row]
  no_control <- df$all_total[row] - df$n_all[row]
  
  df$FE_raqtlvsall[row] <- (df$n_raqtl[row]/df$raqtl_total[row])/(df$n_all[row]/df$all_total[row])
  
  a <- data.frame(
    "hit" = c(raqtl, control),
    "no_hit" = c(no_raqtl, no_control),
    row.names = c("raQTL", "control"),
    stringsAsFactors = FALSE
  )
  colnames(a) <- c("hit", "no hit")
  
  a
  
  mosaicplot(a, main = "Mosaic plot", color = TRUE)
  
  test <- fisher.test(a)
  test
  test$p.value
  df$Fisher_OR_raqtlvsall[row] <- test$estimate
  df$Fisher_p_raqtlvsall[row] <- test$p.value
  
  
  #test2 <- chisq.test(a)
  #test2
  #test2$p.value
  #df$chisq_p_raqtlvsall[row] <- test2$p.value
  
}

df

fwrite(df, file = paste0(path, outname), row.names = F, quote = F, sep = '\t')

