library(data.table)
library(dplyr)

path <- "V:/ddata/CELB/poot/Eline Koornstra/SuRE_general/Results/freeze7/statistics/"
file <- "hnsc_freeze7_for-fisher_10042024.txt"
outname <- "hnsc_freeze7_locations_statistics_10042024.txt"

df <- fread(paste0(path, file))
df$FE_emvarvscontrol <- NA
df$Fisher_OR_emvarvscontrol <- NA
df$Fisher_p_emvarvscontrol <- NA
#df$chisq_p_emvarvscontrol <- NA

df$FE_emvarvsall <- NA
df$Fisher_OR_emvarvsall <- NA
df$Fisher_p_emvarvsall <- NA
#df$chisq_p_emvarvsall <- NA

for (row in 1:nrow(df)) {
  emvar <- df$n_emvar[row]
  no_emvar <- df$emvar_total[row] - df$n_emvar[row]
  control <- df$n_control[row]
  no_control <- df$control_total[row] - df$n_control[row]
  
  df$FE_emvarvscontrol[row] <- (df$n_emvar[row]/df$emvar_total[row])/(df$n_control[row]/df$control_total[row])
  
  a <- data.frame(
    "hit" = c(emvar, control),
    "no_hit" = c(no_emvar, no_control),
    row.names = c("emvar", "control"),
    stringsAsFactors = FALSE
  )
  colnames(a) <- c("hit", "no hit")
  
  a
  
  mosaicplot(a, main = "Mosaic plot", color = TRUE)
  
  test <- fisher.test(a)
  test
  test$p.value
  df$Fisher_OR_emvarvscontrol[row] <- test$estimate
  df$Fisher_p_emvarvscontrol[row] <- test$p.value
 
  
  #test2 <- chisq.test(a)
  #test2
  #test2$p.value
  #df$chisq_p_emvarvscontrol[row] <- test2$p.value
  
}

for (row in 1:nrow(df)) {
  emvar <- df$n_emvar[row]
  no_emvar <- df$emvar_total[row] - df$n_emvar[row]
  control <- df$n_all[row]
  no_control <- df$all_total[row] - df$n_all[row]
  
  df$FE_emvarvsall[row] <- (df$n_emvar[row]/df$emvar_total[row])/(df$n_all[row]/df$all_total[row])
  
  a <- data.frame(
    "hit" = c(emvar, control),
    "no_hit" = c(no_emvar, no_control),
    row.names = c("emvar", "control"),
    stringsAsFactors = FALSE
  )
  colnames(a) <- c("hit", "no hit")
  
  a
  
  mosaicplot(a, main = "Mosaic plot", color = TRUE)
  
  test <- fisher.test(a)
  test
  test$p.value
  df$Fisher_OR_emvarvsall[row] <- test$estimate
  df$Fisher_p_emvarvsall[row] <- test$p.value
  
  
  #test2 <- chisq.test(a)
  #test2
  #test2$p.value
  #df$chisq_p_emvarvsall[row] <- test2$p.value
  
}

df

fwrite(df, file = paste0(path, outname), row.names = F, quote = F, sep = '\t')


