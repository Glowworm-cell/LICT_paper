library(readxl)
library(stringr)
setwd('D:/R/LLMsCellType/240616gc_different_method6/')
pg = read_excel('../240618Mouse_SS_Fibro_different_method/pn_gene - 副本 (2).xlsx')
head(pg)
df_split <- pg
df_split[] <- lapply(df_split, function(x) {
  if (is.character(x)) {
    strsplit(x, ",")
  } else {
    x
  }
})
df_count <- df_split
df_count[] <- lapply(df_count, function(x) {
  if (is.list(x)) {
    sapply(x, length)  # 计算每个列表中元素的数量
  } else {
    x
  }
})
df = df_count[-1,]
df$all_valid = rep(0,nrow(df))
df$model_valid = rep(0,nrow(df))
df$manual_valid = rep(0,nrow(df))
df$all_unvalid = rep(0,nrow(df))


for(i in 1:nrow(df)){
  n = df[i,]
  if(n$...2 >= 4){
    if(n$...11>=4 | n$...14>=4){
      df[i,"all_valid"] = 1
    }else{
      df[i,"manual_valid"] = 1
    }
  }else{
    if(n$...11>=4 | n$...14>=4){
      df[i,"model_valid"] = 1
    }else{
      df[i,"all_unvalid"] = 1
    }
  }
}

pbmc1 = df[,19:22]

pg = read_excel('../240618Mouse_SS_Fibro_different_method_second time/pn_gene - 副本 (2).xlsx')
head(pg)
df_split <- pg
df_split[] <- lapply(df_split, function(x) {
  if (is.character(x)) {
    strsplit(x, ",")
  } else {
    x
  }
})
df_count <- df_split
df_count[] <- lapply(df_count, function(x) {
  if (is.list(x)) {
    sapply(x, length)  # 计算每个列表中元素的数量
  } else {
    x
  }
})
df = df_count[-1,]
df$all_valid = rep(0,nrow(df))
df$model_valid = rep(0,nrow(df))
df$manual_valid = rep(0,nrow(df))
df$all_unvalid = rep(0,nrow(df))


for(i in 1:nrow(df)){
  n = df[i,]
  if(n$...2 >= 4){
    if(n$...11>=4 | n$...14>=4){
      df[i,"all_valid"] = 1
    }else{
      df[i,"manual_valid"] = 1
    }
  }else{
    if(n$...11>=4 | n$...14>=4){
      df[i,"model_valid"] = 1
    }else{
      df[i,"all_unvalid"] = 1
    }
  }
}

pbmc2 = df[,19:22]

pg = read_excel('../240618Mouse_SS_Fibro_different_method_third time/pn_gene - 副本 (2).xlsx')
head(pg)
df_split <- pg
df_split[] <- lapply(df_split, function(x) {
  if (is.character(x)) {
    strsplit(x, ",")
  } else {
    x
  }
})
df_count <- df_split
df_count[] <- lapply(df_count, function(x) {
  if (is.list(x)) {
    sapply(x, length)  # 计算每个列表中元素的数量
  } else {
    x
  }
})
df = df_count[-1,]
df$all_valid = rep(0,nrow(df))
df$model_valid = rep(0,nrow(df))
df$manual_valid = rep(0,nrow(df))
df$all_unvalid = rep(0,nrow(df))


for(i in 1:nrow(df)){
  n = df[i,]
  if(n$...2 >= 4){
    if(n$...11>=4 | n$...14>=4){
      df[i,"all_valid"] = 1
    }else{
      df[i,"manual_valid"] = 1
    }
  }else{
    if(n$...11>=4 | n$...14>=4){
      df[i,"model_valid"] = 1
    }else{
      df[i,"all_unvalid"] = 1
    }
  }
}

pbmc3 = df[,19:22]

pbmc = rbind(pbmc1,pbmc2,pbmc3)
pbmc_sum = colSums(pbmc)
fibro_sum = pbmc_sum
write.csv(fibro_sum,'./fibro_sum_valid.csv')
