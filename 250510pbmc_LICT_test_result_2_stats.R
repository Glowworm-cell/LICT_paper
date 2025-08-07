setwd('D:/R/LICT_paper/pbmc2/')
library(readxl)
df11 = read.csv('./LICT_result_1.csv')
df12 = read.csv('./LICT_result_2.csv')
df1 = rbind(df11,df12)
df2 = read.csv('./avg_log2FC_res_step2_v1.csv')
df3 = read_excel('./pn_gene.xlsx')
df1 = cbind(df1, df2[,34:45])
df1$Free_LLMs_score = mapply(max, df1$Gemini_score, df1$Llama_score)
colnames(df1)[3] = 'ERNIE_annotation'
colnames(df1)[8] = 'Gemini_annotation'
colnames(df1)[13] = 'GPT_annotation'
colnames(df1)[18] = 'Llama_annotation'
colnames(df1)[23] = 'Claude_annotation'


pg = df3[2:nrow(df3),1:3]
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
df = df_count
df$manual_annotation = df2$manual_annotation
df$manual_positive_marker = pg$...2
df$manual_negative_marker = pg$...3
df$manual_reliable = rep('NO',nrow(df))
df$manual_positive_marker
for(i in 1:nrow(df_count)){
  n = df_count[i,]
  if(n$...2 >= 4){
    df[i,'manual_reliable'] = 'YES'
    }
}
df1 = cbind(df1,df[,4:7])
LICT_res = df1


df11 = read.csv('./LICT_feedback_result_1.csv')
df12 = read.csv('./LICT_feedback_result_2.csv')
df1 = rbind(df11,df12)
df2 = read.csv('./avg_log2FC_res_step2_v1_inter1.csv')
df3 = read_excel('./pn_gene.xlsx')
df1 = cbind(df1, df2[,34:45])
df1$Free_LLMs_score = mapply(max, df1$Gemini_score, df1$Llama_score)

df1[1,36:48]
df1[1,c('ERNIE_score','ERNIE_match')][1] == 1
df1[2,c('Llama_score','Llama_match')] = LICT_res[2,c('Llama_score','Llama_match')]

for (i in 1:nrow(df1)){
  n = df1[i,c('ERNIE_score','ERNIE_match')]
  if (n[1] >= LICT_res[i,c('ERNIE_score','ERNIE_match')][1]){
    df1[i,c('ERNIE_score','ERNIE_match')] = df1[i,c('ERNIE_score','ERNIE_match')]
  }else{
    df1[i,c('ERNIE_score','ERNIE_match')] = LICT_res[i,c('ERNIE_score','ERNIE_match')]
  }
  n = df1[i,c('Gemini_score','Gemini_match')]
  if (n[1] >= LICT_res[i,c('Gemini_score','Gemini_match')][1]){
    df1[i,c('Gemini_score','Gemini_match')] = df1[i,c('Gemini_score','Gemini_match')]
  }else{
    df1[i,c('Gemini_score','Gemini_match')] = LICT_res[i,c('Gemini_score','Gemini_match')]
  }
  n = df1[i,c('ChatGPT_score','ChatGPT_match')]
  if (n[1] >= LICT_res[i,c('ChatGPT_score','ChatGPT_match')][1]){
    df1[i,c('ChatGPT_score','ChatGPT_match')] = df1[i,c('ChatGPT_score','ChatGPT_match')]
  }else{
    df1[i,c('ChatGPT_score','ChatGPT_match')] = LICT_res[i,c('ChatGPT_score','ChatGPT_match')]
  }
  n = df1[i,c('Llama_score','Llama_match')]
  if (n[1] >= LICT_res[i,c('Llama_score','Llama_match')][1]){
    df1[i,c('Llama_score','Llama_match')] = df1[i,c('Llama_score','Llama_match')]
  }else{
    df1[i,c('Llama_score','Llama_match')] = LICT_res[i,c('Llama_score','Llama_match')]
  }
  n = df1[i,c('Claude_score','Claude_match')]
  if (n[1] >= LICT_res[i,c('Claude_score','Claude_match')][1]){
    df1[i,c('Claude_score','Claude_match')] = df1[i,c('Claude_score','Claude_match')]
  }else{
    df1[i,c('Claude_score','Claude_match')] = LICT_res[i,c('Claude_score','Claude_match')]
  }
  n = df1[i,c('max_score','max_match')]
  if (n[1] >= LICT_res[i,c('max_score','max_match')][1]){
    df1[i,c('max_score','max_match')] = df1[i,c('max_score','max_match')]
  }else{
    df1[i,c('max_score','max_match')] = LICT_res[i,c('max_score','max_match')]
  }
  n = df1[i,c('Free_LLMs_score')]
  if (n[1] >= LICT_res[i,c('Free_LLMs_score')][1]){
    df1[i,c('Free_LLMs_score')] = df1[i,c('Free_LLMs_score')]
  }else{
    df1[i,c('Free_LLMs_score')] = LICT_res[i,c('Free_LLMs_score')]
  }
}
colnames(df1)[3] = 'ERNIE_annotation'
colnames(df1)[8] = 'Gemini_annotation'
colnames(df1)[13] = 'GPT_annotation'
colnames(df1)[18] = 'Llama_annotation'
colnames(df1)[23] = 'Claude_annotation'


pg = df3[2:nrow(df3),1:3]
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
df = df_count
df$manual_annotation = df2$manual_annotation
df$manual_positive_marker = pg$...2
df$manual_negative_marker = pg$...3
df$manual_reliable = rep('NO',nrow(df))
df$manual_positive_marker
for(i in 1:nrow(df_count)){
  n = df_count[i,]
  if(n$...2 >= 4){
    df[i,'manual_reliable'] = 'YES'
  }
}
df1 = cbind(df1,df[,4:7])
LICT_feedback_res = df1


stats = data.frame('all reliable'= c(0,0,0,0),
                   'only LLMs reliable'= c(0,0,0,0),
                   'only manual reliable' = c(0,0,0,0),
                   'all unreliable'= c(0,0,0,0))
rownames(stats) = c('all LLMs','Free LLMs','all LLMs with feedback','Free LLMs with feedback')
###row1
for(i in 1:nrow(LICT_res)){
  if(LICT_res[i,]$Total_reliable == 'YES'){
    if(LICT_res[i,]$manual_reliable == 'YES'){
      stats[1,]$all.reliable = stats[1,]$all.reliable+1
    }else{
      stats[1,]$only.LLMs.reliable = stats[1,]$only.LLMs.reliable+1
    }
  }else{
    if(LICT_res[i,]$manual_reliable == 'YES'){
      stats[1,]$only.manual.reliable = stats[1,]$only.manual.reliable+1
    }else{
      stats[1,]$all.unreliable = stats[1,]$all.unreliable+1
    }
  }
}
###row2
for(i in 1:nrow(LICT_res)){
  if(LICT_res[i,]$Gemini_reliable == 'YES' | LICT_res[i,]$Llama_reliable == 'YES'){
    if(LICT_res[i,]$manual_reliable == 'YES'){
      stats[2,]$all.reliable = stats[2,]$all.reliable+1
    }else{
      stats[2,]$only.LLMs.reliable = stats[2,]$only.LLMs.reliable+1
    }
  }else{
    if(LICT_res[i,]$Gemini_reliable == 'YES' | LICT_res[i,]$Llama_reliable == 'YES'){
      stats[2,]$only.manual.reliable = stats[2,]$only.manual.reliable+1
    }else{
      stats[2,]$all.unreliable = stats[2,]$all.unreliable+1
    }
  }
}
###row3
for(i in 1:nrow(LICT_feedback_res)){
  if(LICT_feedback_res[i,]$Total_reliable == 'YES'){
    if(LICT_feedback_res[i,]$manual_reliable == 'YES'){
      stats[3,]$all.reliable = stats[3,]$all.reliable+1
    }else{
      stats[3,]$only.LLMs.reliable = stats[3,]$only.LLMs.reliable+1
    }
  }else{
    if(LICT_feedback_res[i,]$manual_reliable == 'YES'){
      stats[3,]$only.manual.reliable = stats[3,]$only.manual.reliable+1
    }else{
      stats[3,]$all.unreliable = stats[3,]$all.unreliable+1
    }
  }
}
###row4
for(i in 1:nrow(LICT_feedback_res)){
  if(LICT_feedback_res[i,]$Gemini_reliable == 'YES' | LICT_feedback_res[i,]$Llama_reliable == 'YES'){
    if(LICT_feedback_res[i,]$manual_reliable == 'YES'){
      stats[4,]$all.reliable = stats[4,]$all.reliable+1
    }else{
      stats[4,]$only.LLMs.reliable = stats[4,]$only.LLMs.reliable+1
    }
  }else{
    if(LICT_feedback_res[i,]$Gemini_reliable == 'YES' | LICT_feedback_res[i,]$Llama_reliable == 'YES'){
      stats[4,]$only.manual.reliable = stats[4,]$only.manual.reliable+1
    }else{
      stats[4,]$all.unreliable = stats[4,]$all.unreliable+1
    }
  }
}
write.csv(LICT_res,'./LICT_res.csv')
write.csv(LICT_feedback_res,'./LICT_feedback_res.csv')
write.csv(stats,'./stats.csv')

#subset mismatch
mismatch_LICT_res = subset(LICT_res, max_score == 0)
mismatch_LICT_feedback_res = subset(LICT_feedback_res, max_score == 0)

mismatch_stats = data.frame('all reliable'= c(0,0,0,0),
                            'only LLMs reliable'= c(0,0,0,0),
                            'only manual reliable' = c(0,0,0,0),
                            'all unreliable'= c(0,0,0,0))
rownames(mismatch_stats) = c('all LLMs','Free LLMs','all LLMs with feedback','Free LLMs with feedback')
###row1
for(i in 1:nrow(mismatch_LICT_res)){
  if(mismatch_LICT_res[i,]$Total_reliable == 'YES'){
    if(mismatch_LICT_res[i,]$manual_reliable == 'YES'){
      mismatch_stats[1,]$all.reliable = mismatch_stats[1,]$all.reliable+1
    }else{
      mismatch_stats[1,]$only.LLMs.reliable = mismatch_stats[1,]$only.LLMs.reliable+1
    }
  }else{
    if(mismatch_LICT_res[i,]$manual_reliable == 'YES'){
      mismatch_stats[1,]$only.manual.reliable = mismatch_stats[1,]$only.manual.reliable+1
    }else{
      mismatch_stats[1,]$all.unreliable = mismatch_stats[1,]$all.unreliable+1
    }
  }
}
###row2
for(i in 1:nrow(mismatch_LICT_res)){
  if(mismatch_LICT_res[i,]$Gemini_reliable == 'YES' | mismatch_LICT_res[i,]$Llama_reliable == 'YES'){
    if(mismatch_LICT_res[i,]$manual_reliable == 'YES'){
      mismatch_stats[2,]$all.reliable = mismatch_stats[2,]$all.reliable+1
    }else{
      mismatch_stats[2,]$only.LLMs.reliable = mismatch_stats[2,]$only.LLMs.reliable+1
    }
  }else{
    if(mismatch_LICT_res[i,]$Gemini_reliable == 'YES' | mismatch_LICT_res[i,]$Llama_reliable == 'YES'){
      mismatch_stats[2,]$only.manual.reliable = mismatch_stats[2,]$only.manual.reliable+1
    }else{
      mismatch_stats[2,]$all.unreliable = mismatch_stats[2,]$all.unreliable+1
    }
  }
}
###row3
for(i in 1:nrow(mismatch_LICT_feedback_res)){
  if(mismatch_LICT_feedback_res[i,]$Total_reliable == 'YES'){
    if(mismatch_LICT_feedback_res[i,]$manual_reliable == 'YES'){
      mismatch_stats[3,]$all.reliable = mismatch_stats[3,]$all.reliable+1
    }else{
      mismatch_stats[3,]$only.LLMs.reliable = mismatch_stats[3,]$only.LLMs.reliable+1
    }
  }else{
    if(mismatch_LICT_feedback_res[i,]$manual_reliable == 'YES'){
      mismatch_stats[3,]$only.manual.reliable = mismatch_stats[3,]$only.manual.reliable+1
    }else{
      mismatch_stats[3,]$all.unreliable = mismatch_stats[3,]$all.unreliable+1
    }
  }
}
###row4
for(i in 1:nrow(mismatch_LICT_feedback_res)){
  if(mismatch_LICT_feedback_res[i,]$Gemini_reliable == 'YES' | mismatch_LICT_feedback_res[i,]$Llama_reliable == 'YES'){
    if(mismatch_LICT_feedback_res[i,]$manual_reliable == 'YES'){
      mismatch_stats[4,]$all.reliable = mismatch_stats[4,]$all.reliable+1
    }else{
      mismatch_stats[4,]$only.LLMs.reliable = mismatch_stats[4,]$only.LLMs.reliable+1
    }
  }else{
    if(mismatch_LICT_feedback_res[i,]$Gemini_reliable == 'YES' | mismatch_LICT_feedback_res[i,]$Llama_reliable == 'YES'){
      mismatch_stats[4,]$only.manual.reliable = mismatch_stats[4,]$only.manual.reliable+1
    }else{
      mismatch_stats[4,]$all.unreliable = mismatch_stats[4,]$all.unreliable+1
    }
  }
}
write.csv(mismatch_LICT_res,'./mismatch_LICT_res.csv')
write.csv(mismatch_LICT_feedback_res,'./mismatch_LICT_feedback_res.csv')
write.csv(mismatch_stats,'./mismatch_stats.csv')

