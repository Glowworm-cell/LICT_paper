setwd('D:/R/LLMsCellType/240616gc_different_method6/')
detach(package:LLMcellidentifier, unload=TRUE)
install.packages("D:/ywj/LLM514/LLMcellidentifier_0.1.0.tar.gz", repos = NULL, type = "source")
library(LLMcellidentifier)
library(reticulate)
library(rols)
library(dplyr)
library(scales)
path <- Sys.which("python")
Sys.setenv(RETICULATE_PYTHON = path)

reticulate::py_run_string("
import os
import openai
ERNIE_api_key = os.environ['ERNIE_api_key']
ERNIE_secret_key = os.environ['ERNIE_secret_key']
GEMINI_api_key = os.environ['GEMINI_api_key']
openai.api_key = os.environ['openai.api_key']
Llama3_api_key = os.environ['Llama3_api_key']
Llama3_secret_key = os.environ['Llama3_secret_key']
ANTHROPIC_API_KEY = os.environ['ANTHROPIC_API_KEY']
")
markers = readRDS('./gc_markers.rds')

gc_res = LLMCelltype(FindAllMarkersResult = markers,
                       species = 'human',
                       topgenenumber = 10,
                       tissuename = 'gastric tumor')

res = gc_res

create_dataframe <- function(res, markers) {
  top_markers <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
  cluster_genes <- top_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(genes = paste(gene, collapse=", ")) %>%
    dplyr::ungroup()
  df1 <- data.frame(
    clusters = res[[1]]$clusters,
    ERNIE = res[['ERNIE']]$cell_type,
    Gemini = res[['Gemini']]$cell_type,
    ChatGPT = res[['GPT']]$cell_type,
    Llama = res[['Llama']]$cell_type,
    Claude = res[['Claude']]$cell_type
  )
  df = cbind(cluster_genes,df1)
  return(df)
}
res1 = create_dataframe(res, markers)
res1$cluster = as.character(res1$cluster)
res1[11,1] = 'B cells'
res1[2,1] = 'Epithelial cells'
write.csv(res1,'./avg_log2FC_step1_res1.csv')
#res1 = read.csv('./avg_log2FC_step1_res1.csv')
#res1 <- res1[ , -1]

cl = res1[,1:2]
OlsSearch(q = 'axial mesoderm', ontology = "cl")
for(j in c('cluster','ERNIE','Gemini','ChatGPT','Llama','Claude')){
  query <- res1[[j]]
  cv <- sapply(query, function(i) {
    if(i == "undefined"|i == "undefined "){
      print('skip undefined')
      c(i, NA, NA)
    }else{
      print(i)
      new_string <- tolower(as.character(i))
      new_string <- sub("s$", "", new_string)
      d <- py$get_class_id_from_ontology("D:/R/LLMsCellType/cl.owl",new_string)
      if(d == c("no match")){
        d <- as(olsSearch(OlsSearch(q = new_string, ontology = "cl")), "data.frame")
        d <- d[grep('CL:', d$obo_id),]
        c(i, d[1, 'label'], d[1, 'obo_id'])
      }else{
        d = gsub("_", ":", d)
        c(i, new_string, d)
      }
    }
  })
  cv <- t(cv)
  cv <- as.data.frame(cv)
  colnames(cv) <- c(j,'CLname','CLID')
  cl = cbind(cl,cv)
}

get_ontology_parents <- function(class_id) {
  print(class_id)
  # 加载Python代码
  py_run_string("
from owlready2 import get_ontology, ThingClass

def get_all_ancestors(ont_class, depth=0):
    # 递归函数，用于获取所有父类的信息
    indent = '  ' * depth  # 用于格式化输出，增加缩进以表示层级关系
    ancestors_info = ''
    for parent in ont_class.is_a:
        if isinstance(parent, ThingClass):  # 确保只考虑本体类，而非约束等其他对象
            label = parent.label[0] if parent.label else 'No Label'
            ancestors_info += f'{indent}- {parent.name}\\n'
            ancestors_info += f'{indent}  Label: {label}\\n'
            ancestors_info += get_all_ancestors(parent, depth + 1)  # 递归调用以获取父类的父类
    return ancestors_info

def get_all_ancestors_of_class(ont_path, class_id):
    try:
        # 加载本体
        onto = get_ontology(ont_path).load()
        # 查找具有指定IRI的类
        specific_class = onto.search_one(iri=f'*{class_id}')
        if not specific_class:
            return 'Class not found.'
        # 构建类的基本信息字符串
        class_info = f'Class: {specific_class.name}\\n'
        if specific_class.label:
            class_info += f'Label: {specific_class.label[0]}\\n'
        # 使用递归函数获取所有父类的信息
        ancestors_info = 'All Ancestors:\\n' + get_all_ancestors(specific_class)
        return ancestors_info
    except Exception as e:
        return f'Failed to load the ontology or find the class: {e}'
")

  # 调用Python函数
  dp <- py$get_all_ancestors_of_class("D:/R/LLMsCellType/cl.owl", class_id)

  # 使用正则表达式删除所有 \n 后面的空格
  dp_cleaned <- gsub("\n\\s+", "\n", dp)

  lines <- str_split(dp_cleaned, "\n")[[1]]

  # 提取ID和标签
  ids <- str_extract(lines, "^\\s*-\\s*[^\\s]+")
  labels <- str_extract(lines, "(?<=Label: ).*")

  ids <- str_remove(ids, "^\\s*-\\s*")
  # 去掉NA值并将ID和标签配对
  ids_na = !is.na(ids)
  labels_na = !is.na(labels)
  id_label_pairs <- data.frame(ID = ids[ids_na], Label = labels[labels_na], stringsAsFactors = FALSE)

  # 创建一个新的数据框，仅包含有效的ID和标签，并去重
  df <- id_label_pairs %>%
    distinct()

  # 输入的ID列表
  excluded_ids <- c('CL_0000197', 'CL_0002321', 'CL_0000039', 'CL_0000413', 'CL_0000255',
                    'CL_0000412', 'CL_0011115', 'CL_0000988', 'CL_0002319', 'CL_0000151',
                    'CL_0002320', 'CL_0000064', 'CL_0000066', 'CL_0000215', 'CL_0000080',
                    'CL_0001035', 'CL_0000325', 'CL_4030031', 'CL_0000183', 'CL_0000188',
                    'CL_0000211', 'CL_0000212', 'CL_0000630', 'CL_0000219', 'CL_0000225',
                    'CL_0002242', 'CL_0000329', 'CL_0000473', 'CL_0000520', 'CL_0000293',
                    'CL_2000021', 'CL_0000349', 'CL_4033054', 'CL_0010017', 'CL_0000415',
                    'CL_0000422', 'CL_0000445', 'CL_0002494', 'CL_0000891', 'CL_0001034',
                    'CL_0000677', 'CL_0000627', 'CL_0000628', 'CL_0000725', 'CL_0001061',
                    'CL_0002559', 'CL_1000497', 'CL_0007001', 'CL_0007005', 'CL_0007022',
                    'CL_0007023', 'CL_0009002', 'CL_0009003', 'CL_0009005', 'CL_0009010',
                    'CL_1000600', 'CL_1000601', 'CL_2000020', 'CL_4029001', 'CL_0000000')

  # 假设 df 是已有的数据框
  # 过滤 ID 以 CL 开头并且不在排除列表中的行
  df_filtered <- df %>%
    filter(str_detect(ID, "^CL") & !ID %in% excluded_ids)

  # 合并ID和Label列
  merged_id <- paste(df_filtered$ID, collapse = ";")
  merged_label <- paste(df_filtered$Label, collapse = ";")

  # 创建新的数据框
  merged_df <- data.frame(ID = merged_id, Label = merged_label, stringsAsFactors = FALSE)

  return(merged_df)

}

res = cl[,1:2]
for(j in c(5,8,11,14,17,20)){
  df = cl[,c(j-2):j]
  print(colnames(df)[[1]])
  query = df[[3]]
  out <- sapply(query, function(i) {
    if(is.na(i)|i == "CL:0000000"){
      print('skip undefined or Cell')
      c(NA, NA)
    }else{
      print(i)
      new_string <- gsub(":", "_", i)
      result <- get_ontology_parents(new_string)
      new_res = gsub("_", ":", result)
      c(c(new_res[2],new_res[1]))
    }
  })
  out <- t(out)
  out <- as.data.frame(out)
  colnames(out) <- c('parent_CLname','parent_CLID')
  out = cbind(df, out)
  res = cbind(res,out)
}

colnames(res) = c("cluster", "genes", "manual_annotation", "manual_CLname", "manual_CLID",
                  "manual_parent_CLname", "manual_parent_CLID", "ERNIE", "ERNIE_CLname", "ERNIE_CLID",
                  "ERNIE_parent_CLname", "ERNIE_parent_CLID", "Gemini", "Gemini_CLname", "Gemini_CLID",
                  "Gemini_parent_CLname", "Gemini_parent_CLID", "ChatGPT", "ChatGPT_CLname", "ChatGPT_CLID",
                  "ChatGPT_parent_CLname", "ChatGPT_parent_CLID", "Llama", "Llama_CLname", "Llama_CLID",
                  "Llama_parent_CLname", "Llama_parent_CLID", "Claude", "Claude_CLname", "Claude_CLID",
                  "Claude_parent_CLname", "Claude_parent_CLID")
res[res == ""] <- NA

library(purrr)
res1 = res
for(i in 1:nrow(res)){
  print(i)
  ERNIE_CLname = res$ERNIE_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  ERNIE_CLname = res$ERNIE_CLname[i]
  ERNIE_parent_CLname = unlist(str_split(res$ERNIE_parent_CLname[i] ,';'))
  if(!is.na(ERNIE_CLname) & ERNIE_CLname != c('cell') & manual_CLname %in% ERNIE_CLname){
    res1$ERNIE_score[i] = 3
    res1$ERNIE_match[i] = manual_CLname
    print('ERNIE')
  }else if(any(!is.na(ERNIE_parent_CLname)) & any(ERNIE_parent_CLname != c('cell')) & manual_CLname %in% ERNIE_parent_CLname){
    res1$ERNIE_score[i] = 2
    res1$ERNIE_match[i] = manual_CLname
    print('ERNIE')
  }else if(!is.na(ERNIE_CLname) & ERNIE_CLname != c('cell') & any(manual_parent_CLname %in% ERNIE_CLname)){
    res1$ERNIE_score[i] = 1
    res1$ERNIE_match[i] = ERNIE_CLname
    print('ERNIE')
  }else if(any(!is.na(ERNIE_parent_CLname)) & any(manual_parent_CLname %in% ERNIE_parent_CLname)){
    res1$ERNIE_score[i] = 1
    res1$ERNIE_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% ERNIE_parent_CLname)]
    print('ERNIE')
  }else{
    res1$ERNIE_score[i] = 0
    res1$ERNIE_match[i] = 'no match'
    print('ERNIE')
  }
  Gemini_CLname = res$Gemini_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  Gemini_CLname = res$Gemini_CLname[i]
  Gemini_parent_CLname = unlist(str_split(res$Gemini_parent_CLname[i] ,';'))
  if(!is.na(Gemini_CLname) & Gemini_CLname != c('cell') & manual_CLname %in% Gemini_CLname){
    res1$Gemini_score[i] = 3
    res1$Gemini_match[i] = manual_CLname
    print('Gemini')
  }else if(any(!is.na(Gemini_parent_CLname)) & any(Gemini_parent_CLname != c('cell')) & manual_CLname %in% Gemini_parent_CLname){
    res1$Gemini_score[i] = 2
    res1$Gemini_match[i] = manual_CLname
    print('Gemini')
  }else if(!is.na(Gemini_CLname) & Gemini_CLname != c('cell') & any(manual_parent_CLname %in% Gemini_CLname)){
    res1$Gemini_score[i] = 1
    res1$Gemini_match[i] = Gemini_CLname
    print('Gemini')
  }else if(any(!is.na(Gemini_parent_CLname)) & any(manual_parent_CLname %in% Gemini_parent_CLname)){
    res1$Gemini_score[i] = 1
    res1$Gemini_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% Gemini_parent_CLname)]
    print('Gemini')
  }else{
    res1$Gemini_score[i] = 0
    res1$Gemini_match[i] = 'no match'
    print('Gemini')
  }
  ChatGPT_CLname = res$ChatGPT_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  ChatGPT_CLname = res$ChatGPT_CLname[i]
  ChatGPT_parent_CLname = unlist(str_split(res$ChatGPT_parent_CLname[i] ,';'))
  if(!is.na(ChatGPT_CLname) & ChatGPT_CLname != c('cell') & manual_CLname %in% ChatGPT_CLname){
    res1$ChatGPT_score[i] = 3
    res1$ChatGPT_match[i] = manual_CLname
    print('ChatGPT')
  }else if(any(!is.na(ChatGPT_parent_CLname)) & any(ChatGPT_parent_CLname != c('cell')) & manual_CLname %in% ChatGPT_parent_CLname){
    res1$ChatGPT_score[i] = 2
    res1$ChatGPT_match[i] = manual_CLname
    print('ChatGPT')
  }else if(!is.na(ChatGPT_CLname) & ChatGPT_CLname != c('cell') & any(manual_parent_CLname %in% ChatGPT_CLname)){
    res1$ChatGPT_score[i] = 1
    res1$ChatGPT_match[i] = ChatGPT_CLname
    print('ChatGPT')
  }else if(any(!is.na(ChatGPT_parent_CLname)) & any(manual_parent_CLname %in% ChatGPT_parent_CLname)){
    res1$ChatGPT_score[i] = 1
    res1$ChatGPT_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% ChatGPT_parent_CLname)]
    print('ChatGPT')
  }else{
    res1$ChatGPT_score[i] = 0
    res1$ChatGPT_match[i] = 'no match'
    print('ChatGPT')
  }
  Llama_CLname = res$Llama_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  Llama_CLname = res$Llama_CLname[i]
  Llama_parent_CLname = unlist(str_split(res$Llama_parent_CLname[i] ,';'))
  if(!is.na(Llama_CLname) & Llama_CLname != c('cell') & manual_CLname %in% Llama_CLname){
    res1$Llama_score[i] = 3
    res1$Llama_match[i] = manual_CLname
    print('Llama')
  }else if(any(!is.na(Llama_parent_CLname)) & any(Llama_parent_CLname != c('cell')) & manual_CLname %in% Llama_parent_CLname){
    res1$Llama_score[i] = 2
    res1$Llama_match[i] = manual_CLname
    print('Llama')
  }else if(!is.na(Llama_CLname) & Llama_CLname != c('cell') & any(manual_parent_CLname %in% Llama_CLname)){
    res1$Llama_score[i] = 1
    res1$Llama_match[i] = Llama_CLname
    print('Llama')
  }else if(any(!is.na(Llama_parent_CLname)) & any(manual_parent_CLname %in% Llama_parent_CLname)){
    res1$Llama_score[i] = 1
    res1$Llama_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% Llama_parent_CLname)]
    print('Llama')
  }else{
    res1$Llama_score[i] = 0
    res1$Llama_match[i] = 'no match'
    print('Llama')
  }
  Claude_CLname = res$Claude_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  Claude_CLname = res$Claude_CLname[i]
  Claude_parent_CLname = unlist(str_split(res$Claude_parent_CLname[i] ,';'))
  if(!is.na(Claude_CLname) & Claude_CLname != c('cell') & manual_CLname %in% Claude_CLname){
    res1$Claude_score[i] = 3
    res1$Claude_match[i] = manual_CLname
    print('Claude')
  }else if(any(!is.na(Claude_parent_CLname)) & any(Claude_parent_CLname != c('cell')) & manual_CLname %in% Claude_parent_CLname){
    res1$Claude_score[i] = 2
    res1$Claude_match[i] = manual_CLname
    print('Claude')
  }else if(!is.na(Claude_CLname) & Claude_CLname != c('cell') & any(manual_parent_CLname %in% Claude_CLname)){
    res1$Claude_score[i] = 1
    res1$Claude_match[i] = Claude_CLname
    print('Claude')
  }else if(any(!is.na(Claude_parent_CLname)) & any(manual_parent_CLname %in% Claude_parent_CLname)){
    res1$Claude_score[i] = 1
    res1$Claude_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% Claude_parent_CLname)]
    print('Claude')
  }else{
    res1$Claude_score[i] = 0
    res1$Claude_match[i] = 'no match'
    print('Claude')
  }
}

score_match = res1[,33:42]
res1$max_score <- do.call(pmax, c(res1[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score")], na.rm = TRUE))

# 定义一个函数来找到每一行最大分数对应的细胞类型
get_max_match <- function(row) {
  # 提取所有分数列
  scores <- c(row["ERNIE_score"], row["Gemini_score"], row["ChatGPT_score"], row["Llama_score"], row["Claude_score"])
  # 提取所有匹配列
  matches <- c(row["ERNIE_match"], row["Gemini_match"], row["ChatGPT_match"], row["Llama_match"], row["Claude_match"])
  # 找到最大分数的索引
  max_indices <- which(scores == max(scores, na.rm = TRUE))
  # 返回第一个最大分数对应的细胞类型
  matches[max_indices[1]]
}
# 应用函数到每一行
res1$max_match <- apply(res1, 1, get_max_match)
score = res1[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score",'max_score')]

write.csv(res1, './avg_log2FC_res_step2_v1.csv')

library(gridExtra)
library(ggsci)
v = pal_npg(palette = c("nrc"), alpha = 0.9)(9)
v[1] = 'orange'
coldf = data.frame(method = c('GTEx','Azimuth', 'BCL','coloncancer', 'TS', 'lungcancer','MCA','HCL', 'literature','cancer'),
                   color = c(v, v[6]),
                   stringsAsFactors = F)

# 将数据框转换为长格式，便于计算和绘图
score_long <- score %>%
  pivot_longer(cols = everything(), names_to = "Model", values_to = "Score")

# 计算每个模型得分的分布情况
score_distribution <- score_long %>%
  group_by(Model, Score) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Model) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 将得分分布补充为 0 分、1 分、2 分、3 分
score_distribution <- score_distribution %>%
  complete(Model, Score = 0:3, fill = list(Count = 0, Percentage = 0))

ggplot(score_distribution, aes(x = Model, y = Percentage, fill = factor(Score))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), position = position_stack(vjust = 0.5)) +
  scale_y_continuous(labels = percent_format(scale = 1), position = "right") +
  labs(x = "Model", y = "Percentage", fill = "Score", title = "gastric tumor Evaluation Scores") +
  theme_bw() +
  scale_fill_manual(values = coldf$color) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  coord_flip()

ggsave(
  filename = './avg_log2FC_res_v1_p.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 5,             # 宽
  height = 3,            # 高
  units = "in",          # 单位
  dpi = 100)

ggplot(score_distribution, aes(x = Model, y = Percentage, fill = factor(Score))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1),position = "right") +
  labs(x = "Model", y = "Percentage", fill = "Score", title = "gastric tumor")+
  theme_bw()+
  scale_fill_manual(values = coldf$color) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  coord_flip()
ggsave(
  filename = './avg_log2FC_res_v1.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 5,             # 宽
  height = 3,            # 高
  units = "in",          # 单位
  dpi = 100)

average_scores <- colMeans(score)
average_scores = as.data.frame(average_scores)
write.csv(average_scores,'./avg_log2FC_average_scores_v1.csv')

next_top_10_markers <- markers %>%
  group_by(cluster) %>%  # 按照cluster分组
  arrange(desc(avg_log2FC)) %>%  # 按照avg_log2FC降序排列
  slice(11:20)

cluster_genes <- next_top_10_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(genes = paste(gene, collapse=",")) %>%
  dplyr::ungroup()
cluster_genes

#get cell type

mg <- character()
cluster = numeric()
row_number <- 1  # 初始化行号

for (i in 1:nrow(res1)) {
  print(i)
  n = res1[i,]  # 获取当前行数据
  if (n$max_score != 3) {
    # 处理并确保以字符形式处理单元格内容，排除NA值
    ct <- unique(c(as.character(n$ERNIE_CLname), as.character(n$Gemini_CLname), as.character(n$ChatGPT_CLname), as.character(n$Claude_CLname), as.character(n$Llama_CLname)))
    ct <- na.omit(ct)  # 排除NA值
    if (length(ct) > 0) {  # 只有在ct不为空的情况下才进行下一步
      ct_text <- paste(ct, collapse=", ")
      mg <- c(mg, paste("row", row_number, ":", ct_text))
    }
    cluster = c(cluster,i)
  } else {
    print('perfect')
  }
  row_number <- row_number + 1  # 增加行号
}

# 将mg的各行合并成一个长字符串，每行之间用换行符分隔
final_text <- paste(mg, collapse="\n")

species = 'human'
# 生成最终输出
output_text <- paste('Provide key marker genes for the following ', species, ' cell types, with 15 key marker genes per cell type. Provide only the abbreviated gene names of key marker genes, full names are not required:\n', final_text, '\nThe format of the final response should be:\n\n1: gene1, gene2, gene3\n2: gene1, gene2, gene3\nN: gene1, gene2, gene3\n\n...where N represents the row number and gene1, gene2, gene3 represent key marker genes.', sep="")

print(output_text)

cluster

markerdata <- data_frame(
  Row = c(3, 5, 8),
  Marker_Gene = c(
    "NCR1, KIR2DL1, KIR2DL3, KIR2DL4, KIR3DL1, KIR3DL2, KIR3DL3, CD16, CD56, CD57, GZMB, PRF1, KLRD1, KLRC1, KLRF1",
    "CD3D, CD3E, CD3G, CD2, CD5, CD7, CD28, CD6, IL7R, TRAC, TRBC1, TRBC2, TRGC1, TRGC2, TRDC, CXCR5, PD1, BCL6, SH2D1A, IL21, SAP, ICOS, CD40LG, CD84, IL6ST, BATF, CD200, CD278, IL4, CXCL13, CCR7, CD45RA, SELL, TCF7, LEF1, CD27, CD62L, BCL11B, LCK, ZAP70",
    "ACTA2, MYH11, TAGLN, CNN1, CALD1, MYLK, TPM2, SM22, LMOD1, ACTG2, MYL9, SMTN, VIM, DES, TPM1"
  )
)

cluster = c(3, 5, 8)

unique(Idents(seurat_obj))
levels_order <- levels(Idents(seurat_obj))
pn_gene <- data.frame(cluters = numeric(), positive_marker = character(), negative_marker = character(), stringsAsFactors = FALSE)

for(i in cluster){
  print(i)
  all_cells <- subset(seurat_obj, idents = levels_order[i])
  marker_gene = unlist(str_split(markerdata[markerdata$Row == i,]$Marker_Gene, ', '))
  positive_marker = character()
  negative_marker = character()
  for(j in marker_gene){
    tryCatch({
      if(any(rownames(all_cells) == j)){
        print(j)
        code <- paste0("positive_cell = subset(all_cells,", j , ">0)")
        eval(parse(text = code))
        percent = nrow(positive_cell@meta.data)/nrow(all_cells@meta.data)
        if (percent>=0.8){
          positive_marker = c(positive_marker,j)
        }else{
          negative_marker = c(negative_marker,j)
        }
      }else{
        print(j)
        print('Can not find marker gene in the data')
      }
    }, error = function(e) {
      # 打印错误消息
      print(paste("Error at marker gene ", j, ":", e$message))
    })
  }
  positive_marker = paste(positive_marker, collapse=",")
  negative_marker = paste(negative_marker, collapse=",")
  pn_gene1 <- as.data.frame(t(c(i, positive_marker,negative_marker)))
  pn_gene = rbind(pn_gene,pn_gene1)
}

colnames(pn_gene) = c('row','positive_marker','negative_marker')
cluster_genes$row = rownames(cluster_genes)
cluster_genes = cluster_genes%>%
  right_join(pn_gene, by = c("row"))

#cluster_genes = m
cluster_genes$positive_marker = paste(cluster_genes$genes, cluster_genes$positive_marker, sep=",")
cluster_genes$positive_marker = gsub(',$','',cluster_genes$positive_marker)
cluster_genes = cluster_genes[,3:5]

# Function to generate the gene text
generate_text <- function(positive_gene, negative_gene) {
  # Create an empty list to store the text for each row
  gene_text <- list()
  negative_gene_text <- list()

  # Get all unique row numbers
  all_rows <- unique(as.numeric(cluster_genes$row))

  # Iterate through all possible row numbers
  for (i in 1:max(all_rows)) {
    # Get the row name
    row_name <- as.numeric(i)

    # Check if the row exists in the data frame
    if (row_name %in% cluster_genes$row) {
      # Get the positive genes for the current row
      positive_genes <- strsplit(cluster_genes$positive_marker[cluster_genes$row == row_name], ",")[[1]]

      # Get the negative genes for the current row
      negative_genes <- strsplit(cluster_genes$negative_marker[cluster_genes$row == row_name], ",")[[1]]

      # Create the text for the current row
      gene_text[[row_name]] <- paste0(
        "row", row_name, " = c(",
        paste0("'", positive_genes, "'", collapse = ","),
        "),"
      )

      # Create the text for the current row
      negative_gene_text[[row_name]] <- paste0(
        "row", row_name, " = c(",
        paste0("'", negative_genes, "'", collapse = ","),
        "),"
      )
    } else {
      # If row doesn't exist, create empty entries
      gene_text[[row_name]] <- paste0("row", row_name, " = c(),")
      negative_gene_text[[row_name]] <- paste0("row", row_name, " = c(),")
    }
  }

  # Create the final text string
  text <- paste0(
    "inter_res_1 = LLM_interect(positive_gene = list(\n",
    paste0(gene_text, collapse = "\n"),
    "),\n",
    "negative_gene = list(\n",
    paste0(negative_gene_text, collapse = "\n"),
    "))"
  )

  # Return the text
  return(text)
}



# Generate the text
text <- generate_text(cluster_genes$positive_marker, cluster_genes$negative_marker)

# Print the text
cat(text)



#########1
inter_res_1 = LLM_interect(positive_gene = list(
  row1 = c(),
  row2 = c(),
  row3 = c('GZMH','CTSW','CD247','KLRC2','CMC1','TYROBP','CLIC3','CCL4','CCL4L2','HOPX','GZMB','KLRD1'),
  row4 = c(),
  row5 = c('CD2','TRAC','PIK3IP1','FYB','TIGIT','CD3E','CD3D','RORA','CREM','PBXIP1','CD3D'),
  row6 = c(),
  row7 = c(),
  row8 = c('NOTCH3','CSRP2','IGFBP7','MCAM','TPM1','BGN','TINAGL1','COL18A1','MYLK','C11orf96','ACTA2','TAGLN','CALD1','TPM2','MYL9','VIM','TPM1')),
  negative_gene = list(
    row1 = c(),
    row2 = c(),
    row3 = c('NCR1','KIR2DL1','KIR2DL3','KIR2DL4','KIR3DL1','KIR3DL2','KIR3DL3','PRF1','KLRC1','KLRF1'),
    row4 = c(),
    row5 = c('CD4','IL2RA','IL7R','CCR7','SELL','CD40LG','IL4','IL5','IL13','IL17A','FOXP3','GATA3','TBX21','RORC','CTLA4','CD3E','CD3G','CD2','CD5','CD7','CD28','CD6','TRAC','TRBC1','TRBC2','TRGC1','TRGC2','CXCR5','BCL6','SH2D1A','IL21','ICOS','CD200'),
    row6 = c(),
    row7 = c(),
    row8 = c('MYH11','CNN1','MYLK','LMOD1','ACTG2','SMTN','DES')))

res = inter_res_1

cluster_genes = readRDS('./gc_cluster_genes_inter1.rds')
res1 = create_dataframe(res, markers)
res1$cluster = as.character(res1$cluster)
res1[11,1] = 'B cells'
res1[2,1] = 'Epithelial cells'
#write.csv(res1,'./avg_log2FC_step1_res1.csv')
res1$row = as.numeric(res1$clusters+1)
cluster_genes$row = as.numeric(cluster_genes$row)
res2 = res1%>%left_join(cluster_genes, by = c('row'))
res3 = read.csv('./avg_log2FC_step1_res1.csv')
res3 <- res3[ , -1]
res = cbind(res3, res2[,10:11],res2[,4:8])

write.csv(res,'./avg_log2FC_step1_res1_inter1.csv')
#res1 = read.csv('./avg_log2FC_step1_res1.csv')
#res1 <- res1[ , -1]

cl = res1[,1:2]
OlsSearch(q = 'axial mesoderm', ontology = "cl")
for(j in c('cluster','ERNIE','Gemini','ChatGPT','Llama','Claude')){
  query <- res1[[j]]
  cv <- sapply(query, function(i) {
    if(i == "undefined"|i == "undefined "){
      print('skip undefined')
      c(i, NA, NA)
    }else{
      print(i)
      new_string <- tolower(as.character(i))
      new_string <- sub("s$", "", new_string)
      d <- py$get_class_id_from_ontology("D:/R/LLMsCellType/cl.owl",new_string)
      if(d == c("no match")){
        d <- as(olsSearch(OlsSearch(q = new_string, ontology = "cl")), "data.frame")
        d <- d[grep('CL:', d$obo_id),]
        c(i, d[1, 'label'], d[1, 'obo_id'])
      }else{
        d = gsub("_", ":", d)
        c(i, new_string, d)
      }
    }
  })
  cv <- t(cv)
  cv <- as.data.frame(cv)
  colnames(cv) <- c(j,'CLname','CLID')
  cl = cbind(cl,cv)
}

res = cl[,1:2]
for(j in c(5,8,11,14,17,20)){
  df = cl[,c(j-2):j]
  print(colnames(df)[[1]])
  query = df[[3]]
  out <- sapply(query, function(i) {
    if(is.na(i)|i == "CL:0000000"){
      print('skip undefined or Cell')
      c(NA, NA)
    }else{
      print(i)
      new_string <- gsub(":", "_", i)
      result <- get_ontology_parents(new_string)
      new_res = gsub("_", ":", result)
      c(c(new_res[2],new_res[1]))
    }
  })
  out <- t(out)
  out <- as.data.frame(out)
  colnames(out) <- c('parent_CLname','parent_CLID')
  out = cbind(df, out)
  res = cbind(res,out)
}

colnames(res) = c("cluster", "genes", "manual_annotation", "manual_CLname", "manual_CLID",
                  "manual_parent_CLname", "manual_parent_CLID", "ERNIE", "ERNIE_CLname", "ERNIE_CLID",
                  "ERNIE_parent_CLname", "ERNIE_parent_CLID", "Gemini", "Gemini_CLname", "Gemini_CLID",
                  "Gemini_parent_CLname", "Gemini_parent_CLID", "ChatGPT", "ChatGPT_CLname", "ChatGPT_CLID",
                  "ChatGPT_parent_CLname", "ChatGPT_parent_CLID", "Llama", "Llama_CLname", "Llama_CLID",
                  "Llama_parent_CLname", "Llama_parent_CLID", "Claude", "Claude_CLname", "Claude_CLID",
                  "Claude_parent_CLname", "Claude_parent_CLID")
res[res == ""] <- NA

library(purrr)
res1 = res
for(i in 1:nrow(res)){
  print(i)
  ERNIE_CLname = res$ERNIE_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  ERNIE_CLname = res$ERNIE_CLname[i]
  ERNIE_parent_CLname = unlist(str_split(res$ERNIE_parent_CLname[i] ,';'))
  if(!is.na(ERNIE_CLname) & ERNIE_CLname != c('cell') & manual_CLname %in% ERNIE_CLname){
    res1$ERNIE_score[i] = 3
    res1$ERNIE_match[i] = manual_CLname
    print('ERNIE')
  }else if(any(!is.na(ERNIE_parent_CLname)) & any(ERNIE_parent_CLname != c('cell')) & manual_CLname %in% ERNIE_parent_CLname){
    res1$ERNIE_score[i] = 2
    res1$ERNIE_match[i] = manual_CLname
    print('ERNIE')
  }else if(!is.na(ERNIE_CLname) & ERNIE_CLname != c('cell') & any(manual_parent_CLname %in% ERNIE_CLname)){
    res1$ERNIE_score[i] = 1
    res1$ERNIE_match[i] = ERNIE_CLname
    print('ERNIE')
  }else if(any(!is.na(ERNIE_parent_CLname)) & any(manual_parent_CLname %in% ERNIE_parent_CLname)){
    res1$ERNIE_score[i] = 1
    res1$ERNIE_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% ERNIE_parent_CLname)]
    print('ERNIE')
  }else{
    res1$ERNIE_score[i] = 0
    res1$ERNIE_match[i] = 'no match'
    print('ERNIE')
  }
  Gemini_CLname = res$Gemini_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  Gemini_CLname = res$Gemini_CLname[i]
  Gemini_parent_CLname = unlist(str_split(res$Gemini_parent_CLname[i] ,';'))
  if(!is.na(Gemini_CLname) & Gemini_CLname != c('cell') & manual_CLname %in% Gemini_CLname){
    res1$Gemini_score[i] = 3
    res1$Gemini_match[i] = manual_CLname
    print('Gemini')
  }else if(any(!is.na(Gemini_parent_CLname)) & any(Gemini_parent_CLname != c('cell')) & manual_CLname %in% Gemini_parent_CLname){
    res1$Gemini_score[i] = 2
    res1$Gemini_match[i] = manual_CLname
    print('Gemini')
  }else if(!is.na(Gemini_CLname) & Gemini_CLname != c('cell') & any(manual_parent_CLname %in% Gemini_CLname)){
    res1$Gemini_score[i] = 1
    res1$Gemini_match[i] = Gemini_CLname
    print('Gemini')
  }else if(any(!is.na(Gemini_parent_CLname)) & any(manual_parent_CLname %in% Gemini_parent_CLname)){
    res1$Gemini_score[i] = 1
    res1$Gemini_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% Gemini_parent_CLname)]
    print('Gemini')
  }else{
    res1$Gemini_score[i] = 0
    res1$Gemini_match[i] = 'no match'
    print('Gemini')
  }
  ChatGPT_CLname = res$ChatGPT_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  ChatGPT_CLname = res$ChatGPT_CLname[i]
  ChatGPT_parent_CLname = unlist(str_split(res$ChatGPT_parent_CLname[i] ,';'))
  if(!is.na(ChatGPT_CLname) & ChatGPT_CLname != c('cell') & manual_CLname %in% ChatGPT_CLname){
    res1$ChatGPT_score[i] = 3
    res1$ChatGPT_match[i] = manual_CLname
    print('ChatGPT')
  }else if(any(!is.na(ChatGPT_parent_CLname)) & any(ChatGPT_parent_CLname != c('cell')) & manual_CLname %in% ChatGPT_parent_CLname){
    res1$ChatGPT_score[i] = 2
    res1$ChatGPT_match[i] = manual_CLname
    print('ChatGPT')
  }else if(!is.na(ChatGPT_CLname) & ChatGPT_CLname != c('cell') & any(manual_parent_CLname %in% ChatGPT_CLname)){
    res1$ChatGPT_score[i] = 1
    res1$ChatGPT_match[i] = ChatGPT_CLname
    print('ChatGPT')
  }else if(any(!is.na(ChatGPT_parent_CLname)) & any(manual_parent_CLname %in% ChatGPT_parent_CLname)){
    res1$ChatGPT_score[i] = 1
    res1$ChatGPT_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% ChatGPT_parent_CLname)]
    print('ChatGPT')
  }else{
    res1$ChatGPT_score[i] = 0
    res1$ChatGPT_match[i] = 'no match'
    print('ChatGPT')
  }
  Llama_CLname = res$Llama_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  Llama_CLname = res$Llama_CLname[i]
  Llama_parent_CLname = unlist(str_split(res$Llama_parent_CLname[i] ,';'))
  if(!is.na(Llama_CLname) & Llama_CLname != c('cell') & manual_CLname %in% Llama_CLname){
    res1$Llama_score[i] = 3
    res1$Llama_match[i] = manual_CLname
    print('Llama')
  }else if(any(!is.na(Llama_parent_CLname)) & any(Llama_parent_CLname != c('cell')) & manual_CLname %in% Llama_parent_CLname){
    res1$Llama_score[i] = 2
    res1$Llama_match[i] = manual_CLname
    print('Llama')
  }else if(!is.na(Llama_CLname) & Llama_CLname != c('cell') & any(manual_parent_CLname %in% Llama_CLname)){
    res1$Llama_score[i] = 1
    res1$Llama_match[i] = Llama_CLname
    print('Llama')
  }else if(any(!is.na(Llama_parent_CLname)) & any(manual_parent_CLname %in% Llama_parent_CLname)){
    res1$Llama_score[i] = 1
    res1$Llama_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% Llama_parent_CLname)]
    print('Llama')
  }else{
    res1$Llama_score[i] = 0
    res1$Llama_match[i] = 'no match'
    print('Llama')
  }
  Claude_CLname = res$Claude_CLname[i]
  manual_CLname = res$manual_CLname[i]
  manual_parent_CLname = unlist(str_split(res$manual_parent_CLname[i],';'))
  Claude_CLname = res$Claude_CLname[i]
  Claude_parent_CLname = unlist(str_split(res$Claude_parent_CLname[i] ,';'))
  if(!is.na(Claude_CLname) & Claude_CLname != c('cell') & manual_CLname %in% Claude_CLname){
    res1$Claude_score[i] = 3
    res1$Claude_match[i] = manual_CLname
    print('Claude')
  }else if(any(!is.na(Claude_parent_CLname)) & any(Claude_parent_CLname != c('cell')) & manual_CLname %in% Claude_parent_CLname){
    res1$Claude_score[i] = 2
    res1$Claude_match[i] = manual_CLname
    print('Claude')
  }else if(!is.na(Claude_CLname) & Claude_CLname != c('cell') & any(manual_parent_CLname %in% Claude_CLname)){
    res1$Claude_score[i] = 1
    res1$Claude_match[i] = Claude_CLname
    print('Claude')
  }else if(any(!is.na(Claude_parent_CLname)) & any(manual_parent_CLname %in% Claude_parent_CLname)){
    res1$Claude_score[i] = 1
    res1$Claude_match[i] = manual_parent_CLname[c(manual_parent_CLname %in% Claude_parent_CLname)]
    print('Claude')
  }else{
    res1$Claude_score[i] = 0
    res1$Claude_match[i] = 'no match'
    print('Claude')
  }
}

score_match = res1[,33:42]
res1$max_score <- do.call(pmax, c(res1[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score")], na.rm = TRUE))

# 应用函数到每一行
res1$max_match <- apply(res1, 1, get_max_match)
score = res1[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score",'max_score')]

write.csv(res1, './avg_log2FC_res_step2_v1_inter1.csv')

# 将数据框转换为长格式，便于计算和绘图
score_long <- score %>%
  pivot_longer(cols = everything(), names_to = "Model", values_to = "Score")

# 计算每个模型得分的分布情况
score_distribution <- score_long %>%
  group_by(Model, Score) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Model) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 将得分分布补充为 0 分、1 分、2 分、3 分
score_distribution <- score_distribution %>%
  complete(Model, Score = 0:3, fill = list(Count = 0, Percentage = 0))

ggplot(score_distribution, aes(x = Model, y = Percentage, fill = factor(Score))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), position = position_stack(vjust = 0.5)) +
  scale_y_continuous(labels = percent_format(scale = 1), position = "right") +
  labs(x = "Model", y = "Percentage", fill = "Score", title = "gastric tumor Evaluation Scores") +
  theme_bw() +
  scale_fill_manual(values = coldf$color) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  coord_flip()

ggsave(
  filename = './avg_log2FC_res_v1_inter1_p.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 5,             # 宽
  height = 3,            # 高
  units = "in",          # 单位
  dpi = 100)

ggplot(score_distribution, aes(x = Model, y = Percentage, fill = factor(Score))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1),position = "right") +
  labs(x = "Model", y = "Percentage", fill = "Score", title = "gastric tumor")+
  theme_bw()+
  scale_fill_manual(values = coldf$color) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  coord_flip()
ggsave(
  filename = './avg_log2FC_res_inter1_v1.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 5,             # 宽
  height = 3,            # 高
  units = "in",          # 单位
  dpi = 100)

average_scores <- colMeans(score)
average_scores = as.data.frame(average_scores)
write.csv(average_scores,'./avg_log2FC_average_scores_inter_v1.csv')

res1 = read.csv('./avg_log2FC_res_step2_v1_inter1.csv')

#get cell type

MarkerPrompt = function(data,model,species){
  mg <- character()
  cluster = numeric()
  row_number <- 1  # 初始化行号
  res1 = data
  for (i in 1:nrow(res1)) {
    n = res1[i,]# 获取当前行数据
    # 处理并确保以字符形式处理单元格内容，排除NA值
    ct <- as.character(n[[model]])
    ct <- na.omit(ct)  # 排除NA值
    if (length(ct) > 0) {  # 只有在ct不为空的情况下才进行下一步
      ct_text <- paste(ct, collapse=", ")
      mg <- c(mg, paste("row", row_number, ":", ct_text))
    }
    cluster = c(cluster,i)
    row_number <- row_number + 1  # 增加行号
  }

  # 将mg的各行合并成一个长字符串，每行之间用换行符分隔
  final_text <- paste(mg, collapse="\n")

  specie = species
  # 生成最终输出
  output_text <- paste('Provide key marker genes for the following ', specie, ' cell types, with 15 key marker genes per cell type. Provide only the abbreviated gene names of key marker genes, full names are not required:\n', final_text, '\nThe format of the final response should be:\n\n1: gene1, gene2, gene3\n2: gene1, gene2, gene3\nN: gene1, gene2, gene3\n\n...where N represents the row number and gene1, gene2, gene3 represent key marker genes.', sep="")
  print(output_text)
}

res1$ERNIE_CLname
res1$Claude_CLname
MarkerPrompt(res1, 'manual_CLname', species = "human")
MarkerPrompt(res1, 'ERNIE_CLname', species = "human")
MarkerPrompt(res1, 'ChatGPT_CLname', species = "human")
MarkerPrompt(res1, 'Gemini_CLname', species = "human")
MarkerPrompt(res1, 'Llama_CLname', species = "human")
MarkerPrompt(res1, 'Claude_CLname', species = "human")
