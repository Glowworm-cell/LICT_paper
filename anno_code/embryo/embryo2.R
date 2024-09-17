library(LLMsCellIdentifier)
library(reticulate)
library(rols)
library(dplyr)
library(scales)
library(Seurat)
path <- Sys.which("python")
Sys.setenv(RETICULATE_PYTHON = path)
Sys.setenv(Llama3_api_key = 'xxx')
Sys.setenv(Llama3_secret_key = 'xxx')
Sys.setenv(ERNIE_api_key = 'xxx')
Sys.setenv(ERNIE_secret_key = 'xxx')
Sys.setenv(GEMINI_api_key = 'xxx')
Sys.setenv(openai.api_key = 'xxx')
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
ANTHROPIC_API_KEY = 'xxx'
")
seurat_obj = readRDS('../human.rds')
seurat_obj = FindClusters(seurat_obj, res = 2)
Idents(seurat_obj) = factor(seurat_obj@meta.data$cluster_id)
markers <- FindAllMarkers(object = seurat_obj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
LLMCelltyperes = LLMCelltype(FindAllMarkersResult = markers,
                             species = 'human',
                             topgenenumber = 10,
                             tissuename = 'embryo')

res = LLMCelltyperes
res[['Gemini']] = Geminires


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


#res1 = read.csv('./embryo_res/embryo_res1.csv')
#res1 <- res1[ , -1]
res1$cluster = as.character(res1$cluster)
res1[11,1] = 'Extraembryonic mesenchymal cells'
res1[1,1] = 'Advanced Mesoderm cells'
res1[2,1] = 'axial mesoderm cells'
res1[3,1] = 'surface ectodermal cell'
res1[4,1] = 'Emergent Mesoderm cells'
res1[5,1] = 'Endoderm cells'
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
  Row = c(1, 2, 3, 5, 8, 9, 10),
  Marker_Gene = c(
    "GATA4, GATA6, SOX7, HAND1, HAND2, NKX2-5, MYH6, MYH7, CDX2, TEAD4, CK, HLA-G, PDGFRB, KDR, ESRRB",
    "FOXJ1, DNAI1, DNAH5, DNAH9, TUBB4B, CCDC39, CCDC40, SPEF2, DYX1C1, RSPH4A, RSPH9, CFAP298, MUC1, MUC4, MUC16",
    "SOX10, SNAI2, FOXD3, NGFR, NCAM1, FCGR3A, KIR2DL1, KIR3DL1, VWF, PECAM1, KDR, PECAM1, CDH5, ITGA2, ITGB1",
    "ALB, AAT, TTR, AFP, HNF4A, HNF1A, CYP3A4, CYP2E1, GATA4, SOX17, FOXA2, FOXA3, HHEX, FABP1, CK",
    "PPBP, PF4, ITGA2B, ITGB3, GP1BA, GP9, VWF, MPL, GATA1, THPO, ENG, CD42b, CD61, TBXAS1, FGA",
    "VIM, WT1, PAX2, LHX1, TBXT, SHH, NKX6-1, MAP2, NEFL, SYP, GAP43, RBFOX3, NGN2, ISL1, MUSASHI",
    "CDX2, TEAD4, CK, CGA, CGB5, HLA-G, EPCAM, CDH1, CK, CK, CK, TP63, IVL, KLF4, P63"
  )
)

cluster = c(1, 2, 3, 5, 8, 9, 10)

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
  row1 = c('PMP22','RGS2','COL3A1','CYP27A1','EMP2','BMPER','VIM','SNAI2','KRT8','SLC9A3R1','GATA6','HAND1'),
  row2 = c('SPAG6','DNAH9','DRC1','DNAH5','DNALI1','C9orf24','TEKT1','CETN2','NRAV','OSCP1','DNAH5','DNAH9','RSPH4A','RSPH9','CFAP298'),
  row3 = c('MSX2','HERPUD1','TPM1','GLB1','ABCG2','KRT18','KRT7','TMEM54','DSP','DLX5','ITGB1'),
  row4 = c(),
  row5 = c('RBP4','APOA4','AGT','SERPINA1','FLRT3','CCKBR','IGFBP6','FGA','GPC3','APOE'),
  row6 = c(),
  row7 = c(),
  row8 = c('THBS1','LCP1','PECAM1','AIF1','HLA.E','GPX1','B2M','MEF2C','S100A4','CCL3'),
  row9 = c('KCNK17','LINC01356','RSPO3','ABLIM1','TBXT','TMSB15A','CKS1B','PPIAP46','HMGB3','HVCN1'),
  row10 = c('SCG3','POU5F1','TBXT','L1TD1','FABP5','EFEMP1','ATP5PD','TSTD1','UCHL1','PTPRZ1','EPCAM')),
  negative_gene = list(
    row1 = c('GATA4','SOX7','HAND2','CDX2','TEAD4','PDGFRB','KDR'),
    row2 = c('FOXJ1','DNAI1','TUBB4B','CCDC39','CCDC40','SPEF2','MUC4'),
    row3 = c('SNAI2','NGFR','PECAM1','KDR','PECAM1','ITGA2'),
    row4 = c(),
    row5 = c('ALB','TTR','AFP','HNF4A','HNF1A','CYP2E1','GATA4','SOX17','FOXA2','FOXA3','HHEX','FABP1'),
    row6 = c(),
    row7 = c(),
    row8 = c('PPBP','PF4','ITGA2B','ITGB3','GP1BA','GP9','VWF','MPL','GATA1','THPO','ENG','TBXAS1','FGA'),
    row9 = c('VIM','LHX1','TBXT','MAP2','NEFL','SYP','GAP43','RBFOX3','ISL1'),
    row10 = c('CDX2','TEAD4','CGA','CGB5','CDH1','KLF4')))


res = inter_res_1
res[['Gemini']] = Geminires


#cluster_genes = readRDS('./gc_cluster_genes_inter1.rds')
res1 = create_dataframe(res, markers)
res1$cluster = as.character(res1$cluster)
res1[11,1] = 'Extraembryonic mesenchymal cells'
res1[1,1] = 'Advanced Mesoderm cells'
res1[2,1] = 'axial mesoderm cells'
res1[3,1] = 'surface ectodermal cell'
res1[4,1] = 'Emergent Mesoderm cells'
res1[5,1] = 'Endoderm cells'
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

res1$manual_annotation
MarkerPrompt(res1, 'manual_annotation', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "MIXL1, EOMES, TBXT, GSC, FGF8, PDGFRA, NKX2-5, MESP1, WNT3A, TBX6, GATA4, ISL1, BMP4, EVX1, FOXF1",
    "TBXT, TBXT, SHH, FOXA2, CDX2, HOXB1, NOTCH1, LHX1, FGF8, WNT3A, GDF11, SOX2, EOMES, PITX2, BMP4",
    "TP63, KRT5, KRT14, KRT18, GRHL3, OCLN, CDH1, IVL, FLG, JUP, DSP, PKP1, KRT15, DSG1, PERP",
    "MIXL1, TBXT, EOMES, PDGFRA, KDR, CDX2, GSC, FGF8, TBX6, MESP1, FOXF1, NKX2-5, WNT3A, BMP4, GATA4",
    "SOX17, FOXA2, GATA4, GATA6, CXCR4, EPCAM, HNF4A, KRT19, AFP, CER1, PDX1, CDH1, SERPINA1, SOX7, HHEX",
    "NANOG, POU5F1, SOX2, KLF4, ZFP42, TDGF1, LIN28A, DPPA4, GDF3, SALL4, DPPA5, UTF1, TFCP2L1, LEFTY1, NODAL",
    "GYPA, HBA1, HBB, KLF1, EPOR, GATA1, ALAS2, HEMGN, BCL11A, LMO2, ANK1, SPTA1, AHSP, CA1, FECH",
    "KDR, CD34, CDH5, TAL1, GFI1B, RUNX1, SPI1, LMO2, TEK, GATA2, CD144, FLT3, PROM1, HAND1, ENG",
    "MIXL1, TBXT, EOMES, MESP1, PDGFRA, WNT3A, GSC, FGF8, BMP4, CDX2, TBX6, FOXF1, GATA4, EVX1, NKX2-5",
    "TBXT, MIXL1, EOMES, GSC, FGF8, WNT3A, TDGF1, MESP1, EVX1, CER1, NODAL, LEFTY1, TDGF1, CDX2, BMP4",
    "HAND1, GATA2, GATA3, CDX2, EOMES, ISL1, PDGFRA, TCF15, WNT5A, FOXF1, SOX7, BMP4, MSX2, TBX18, FGF10"
  )
)


cluster = markerdata$Row


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
write.csv(pn_gene,'./embyro_5_manual_pn_gene.csv')

MarkerPrompt(res1, 'ERNIE', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "APOA1, APOE, TTR, SOX17, AFP, FOXA2, GATA4, GATA6, HNF4A, CUBN, LRP2, DAB2, PDGFRA, KRT19, HNF1B",
    "FOXJ1, DNAI1, DNAH5, RSPH4A, CCDC39, CCDC40, CFAP300, SPAG6, HYDIN, MCIDAS, DNAH11, CCNO, LRRC6, TEKT1, RSPH9",
    "SOX10, PAX3, SNAI2, FOXD3, TFAP2A, MITF, NGFR, SOX9, EDN3, KIT, PHOX2B, WNT1, CRABP1, ZIC1, MSX1",
    "ALB, AFP, HNF4A, TTR, AAT, CYP3A4, KRT18, KRT19, GPC3, ASGR1, HAMP, CEBPA, TDO2, SULT2A1, FABP1",
    "NANOG, POU5F1, SOX2, LIN28A, KLF4, SALL4, DPPA4, ZFP42, TDGF1, UTF1, DNMT3B, TERT, GDF3, ESRRB, SSEA4",
    "GYPA, HBA1, HBB, KLF1, EPOR, GATA1, ALAS2, BCL11A, LMO2, CD36, AHSP, ANK1, SPTA1, HEMGN, CA1",
    "PECAM1, VWF, CD34, KDR, TEK, ENG, FLT1, NOS3, CADHERIN5, ROBO4, CXCR4, FLT4, NOTCH1, CDH5, VEPTP",
    "TBXT, MIXL1, MESP1, EOMES, GSC, WNT3A, PDGFRA, KDR, CDX2, EVX1, TBX6, BMP4, GATA4, FGF8, FOXF1",
    "POU5F1, SOX2, NANOG, ZFP42, OCT4, KLF4, LIN28A, FOXD3, NES, TDGF1, DPPA4, LEFTY1, NODAL, CRIPTO, FGF5",
    "CD44, CD73, CD90, CD105, CD166, CD29, CD71, SCA1, CD106, CD146, PDGFRB, SPP1, LIFR, NESTIN, OCT4"
  )
)


cluster = markerdata$Row


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
write.csv(pn_gene,'./embyro_5_ERNIE_pn_gene.csv')



MarkerPrompt(res1, 'ChatGPT', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "GATA4, GATA6, SOX17, PDGFRA, FOXA2, SOX7, LAMA1, LAMC1, DAB2, CXCR4, FGF5, PRL3D1, HNF1B, KRT8, KRT18",
    "FOXJ1, DNAI1, DNAH5, RSPH4A, CCDC39, CCDC40, CFAP300, SPAG6, HYDIN, MCIDAS, DNAH11, CCNO, LRRC6, TEKT1, RSPH9",
    "EPCAM, CDH1, OCLN, CD24, KRT8, KRT18, KRT19, MUC1, TJP1, CLDN1, CLDN7, ITGA6, ITGB4, TP63, DSP",
    "TBXT, MIXL1, MESP1, EOMES, GSC, WNT3A, PDGFRA, KDR, CDX2, EVX1, TBX6, BMP4, GATA4, FGF8, FOXF1",
    "GATA4, GATA6, SOX17, PDGFRA, FOXA2, SOX7, LAMA1, LAMC1, DAB2, CXCR4, FGF5, PRL3D1, HNF1B, KRT8, KRT18",
    "NANOG, POU5F1, SOX2, LIN28A, KLF4, SALL4, DPPA4, ZFP42, TDGF1, UTF1, DNMT3B, TERT, GDF3, ESRRB, SSEA4",
    "HBA1, HBB, GYPA, EPB42, ANK1, SPTA1, SPTB, HBG1, HBG2, HBE1, HBZ, HBD, AHSP, CA1, ALAS2",
    "CD45, CD3, CD19, CD56, CD4, CD8, CD14, CD16, CD20, CD57, HLA-DR, CD11b, CD11c, CD33, CD34",
    "MESP2, MSGN1, TBX6, RIPPLY2, LHX1, PAX3, CDX1, CDX2, WNT3A, FGF8, TBXT, HES7, DUSP6, EVX1, MYOD1",
    "TBXT, SHH, FOXA2, NOG, BRACHYURY, CDX1, CDX2, SOX5, SOX6, LAMC1, NKX3-2, ACAN, COL2A1, KRT19, HNF3B",
    "CD44, CD73, CD90, CD105, CD166, CD29, CD71, SCA1, CD106, CD146, PDGFRB, SPP1, LIFR, NESTIN, FAP"
  )
)


cluster = markerdata$Row


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
write.csv(pn_gene,'./embyro_5_ChatGPT_pn_gene.csv')


MarkerPrompt(res1, 'Gemini', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "NKX2-5, TBX5, GATA4, MEF2C, ISL1, HAND1, HAND2, MYH6, TNNT2, TNNI3, ACTC1, ACTN2, MYL2, MYL3, GJA1",
    "FOXJ1, DNAI1, DNAH5, RSPH4A, CCDC39, CCDC40, CFAP300, SPAG6, HYDIN, MCIDAS, DNAH11, CCNO, LRRC6, TEKT1, RSPH9",
    "SOX10, PAX3, SNAI2, FOXD3, TFAP2A, MITF, NGFR, SOX9, EDN3, KIT, PHOX2B, WNT1, CRABP1, ZIC1, MSX1",
    "T, TBXT, SHH, FOXA2, CDX2, HOXB1, NOTCH1, LHX1, FGF8, WNT3A, GDF11, SOX2, EOMES, PITX2, BMP4",
    "LAMININ, COL4A1, COL4A2, ENTACTIN, TRANSFERRIN, APOH, FGG, FGB, FGA, PLG, C3, C8A, A2M, SERPINC1, ALB",
    "NANOG, POU5F1, SOX2, LIN28A, KLF4, SALL4, DPPA4, ZFP42, TDGF1, UTF1, DNMT3B, TERT, GDF3, ESRRB, SSEA4",
    "GYPA, HBA1, HBB, KLF1, EPOR, GATA1, ALAS2, BCL11A, LMO2, CD36, AHSP, ANK1, SPTA1, HEMGN, CA1",
    "PECAM1, VWF, CD34, KDR, TEK, ENG, FLT1, NOS3, CADHERIN5, ROBO4, CXCR4, FLT4, NOTCH1, CDH5, VEPTP",
    "MESP2, MSGN1, TBX6, RIPPLY2, LHX1, PAX3, CDX1, CDX2, WNT3A, FGF8, TBXT, HES7, DUSP6, EVX1, MYOD1",
    "KRT7, KRT18, GATA6, GATA4, BMP4, FOXO1, PAX8, FGF2, AMOT, OTX2, ISL1, HAND1, VIM, CD44, PDGFRA",
    "T, MIXL1, GATA2, EOMES, HAND1, PDGFRA, WNT3, GSC, CER1, LEFTY1, NODAL, BMP4, FGF8, EVX1, MESP1"
  )
)


cluster = markerdata$Row


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
write.csv(pn_gene,'./embyro_5_Gemini_pn_gene.csv')


MarkerPrompt(res1, 'Llama', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "GATA3, KRT7, TEAD4, CDX2, HAND1, TFAP2C, HLA-G, ITGA6, ITGB1, ELF5, CSH1, GCM1, ERVW-1, EGFR, TP63",
    "FOXJ1, DNAI1, DNAH5, RSPH4A, CCDC39, CCDC40, CFAP300, SPAG6, HYDIN, MCIDAS, DNAH11, CCNO, LRRC6, TEKT1, RSPH9",
    "PECAM1, VWF, CD34, KDR, TEK, ENG, FLT1, NOS3, CADHERIN5, ROBO4, CXCR4, FLT4, NOTCH1, CDH5, VEPTP",
    "T, MIXL1, MESP1, EOMES, GSC, WNT3A, PDGFRA, KDR, CDX2, EVX1, TBX6, BMP4, GATA4, FGF8, FOXF1",
    "ALB, AFP, HNF4A, TTR, AAT, CYP3A4, KRT18, KRT19, GPC3, ASGR1, HAMP, CEBPA, TDO2, SULT2A1, FABP1",
    "NANOG, POU5F1, SOX2, LIN28A, KLF4, SALL4, DPPA4, ZFP42, TDGF1, UTF1, DNMT3B, TERT, GDF3, ESRRB, SSEA4",
    "GYPA, HBA1, HBB, KLF1, EPOR, GATA1, ALAS2, BCL11A, LMO2, CD36, AHSP, ANK1, SPTA1, HEMGN, CA1",
    "ITGA2B, ITGB3, GP1BA, GP1BB, GP9, PF4, PPBP, VWF, CD41, CD61, CD42a, CD42b, THBS1, F13A1, MPL",
    "MAP2, NEUN, TUBB3, SYN1, GAP43, SYP, PSD95, NCAM1, VGLUT1, GABA, NMDAR1, AMPAR1, CHAT, GFAP, NFL",
    "NANOG, POU5F1, SOX2, LIN28A, KLF4, SALL4, DPPA4, ZFP42, TDGF1, UTF1, DNMT3B, TERT, GDF3, ESRRB, SSEA4",
    "CD44, CD73, CD90, CD105, CD166, CD29, CD71, SCA1, CD106, CD146, PDGFRB, SPP1, LIFR, NESTIN, FAP"
  )
)


cluster = markerdata$Row


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
write.csv(pn_gene,'./embyro_5_Llama_pn_gene.csv')


MarkerPrompt(res1, 'Claude', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "GATA2, EOMES, T, MIXL1, HAND1, PDGFRA, ENG, KDR, CDX2, NODAL, LEFTY2, CER1, WNT3, FGF8, BMP4",
    "FOXJ1, DNAI1, DNAH9, DNAH5, CCDC39, DNAI2, RSPH4A, RSPH1, DYNLRB2, HYDIN, CCNO, MCIDAS, FOXI1, SPEF2, LRRC6",
    "CGB, CGA, KRT7, GCM1, ERVW-1, TEAD4, TFAP2C, CDX2, HLA-G, TP63, ITGA6, ITGB1, SDC1, MMP9, HAND1",
    "MIXL1, T, EOMES, GSC, FGF8, PDGFRA, MESP1, SNAIL, WNT3A, LEFTY1, TDGF1, CAD, EVX1, FOXH1, LHX1",
    "SOX17, FOXA2, GATA4, GATA6, HHEX, APOA1, AFP, CER1, LAMB1, PDGFRA, CXCR4, KRT19, FOXQ1, ALB, FABP1",
    "SOX2, NANOG, POU5F1, TDGF1, ZFP42, DPPA4, DPPA2, SALL4, KLF4, FGF2, T, FGF5, UTF1, LIN28, GBX2",
    "GYPA, HBB, HBA, ANK1, EPB42, SPTA1, SPTB, CD47, RHAG, KLF1, HEMGN, GATA1, UROS, ALAS2, CA1",
    "CD68, CD163, CD64, ITGAM, CSF1R, CX3CR1, CD14, MRC1, CCR2, FCER1G, MARCO, TNF, IL10, SPI1, TLR4",
    "NOG, SHH, FOXA2, T, FGF8, ACTA1, ACAN, COL2A1, BRACHYURY, NKX3-2, FOXA1, CHRD, LMX1A, PAX1, TBX18",
    "PAX6, SOX1, SOX2, NEUROG1, NEUROG2, NES, ASCL1, FOXG1, SIX3, HESX1, OTX2, LHX2, RGMA, NCAM1, TUBB3",
    "VIM, CD44, FN1, SNAI1, TWIST1, ACTA2, MMP2, ENG, PDGFRB, THY1, FAP, ITGA11, ZEB1, COL1A1, CD105")
)


cluster = markerdata$Row


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
write.csv(pn_gene,'./embyro_5_Claude_pn_gene.csv')
