setwd('D:/R/LLMsCellType/240617embryo_different_method6/')
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
seurat_obj = readRDS('D:/R/wb_2/human/human.rds')
seurat_obj = FindClusters(seurat_obj, res = 2)
Idents(seurat_obj) = factor(seurat_obj@meta.data$cluster_id)
markers <- FindAllMarkers(object = seurat_obj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
LLMCelltyperes = LLMCelltype(FindAllMarkersResult = markers,
                             species = 'human',
                             topgenenumber = 10,
                             tissuename = 'embryo')

Geminires = LLMCelltyperes$Gemini
Geminires$cell_type <- gsub("\\*\\*", "", Geminires$cell_type)
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
  Row = c(1, 2, 3, 5, 6, 8, 9, 10, 11),
  Marker_Gene = c(
    "GATA4, GATA6, PDGFRB, CK, HLA-G, ENG, TEAD4, CDX2, CGB5, ITGA5, ITGB1, MMP2, MMP9, VIM, ESRRB",
    "FOXJ1, DNAI1, DNAH5, DNAH9, TUBB4B, CCDC39, CCDC40, SPEF2, DYX1C1, RSPH4A, RSPH9, CFAP298, MUC1, MUC4, MUC16",
    "SOX10, PAX3, NGFR, MPZ, SNAI2, FOXD3, VWF, PECAM1, KDR, CD34, CDH1, CDH5, FOXJ1, DNAI2, TUBA1A",
    "PDGFRB, COL1A1, GFAP, LXN, ACTA2, ALB, AAT (SERPINA1), TTR, AFP, HNF4A, CYP3A4, CYP2E1, GATA4, HHEX, CK",
    "SOX2, POU5F1, NANOG, LIN28A, DAZL, VASA, DDX4, KLF4, SSEA4, TRA-1-60, TRA-1-81, FGF4, TDGF1, SALL4, ZFP42",
    "PPBP, PF4, ITGA2B, ITGB3, GP1BA, GP9, VWF, MPL, GATA1, THPO, ENG, CD42b, CD61, TBXAS1, FGA",
    "VIM, WT1, PAX2, T, SNAI1, SNAI2, FOXD3, SOX10, BRACHYURY, SHH, NKX6-1, NKX2-2, SOX9, NOG, BMP4",
    "CSH1, PLAC1, CGA, CK, CDX2, TEAD4, SOX2, NESTIN, PAX6, SOX9, FABP7, MSI1, VIM, NANOG, EOMES",
    "COL1A1, COL1A2, COL3A1, FAP, PDGFRB, THY1, CD44, ENG, CD73, VIM, S100A4, ACTA2, DDR2, TWIST1, THY1"
  )
)

cluster = c(1, 2, 3, 5, 6, 8, 9, 10, 11)

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
  row1 = c('PMP22','RGS2','COL3A1','CYP27A1','EMP2','BMPER','VIM','SNAI2','KRT8','SLC9A3R1','GATA6','ITGB1','VIM'),
  row2 = c('SPAG6','DNAH9','DRC1','DNAH5','DNALI1','C9orf24','TEKT1','CETN2','NRAV','OSCP1','DNAH5','DNAH9','RSPH4A','RSPH9','CFAP298'),
  row3 = c('MSX2','HERPUD1','TPM1','GLB1','ABCG2','KRT18','KRT7','TMEM54','DSP','DLX5','CDH1','TUBA1A'),
  row4 = c(),
  row5 = c('RBP4','APOA4','AGT','SERPINA1','FLRT3','CCKBR','IGFBP6','FGA','GPC3','APOE'),
  row6 = c('SNRPN','CTSC','SFRP2','GPRC5B','CYP2S1','KPNA2','UGP2','PHGDH','CDC20','DPPA4','POU5F1','LIN28A','TDGF1','ZFP42'),
  row7 = c(),
  row8 = c('THBS1','LCP1','PECAM1','AIF1','HLA.E','GPX1','B2M','MEF2C','S100A4','CCL3'),
  row9 = c('KCNK17','LINC01356','RSPO3','ABLIM1','TBXT','TMSB15A','CKS1B','PPIAP46','HMGB3','HVCN1'),
  row10 = c('SCG3','POU5F1','TBXT','L1TD1','FABP5','EFEMP1','ATP5PD','TSTD1','UCHL1','PTPRZ1'),
  row11 = c('SPARC','COLEC11','TGFBI','RGS16','IGF2','CXCL14','ANXA8','AQP1','DLK1','STX18','COL1A1','COL1A2','COL3A1','VIM')),
  negative_gene = list(
    row1 = c('GATA4','PDGFRB','ENG','TEAD4','CDX2','CGB5','ITGA5','MMP2','MMP9'),
    row2 = c('FOXJ1','DNAI1','TUBB4B','CCDC39','CCDC40','SPEF2','MUC4'),
    row3 = c('PAX3','NGFR','SNAI2','PECAM1','KDR','CD34'),
    row4 = c(),
    row5 = c('PDGFRB','COL1A1','GFAP','LXN','ACTA2','ALB','TTR','AFP','HNF4A','CYP2E1','GATA4','HHEX'),
    row6 = c('SOX2','NANOG','DAZL','KLF4','FGF4','SALL4'),
    row7 = c(),
    row8 = c('PPBP','PF4','ITGA2B','ITGB3','GP1BA','GP9','VWF','MPL','GATA1','THPO','ENG','TBXAS1','FGA'),
    row9 = c('VIM','SNAI1','SNAI2','FOXD3','SOX9','NOG','BMP4'),
    row10 = c('PLAC1','CGA','CDX2','TEAD4','SOX2','PAX6','SOX9','FABP7','MSI1','VIM','NANOG','EOMES'),
    row11 = c('FAP','PDGFRB','THY1','CD44','ENG','S100A4','ACTA2','DDR2','TWIST1','THY1')))

Geminires = inter_res_1[['Gemini']]
Geminires$cell_type <- gsub("\\*\\*", "", Geminires$cell_type)
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
write.csv(pn_gene,'./embyro_6_manual_pn_gene.csv')

MarkerPrompt(res1, 'ERNIE', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "APOA1, APOA2, APOB, FABP1, ALB, AFP, TTR, HNF4A, SOX17, FOXA2, CYP3A4, GATA6, HNF1A, CER1, GATA4",
    "FOXJ1, DNAI1, DNAH5, DNAH9, CCNO, MCIDAS, RSPH4A, RSPH9, DYNC2LI1, CFAP298, CFAP300, CFAP53, CFAP221, HYDIN, SPEF2",
    "SOX10, SNAI2, PAX3, FOXD3, EDN3, MITF, TFAP2A, SOX9, KIT, NGFR, MPZ, PHOX2B, RET, PLP1, ROBO1",
    "AFP, ALB, TAT, CYP3A4, HNF4A, FABP1, G6PC, HSD17B13, APOC3, SERPINA1, CYP7A1, CEBPA, ASGR1, ABCC2, BSEP",
    "SOX2, NANOG, POU5F1, LIN28A, CRIPTO, KLF4, DPPA4, TDGF1, SALL4, ZFP42, UTF1, LEFTY1, DPPA2, LIN28, GDF3",
    "HBA, HBB, KLF1, GYPA, EPB42, HBG1, GATA1, AHSP, ALAS2, UROS, ANK1, FECH, SPTA1, SPTB, HEMGN",
    "CD68, CD163, ITGAM, CSF1R, IL10, CCR2, CX3CR1, CD64, MRC1, CD14, TNF, FCER1G, MARCO, TLR4, SPI1",
    "MIXL1, EOMES, GSC, SOX17, T, FGF8, LEFTY1, NODAL, CER1, TDGF1, WNT3, TCF15, LHX1, TBXT, PDGFRB",
    "CSH1, CSH2, PL1, PL2, KRT7, EGFR, MMP2, MMP9, MMP12, CGA, CGB, LGALS14, TEAD4, HAND1, PSG1",
    "COL1A1, COL3A1, FAP, VIM, FN1, CD44, LUM, SPARC, PDGFRB, ACTA2, MMP2, MMP9, CTGF, S100A4, THY1"
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
write.csv(pn_gene,'./embyro_6_ERNIE_pn_gene.csv')



MarkerPrompt(res1, 'ChatGPT', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "VIM, ENG, FN1, CD44, MMP2, SNAI1, SNAI2, TWIST1, SOX9, FSP1, ACTA2, TGFBI, PDGFRB, LTBP3, THBS1",
    "FOXJ1, DNAI1, DNAH5, DNAH9, RSPH4A, CCNO, MCIDAS, RSPH9, CCDC39, CCDC40, HYDIN, CFAP53, CFAP298, SPEF2, LRRC6",
    "KRT8, KRT18, KRT19, CDH1, MUC1, TACSTD2, GATA3, ESR1, PGR, FOXA1, SPP1, AREG, BMP4, WNT7A, GREB1",
    "TBXT, MIXL1, GSC, MESP1, SNAIL, FGF8, EOMES, PDGFRA, WNT3, NT5E, HAND1, KDR, BMP4, EVX1",
    "SOX17, FOXA2, GATA4, GATA6, HNF1B, HNF4A, AFP, CXCR4, FABP1, APOA1, CER1, HHEX, SHH, PDGFRA, SFTA3",
    "DAZL, VASA, IFITM1, PRDM1, TFAP2C, NANOG, DND1, PLZF, OCT4, SALL4, KIT, LIN28, MVH, GDF3, SCP3",
    "HBA, HBB, GYPA, ANK1, EPB42, SPTA1, SPTB, CD47, RHAG, KLF1, GATA1, UROS, ALAS2, HBQ1, HBE1",
    "ITGAM, CD14, CD33, S100A9, MPO, CSF1R, CSF2RB, LYZ, CD11c, CD16, CD64, HLA-DRA, FLT3, LILRB4, CLEC4D",
    "NOG, SHH, FOXA2, TBXT, FGF8, ACTA1, ACAN, COL2A1, BRACHYURY, NKX3-2, FOXA1, CHRD, LMX1A, PAX1, TBX18",
    "NANOG, SOX2, EPIPL1, DPPA4, OTX2, TBXT, UTF1, POU5F1, ZFP42, KLF4, YAP1, NODAL, SALL4, CDH2, LEFTY1",
    "COL1A1, COL3A1, FAP, VIM, FN1, CD44, PDPN, PDGFRB, S100A4, DDR2, SPARC, ACTA2, BSP, LUM, MMP2"
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
write.csv(pn_gene,'./embyro_6_ChatGPT_pn_gene.csv')


MarkerPrompt(res1, 'Gemini', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "HLA-G, LGALS14, HSD3B1, CYP11A1, GCM1, KRT7, TACSTD2, MMP2, MMP9, ENG, VEGFA, FLT1, ACVR1B, PDGFA, EGFR",
    "FOXJ1, DNAI2, DNAH5, DNAH9, CCNO, MCIDAS, RSPH4A, CFAP53, CFAP298, CCDC39, HYDIN, SPAG6, STOM, PIH1D3, TTC25",
    "MITF, TYR, DCT, MC1R, PMEL, TYRP1, MLANA, KIT, GPR143, SOX10, SLC45A2, SLC24A5, RAB27A, SILV, EDNRB",
    "T, MESP1, MIXL1, FGF8, GSC, WNT3A, TBXT, FOXA2, LHX1, EOMES, PDGFRB, SNAIL1, FOXC2, NODAL, ENG",
    "SOX17, FOXA2, CXCR4, GATA4, GATA6, EPCAM, CER1, HNF1B, HNF4A, PDX1, KDR, CFC1, NODAL, LAMA1, SERPINA1",
    "DAZL, VASA, IFITM1, PRDM1, TFAP2C, NANOG, PLZF, SCP3, LIN28, OCT4, DND1, MVH, GDF3, TIAL1, KIT",
    "GATA1, EPOR, KLF1, HBZ, ALAS2, AHSP, HBA1, HBA2, HBB, HBG2, HBG1, ANK1, SPTA1, UROS, HBE1",
    "CD14, CD16, CD64, CCR2, S100A8, FCGR3A, CSF1R, LY86, ITGAM, HLA-DR, LYZ, MPO, CSF1, CD86, CD115",
    "BRACHYURY, NOG, NKX3-2, FOXA2, SHH, LMX1A, CHRD, ACAN, COL2A1, MMP2, CDH2, SPP1, GFAP, CTNNB1, FOXC2",
    "NANOG, SOX2, POU5F1, TDGF1, ZFP42, DPPA4, DPPA5, UTF1, SALL4, EOMES, FGF5, LEFTY1, NODAL, GDF3, LIN28",
    "VIM, FAP, CD90, PDGFRB, CD44, CXCL12, DKK3, LUM, ACTA2, MMP2, PDPN, THY1, ENG, SPARC, SPP1"
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
write.csv(pn_gene,'./embyro_6_Gemini_pn_gene.csv')


MarkerPrompt(res1, 'Llama', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "VIM, ENG, FN1, CD44, MMP2, SNAI1, SNAI2, TWIST1, SOX9, FSP1, ACTA2, TGFBI, PDGFRB, LTBP3, THBS1",
    "FOXJ1, DNAI1, DNAH5, DNAH9, CCNO, MCIDAS, RSPH4A, RSPH9, DYNC2LI1, CFAP298, CFAP300, CFAP53, CFAP221, HYDIN, SPEF2",
    "KRT8, KRT18, KRT19, CDH1, MUC1, TACSTD2, GATA3, ESR1, PGR, FOXA1, SPP1, AREG, BMP4, WNT7A, GREB1",
    "ALB, AFP, CYP3A4, HNF4A, TAT, TBIL, AAT, HNF1A, CYP2E1, CYP2C9, GGT1, FABP1, HSD17B4, CK8, CK18",
    "NANOG, OCT4, SOX2, KLF4, TDGF1, ZFP42, SALL4, LIN28, DPPA4, UTF1, LEFTY1, FGF4, GDF3, NANOS3, ESRRB",
    "GATA1, EPOR, KLF1, HBB, HBA1, ANK1, SPTA1, BAND3, GYPC, GYPB, HBE1, HBG1, HBG2, CABIN1, ALAS2",
    "CD41, CD61, MPL, PF4, CD42b, CD42c, CD42d, VWF, GPIX, F13A1, THPO, CLEC2, FGA, FGB, FGG",
    "SOX10, TFAP2A, SNAI2, EDNRB, FOXD3, PAX3, MITF, KIT, NGFR, MPZ, S100, PHOX2B, RET, SOX9, WNT1",
    "CGB, CGA, KRT7, GCM1, ERVW-1, TEAD4, TFAP2C, CDX2, HLA-G, TP63, ITGA6, ITGB1, SDC1, MMP9, HAND1",
    "COL1A1, COL3A1, S100A4, VIM, FN1, CD44, FAP, PDGFRB, SPARC, DDR2, FSP1, LUM, ACTA2, MMP2, CD90"
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
write.csv(pn_gene,'./embyro_6_Llama_pn_gene.csv')


MarkerPrompt(res1, 'Claude', species = "human")
markerdata <- data.frame(
  Row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  Marker_Genes = c(
    "HLA-G, ITGA5, ITGB1, MMP2, MMP9, CGB, ENG, INHA, KRT7, TACSTD2, TIMP1, VIM, LGALS14, DKK1, GCM1",
    "FOXJ1, DNAI1, DNAH5, DNAH9, TEKT1, CCNO, MCIDAS, RSPH4A, RSPH9, HYDIN, CFAP298, CFAP53, CCDC39, CCDC40, LRRC6",
    "CGB, CGA, KRT7, CDH1, TFAP2C, CYP19A1, ERVW-1, GCM1, EGFR, TEAD4, EPCAM, HAND1, ITGA6, HLA-G, CER1",
    "NANOG, SOX2, POU5F1, TDGF1, ZFP42, DPPA4, DPPA5, SALL4, GDF3, UTF1, LEFTY1, FGF5, LIN28, ESRRB, KLF4",
    "SOX7, SOX17, GATA4, GATA6, FOXA2, AFP, KRT19, HNF4A, KDR, TTR, CUBN, DAB2, PDGFRA, LAMA1, VEGFR2",
    "NANOG, SOX2, OCT4, LIN28, KLF4, SALL4, TDGF1, DPPA4, DPPA5, UTF1, LEFTY1, GDF3, FGF4, REX1, ESRRB",
    "GYPA, HBA1, HBA2, HBB, HBG1, HBZ, ANK1, EPB42, KLF1, GATA1, CD71, UROS, ALAS2, EPO, SPTA1",
    "CD68, ITGAM, CD163, CSF1R, CX3CR1, CD14, MRC1, CD64, TNF, IL10, TLR4, CCR2, SPI1, MARCO, MSR1",
    "T, TBXT, MIXL1, GSC, FOXA2, EOMES, BRACHYURY, NODAL, WNT3, CER1, SNAIL1, ENG, LHX1, FGF8, PDGFRB",
    "SOX2, T, TBXT, FGF8, NOTCH1, WNT3A, MSGN1, MEOX1, SNAI1, FOXC1, LHX1, PAX6, CD24, GDF11, CDH2",
    "VIM, ENG, CD44, CD90, S100A4, PDGFRB, CD105, PDGFRA, FAP, CXCL12, ACTA2, MMP2, LEP, DDR2, BMP2"
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
write.csv(pn_gene,'./embyro_6_Claude_pn_gene.csv')
