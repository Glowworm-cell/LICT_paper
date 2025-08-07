setwd('D:/R/LLMsCellType/240618Mouse_SS_Fibro_different_method/')
detach(package:LLMTissueIdentifier, unload=TRUE)
install.packages("D:/R/LLMTissueIdentifier_0.1.0.tar.gz", repos = NULL, type = "source")
library(LLMTissueIdentifier)
library(reticulate)
library(rols)
library(dplyr)
library(scales)
library(Seurat)
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
markers = readRDS('./Mouse_SS_Fibro_markers.rds')

top10_pvalue_adj <- markers %>%
  group_by(cluster) %>%  # 按照cluster分组
  arrange(p_val_adj) %>%  # 按照p_val_adj降序排列
  slice(1:10)

Fibro_res = LLMtissuetype(FindAllMarkersResult = top10_pvalue_adj,
                       species = 'mouse',
                       topgenenumber = 10,
                       celltype = 'fibroblast')

res = Fibro_res

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
res1 = create_dataframe(res, top10_pvalue_adj)
write.csv(res1,'./pvalue_adj_step1_res1.csv')



#######avg_log2FC
Fibro_res = LLMtissuetype(FindAllMarkersResult = markers,
                       species = 'mouse',
                       topgenenumber = 10,
                       celltype = 'fibroblast')

#Fibro_res[['Gemini']]$cell_type <- gsub("\\*\\*", "", Fibro_res[['Gemini']]$cell_type)
res = Fibro_res
#res[['Gemini']] = Geminires

res1 = create_dataframe(res, markers)
write.csv(res1,'./avg_log2FC_step1_res1.csv')

#######Dvalue
markers$Dvalue = markers$pct.1-markers$pct.2
top10_Dvalue <- markers %>%
  group_by(cluster) %>%  # 按照cluster分组
  arrange(desc(Dvalue)) %>%  # 按照p_val_adj降序排列
  slice(1:10)
Fibro_res = LLMtissuetype(FindAllMarkersResult = top10_Dvalue,
                       species = 'mouse',
                       topgenenumber = 10,
                       celltype = 'fibroblast')

res = Fibro_res
#res[['Gemini']] = Geminires

res1 = create_dataframe(res, top10_Dvalue)
write.csv(res1,'./Dvalue_step1_res1.csv')

#######pvalue
top10_p_val <- markers %>%
  group_by(cluster) %>%  # 按照cluster分组
  arrange(p_val) %>%  # 按照p_val_adj降序排列
  slice(1:10)
Fibro_res = LLMtissuetype(FindAllMarkersResult = top10_p_val,
                       species = 'mouse',
                       topgenenumber = 10,
                       celltype = 'fibroblast')


res = Fibro_res
#res[['Gemini']] = Geminires

res1 = create_dataframe(res, top10_p_val)
write.csv(res1,'./pvalue_step1_res1.csv')

