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

markers = readRDS('../markers/Mouse_SS_Fibro_markers.rds')

#######avg_log2FC
Fibro_res = LLMtissuetype(FindAllMarkersResult = markers,
                       species = 'mouse',
                       topgenenumber = 10,
                       celltype = 'fibroblast')

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
#Fibro_res[['Gemini']]$cell_type <- gsub("\\*\\*", "", Fibro_res[['Gemini']]$cell_type)
res = Fibro_res
#res[['Gemini']] = Geminires

res1 = create_dataframe(res, markers)
write.csv(res1,'./avg_log2FC_step1_res1.csv')
res1 = read.csv('./avg_log2FC_step1_res1.csv')

score = res1[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score",'all_model_score')]

average_scores <- colMeans(score)
average_scores = as.data.frame(average_scores)
write.csv(average_scores,'./avg_log2FC_average_scores_v1.csv')
res1 = res1[,-1]

next_top_10_markers <- markers %>%
  group_by(cluster) %>%  # 按照cluster分组
  arrange(desc(avg_log2FC)) %>%  # 按照avg_log2FC降序排列
  slice(11:20)

cluster_genes <- next_top_10_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(genes = paste(gene, collapse=",")) %>%
  dplyr::ungroup()
cluster_genes

mg <- character()
cluster = numeric()
row_number <- 1  # 初始化行号

for (i in 1:nrow(res1)) {
  print(i)
  n = res1[i,]  # 获取当前行数据
  if (n$all_model_score == 1 | n$all_model_score == 0){
    # 处理并确保以字符形式处理单元格内容，排除NA值
    ct <- unique(c(as.character(n$ERNIE), as.character(n$Gemini), as.character(n$ChatGPT), as.character(n$Claude), as.character(n$Llama)))
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

species = 'mouse'
# 生成最终输出
output_text <- paste('Provide key marker genes for the following ', species, ' cell types, with 15 key marker genes per cell type. Provide only the abbreviated gene names of key marker genes, full names are not required:\n', final_text, '\nThe format of the final response should be:\n\n1: gene1, gene2, gene3\n2: gene1, gene2, gene3\nN: gene1, gene2, gene3\n\n...where N represents the row number and gene1, gene2, gene3 represent key marker genes.', sep="")

print(output_text)

cluster

inter_res_1 = LLM_interect(positive_gene = list(
  row1 = c('Slc6a6','Igfbp4','Runx1','Pappa','Grem1','Tnc','Kitl','1500009L16Rik','Serpina12','Gas6','Spp1'),
  row2 = c('Col4a1','Nmb','Ctsh','Gdf10','Adamts5','Bmper','Gsn','Arl4a','Vwa1','Steap4','Fabp4','Lpl'),
  row3 = c('Smpd3','Has1','Il6','Hp','Thbd','Phlda1','Upk3b','Cebpb','Mt1','Ugdh'),
  row4 = c('Rgs5','Spon2','Lrat','Kcnk2','Ces1d','Grem1','Fmo2','Aldh1a3','Icam1','P2ry14','Klf4'),
  row5 = c('Tpt1','Tcf4','Crispld2','Imp3','Rps21','Fibin','Abca8a','Ccdc80','Lamc1','Rbp1','Bgn','Col1a1'),
  row6 = c('Klf4','Ly6a','Serpinf1','Clec3b','Fabp4','Meg3','Col1a2','Tnfaip6','Sparc','Sfrp2'),
  row7 = c('Lrat','Tgfbi','Vipr1','Ramp1','Ehd3','Gdf2','Rspo3','Ednrb','Abcc9','Apoc1','Pdgfrb','Rhoa'),
  row8 = c('Srgn','Csf2rb','Cpxm1','C3','Cxcl13','Ifitm1','H2-Ab1','Stmn2','Bst1','Dnajb1'),
  row9 = c('Mgp','Csrp1','Fmo2','Ppp1r14a','Rgs2','Lbh','Eln','Sod3','Selenbp1','Mamdc2','Sftpc'),
  row10 = c('Ccl2','Ifi203','Emp1','Rps27rt','Rps28','Nfkb1','Tsc22d2','Penk','Meg3','Trib1'),
  row11 = c('Ifi205','Fgl2','Ccl11','Emp1','Cilp','Serpinf1','Igfbp6','Itm2a','Hspa8','Tppp3'),
  row12 = c('Tmem100','Rps28','Rpl23a','Thy1','Prss23','Maf','Rps27','Rpl27','Rpl35','Penk'),
  row13 = c('Ly6c1','Slfn5','Flrt2','Tgfbr2','Marcks','Akr1c14','Socs1','Ly6a','Myoc','Ifi205'),
  row14 = c('Rdh10','Tlx1','Ifi27l2a','P2ry14','Stbd1','Mfge8','Clic5','Cd36','Frzb','Itga1','Col1a1','Timp3'),
  row15 = c('Lars2','Fmod','S100a4','Thbs4','Spp1','Col1a2','Cilp','Fn1','Col11a1','Anxa8'),
  row16 = c('Acta2','Fbln5','Wisp2','Mfap4','Cst3','Cd200','Cfh','Eln','Fmod','Pdlim3')),
  negative_gene = list(
    row1 = c('Col1a1','Bglap','Sp7','Runx2','Alpl','Ibsp','Acp5','Ctsk','Pth1r','Tnfrsf11b','Mmp13','Sox9','Tnfsf11'),
    row2 = c('Adipoq','Pparg','Plin1','Cebpa','Lep','Retn','Cd36','Fasn','Srebf1','Slc2a4','Npy1r','Adipor1'),
    row3 = c('Aqp2','Clcn5','Aqp1','Podxl','Lrp2'),
    row4 = c('Vil1','Muc2','Lyz1','Lgr5','Fabp2','Smoc2','Chga','Cldn3'),
    row5 = c('Col2a1','Acan','Sox9','Comp','Matn3','Col9a1','Col11a2','Col27a1','Fmod','Cilp'),
    row6 = c('Krt14','Krt5','Krt10','Krt1','Cdh3','Mlph','Sox9','Itga6'),
    row7 = c('Acta2','Myh11','Tagln','Col18a1','Vwf','Pecam1','Endog','Gja4','Egln3','Klf4','Nfatc1','Notch3'),
    row8 = c('Cd3e','Cd4','Cd8a','Cd19','Cd40','Foxp3','Ccr7','Bcl6','S1pr1','Lta','Ltb','Cxcl13','Cd79a'),
    row9 = c('Sftpb','Abca3','Scgb1a1','Foxj1','Aqp5','Muc5ac','Hif1a','Surf4','Ttf1','Scgb3a2','Pdpn','Napsa','Cdh1'),
    row10 = c('Alb','Cyp2e1','Fah','Rbp4','Ass1'),
    row11 = c('Acta1','Myod1','Dmd','Tnni2','Tnnt1','Tnnt3','Des','Ckm','Srf','Mstn','Myf5'),
    row12 = c('Ins1','Ins2','Cpa1','Cela3b','Spp1','Gata4'),
    row13 = c('Aqp2','Slc12a1','Slc13a3','Clcn5','Aqp1','Podxl','Calb1','Cubn'),
    row14 = c('Nfatc1','Gata4','Sox9','Runx2','Has2','Ctgf','Vcan','Smad6','Bmp4','Acta2','Postn'),
    row15 = c('Cd34','Cd38','Gata2','Runx1','Meis1','Cd44','Mpo'),
    row16 = c('Rho','Opn1mw','Gnat2','Cngb1','Sag','Gfap')))

Claude_input_result = c('1: Bone\n2: Adipose\n3: Kidney\n4: Lung\n5: Mammary\n6: Skin\n7: Bladder\n8: Lymph Node\n9: Lung\n10: Liver\n11: Muscle\n12: Pancreas\n13: Kidney\n14: Ligament\n15: Joint\n16: Eye')
rows <- strsplit(Claude_input_result, "\n")[[1]]
data <- sapply(rows, function(row) strsplit(row, ": ")[[1]])
df8 <- data.frame(clusters = as.integer(gsub(">", "", data[1, ]))-1, cell_type = data[2, ])
inter_res_1$Claude = df8
res = Fibro_res
Fibro_res_inter1 = create_dataframe(inter_res_1, markers)
Fibro_res_inter1$cluster = as.character(Fibro_res_inter1$cluster)

#write.csv(res1,'./avg_log2FC_step1_res1.csv')
Fibro_res_inter1$row = as.numeric(Fibro_res_inter1$clusters+1)
cluster_genes = readRDS('./cluster_genes.rds')
cluster_genes$row = as.numeric(cluster_genes$row)
Fibro_res_inter2 = Fibro_res_inter1%>%left_join(cluster_genes, by = c('row'))
Fibro_res_inter3 = read.csv('./avg_log2FC_step1_res1.csv')
Fibro_res_inter3 <- Fibro_res_inter3[ , -1]
Fibro_res_inter = cbind(Fibro_res_inter3, Fibro_res_inter2[,10:11], Fibro_res_inter2[,4:8])

write.csv(Fibro_res_inter,'./avg_log2FC_step1_res1_inter1.csv')




res1 = read.csv('../240618Mouse_SS_Fibro_different_method_third time/avg_log2FC_step1_res1_inter1.csv')
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
  output_text <- paste('Provide key marker genes for the following ', specie, ' tissue, with 15 key marker genes per cell type. Provide only the abbreviated gene names of key marker genes, full names are not required:\n', final_text, '\nThe format of the final response should be:\n\n1: gene1, gene2, gene3\n2: gene1, gene2, gene3\nN: gene1, gene2, gene3\n\n...where N represents the row number and gene1, gene2, gene3 represent key marker genes.', sep="")
  print(output_text)
}

res1$ERNIE_CLname
res1$Claude_CLname
MarkerPrompt(res1, 'cluster', species = "mouse")
MarkerPrompt(res1, 'ERNIE.1', species = "mouse")
MarkerPrompt(res1, 'ChatGPT.1', species = "mouse")
MarkerPrompt(res1, 'Gemini.1', species = "mouse")
MarkerPrompt(res1, 'Llama.1', species = "mouse")
MarkerPrompt(res1, 'Claude.1', species = "mouse")
