
#clear Global Environment
rm(list=ls())
library(ggprism)
library(reshape)
library(ggalluvial)

######fig2.b
library(gridExtra)
library(ggsci)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
v = pal_npg(palette = c("nrc"), alpha = 0.9)(9)
v[1] = 'orange'
coldf = c(v, v[6])
setwd('D:/R/LLMsCellType/')
#write.csv(res1, './240615pbmc_different_method/pvalue_adj_res_step2_v1.csv')
res1 = read.csv('./240615pbmc_different_method4/avg_log2FC_res_step2_v1.csv')
res2 = read.csv('./240615pbmc_different_method5/avg_log2FC_res_step2_v1.csv')
res3 = read.csv('./240615pbmc_different_method6/avg_log2FC_res_step2_v1.csv')
sc1 = res1[,34:45]
sc2 = res2[,34:45]
sc3 = res3[,34:45]
sc = rbind(sc1,sc2,sc3)
# 定义一个函数来找到每一行最大分数对应的细胞类型
sc$free_max_score <- do.call(pmax, c(sc[c("Gemini_score","Llama_score")], na.rm = TRUE))
get_Llama_Gemini_max_match <- function(row) {
  # 提取所有分数列
  scores <- c(row["Gemini_score"],row["Llama_score"])
  # 提取所有匹配列
  matches <- c(row["Gemini_match"], row["Llama_match"])
  # 找到最大分数的索引
  max_indices <- which(scores == max(scores, na.rm = TRUE))
  # 返回第一个最大分数对应的细胞类型
  matches[max_indices[1]]
}
# 应用函数到每一行
sc$free_max_match <- apply(res1, 1, get_Llama_Gemini_max_match)
score = sc[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score",'max_score','free_max_score')]
score[score == 1] = 0.5
score[score == 2] = 0.5
score[score == 3] = 1
score = score[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score")]
pbmc = colMeans(score)
#write.csv(sc, './pvalue_res_score_v1.csv')
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
  complete(Model, Score = 0:1, fill = list(Count = 0, Percentage = 0))

# 确保数据框和变量已经按前面的步骤准备好
score_distribution <- score_distribution %>%
  mutate(Label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), ""))  

score_distribution_summary <- score_distribution %>%
  group_by(Model) %>%
  summarise(
    Percentage0 = sum(Percentage[Score == 0]),
    Percentage3 = sum(Percentage[Score == 1])
  )

# 根据Score为0的百分比进行主要排序，如果相同，则按Score为3的百分比排序
ordering <- score_distribution_summary %>%
  arrange(Percentage3,desc(Percentage0)) %>%
  pull(Model)

# 使用factor函数对Model列重新设置级别，按Score为0的百分比排序
score_distribution$Model <- factor(score_distribution$Model, levels = ordering)

ggplot(score_distribution, aes( x = Model, y=Percentage, fill = factor(Score),
                 stratum = factor(Score), alluvium = factor(Score)))+
  geom_col(width = 0.65,
           color = '#8e8e8e', size = 1.5) +
  geom_flow(width = 0.65, alpha = 0.75, knot.pos = 0.3,
            color = '#8e8e8e', size = 1.5)+
  coord_flip()+
  scale_fill_manual(values = c("#efefef","#83dbca","#12846e"))+
  scale_y_continuous(expand = c(0,0),name="",
                     label=c("0%","25%","50%","75%","100%"))+
  scale_x_discrete(expand = c(0,0),name="")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color="black",size=13),
        axis.ticks.length.x = unit(0.1,"cm"),
        plot.margin = margin(10, 10, 10, 10)
  )+
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4)


ggsave(
  filename = './fig/figure2/pbmc1.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # 宽
  height = 4,            # 高
  units = "in",          # 单位
  dpi = 100)



######fig2.b

#write.csv(res1, './240615pbmc_different_method/pvalue_adj_res_step2_v1.csv')
res1 = read.csv('./240616gc_different_method4/avg_log2FC_res_step2_v1.csv')
res2 = read.csv('./240616gc_different_method5/avg_log2FC_res_step2_v1.csv')
res3 = read.csv('./240616gc_different_method6/avg_log2FC_res_step2_v1.csv')
sc1 = res1[,34:45]
sc2 = res2[,34:45]
sc3 = res3[,34:45]
sc = rbind(sc1,sc2,sc3)
# 定义一个函数来找到每一行最大分数对应的细胞类型
sc$free_max_score <- do.call(pmax, c(sc[c("Gemini_score","Llama_score")], na.rm = TRUE))
get_Llama_Gemini_max_match <- function(row) {
  # 提取所有分数列
  scores <- c(row["Gemini_score"],row["Llama_score"])
  # 提取所有匹配列
  matches <- c(row["Gemini_match"], row["Llama_match"])
  # 找到最大分数的索引
  max_indices <- which(scores == max(scores, na.rm = TRUE))
  # 返回第一个最大分数对应的细胞类型
  matches[max_indices[1]]
}

# 应用函数到每一行
sc$free_max_match <- apply(res1, 1, get_Llama_Gemini_max_match)
score = sc[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score",'max_score','free_max_score')]
score[score == 1] = 0.5
score[score == 2] = 0.5
score[score == 3] = 1
score = score[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score")]
gc = colMeans(score)
#write.csv(sc, './pvalue_res_score_v1.csv')
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
  complete(Model, Score = 0:1, fill = list(Count = 0, Percentage = 0))

# 确保数据框和变量已经按前面的步骤准备好
score_distribution <- score_distribution %>%
  mutate(Label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), ""))  

score_distribution_summary <- score_distribution %>%
  group_by(Model) %>%
  summarise(
    Percentage0 = sum(Percentage[Score == 0]),
    Percentage3 = sum(Percentage[Score == 1])
  )

# 根据Score为0的百分比进行主要排序，如果相同，则按Score为3的百分比排序
ordering <- score_distribution_summary %>%
  arrange(Percentage3,desc(Percentage0)) %>%
  pull(Model)

# 使用factor函数对Model列重新设置级别，按Score为0的百分比排序
score_distribution$Model <- factor(score_distribution$Model, levels = ordering)

ggplot(score_distribution, aes( x = Model, y=Percentage, fill = factor(Score),
                                stratum = factor(Score), alluvium = factor(Score)))+
  geom_col(width = 0.65,
           color = '#8e8e8e', size = 1.5) +
  geom_flow(width = 0.65, alpha = 0.75, knot.pos = 0.3,
            color = '#8e8e8e', size = 1.5)+
  coord_flip()+
  scale_fill_manual(values = c("#efefef","#83dbca","#12846e"))+
  scale_y_continuous(expand = c(0,0),name="",
                     label=c("0%","25%","50%","75%","100%"))+
  scale_x_discrete(expand = c(0,0),name="")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color="black",size=13),
        axis.ticks.length.x = unit(0.1,"cm"),
        plot.margin = margin(10, 10, 10, 10)
  )+
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4)


ggsave(
  filename = './fig/figure2/gc1.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # 宽
  height = 4,            # 高
  units = "in",          # 单位
  dpi = 100)

######fig2.c

#write.csv(res1, './240615pbmc_different_method/pvalue_adj_res_step2_v1.csv')
res1 = read.csv('./240617embryo_different_method4/avg_log2FC_res_step2_v1.csv')
res2 = read.csv('./240617embryo_different_method5/avg_log2FC_res_step2_v1.csv')
res3 = read.csv('./240617embryo_different_method6/avg_log2FC_res_step2_v1.csv')
sc1 = res1[,34:45]
sc2 = res2[,34:45]
sc3 = res3[,34:45]
sc = rbind(sc1,sc2,sc3)
# 定义一个函数来找到每一行最大分数对应的细胞类型
sc$free_max_score <- do.call(pmax, c(sc[c("Gemini_score","Llama_score")], na.rm = TRUE))
get_Llama_Gemini_max_match <- function(row) {
  # 提取所有分数列
  scores <- c(row["Gemini_score"],row["Llama_score"])
  # 提取所有匹配列
  matches <- c(row["Gemini_match"], row["Llama_match"])
  # 找到最大分数的索引
  max_indices <- which(scores == max(scores, na.rm = TRUE))
  # 返回第一个最大分数对应的细胞类型
  matches[max_indices[1]]
}

# 应用函数到每一行
sc$free_max_match <- apply(res1, 1, get_Llama_Gemini_max_match)
score = sc[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score",'max_score','free_max_score')]
score[score == 1] = 0.5
score[score == 2] = 0.5
score[score == 3] = 1
score = score[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score")]
embryo = colMeans(score)
#write.csv(sc, './pvalue_res_score_v1.csv')
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
  complete(Model, Score = 0:1, fill = list(Count = 0, Percentage = 0))

# 确保数据框和变量已经按前面的步骤准备好
score_distribution <- score_distribution %>%
  mutate(Label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), ""))  

score_distribution_summary <- score_distribution %>%
  group_by(Model) %>%
  summarise(
    Percentage0 = sum(Percentage[Score == 0]),
    Percentage3 = sum(Percentage[Score == 1])
  )

# 根据Score为0的百分比进行主要排序，如果相同，则按Score为3的百分比排序
ordering <- score_distribution_summary %>%
  arrange(Percentage3,desc(Percentage0)) %>%
  pull(Model)

# 使用factor函数对Model列重新设置级别，按Score为0的百分比排序
score_distribution$Model <- factor(score_distribution$Model, levels = ordering)

ggplot(score_distribution, aes( x = Model, y=Percentage, fill = factor(Score),
                                stratum = factor(Score), alluvium = factor(Score)))+
  geom_col(width = 0.65,
           color = '#8e8e8e', size = 1.5) +
  geom_flow(width = 0.65, alpha = 0.75, knot.pos = 0.3,
            color = '#8e8e8e', size = 1.5)+
  coord_flip()+
  scale_fill_manual(values = c("#efefef","#83dbca","#12846e"))+
  scale_y_continuous(expand = c(0,0),name="",
                     label=c("0%","25%","50%","75%","100%"))+
  scale_x_discrete(expand = c(0,0),name="")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color="black",size=13),
        axis.ticks.length.x = unit(0.1,"cm"),
        plot.margin = margin(10, 10, 10, 10)
  )+
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4)


ggsave(
  filename = './fig/figure2/embyro1.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # 宽
  height = 4,            # 高
  units = "in",          # 单位
  dpi = 100)



######fig2.d

#write.csv(res1, './240615pbmc_different_method/pvalue_adj_res_step2_v1.csv')
res1 = read.csv('./240618Mouse_SS_Fibro_different_method/avg_log2FC_final_v1.csv')
res2 = read.csv('./240618Mouse_SS_Fibro_different_method_second time/avg_log2FC_final_v1.csv')
res3 = read.csv('./240618Mouse_SS_Fibro_different_method_third time/avg_log2FC_final_v1.csv')
sc1 = res1[,10:16]
sc2 = res2[,10:16]
sc3 = res3[,10:16]
sc = rbind(sc1,sc2,sc3)
sc[sc==3]=1
# 定义一个函数来找到每一行最大分数对应的细胞类型
sc$free_max_score <- do.call(pmax, c(sc[c("Gemini_score","Llama_score")], na.rm = TRUE))
score = sc[c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score")]
fibro = colMeans(score)
#write.csv(sc, './pvalue_res_score_v1.csv')
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
  complete(Model, Score = 0:1, fill = list(Count = 0, Percentage = 0))

# 确保数据框和变量已经按前面的步骤准备好
score_distribution <- score_distribution %>%
  mutate(Label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), ""))  

score_distribution_summary <- score_distribution %>%
  group_by(Model) %>%
  summarise(
    Percentage0 = sum(Percentage[Score == 0]),
    Percentage3 = sum(Percentage[Score == 1])
  )

# 根据Score为0的百分比进行主要排序，如果相同，则按Score为3的百分比排序
ordering <- score_distribution_summary %>%
  arrange(Percentage3,desc(Percentage0)) %>%
  pull(Model)

# 使用factor函数对Model列重新设置级别，按Score为0的百分比排序
score_distribution$Model <- factor(score_distribution$Model, levels = ordering)

ggplot(score_distribution, aes( x = Model, y=Percentage, fill = factor(Score),
                                stratum = factor(Score), alluvium = factor(Score)))+
  geom_col(width = 0.65,
           color = '#8e8e8e', size = 1.5) +
  geom_flow(width = 0.65, alpha = 0.75, knot.pos = 0.3,
            color = '#8e8e8e', size = 1.5)+
  coord_flip()+
  scale_fill_manual(values = c("#efefef","#83dbca","#12846e"))+
  scale_y_continuous(expand = c(0,0),name="",
                     label=c("0%","25%","50%","75%","100%"))+
  scale_x_discrete(expand = c(0,0),name="")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color="black",size=13),
        axis.ticks.length.x = unit(0.1,"cm"),
        plot.margin = margin(10, 10, 10, 10)
  )+
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4)


ggsave(
  filename = './fig/figure2/fibro1.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # 宽
  height = 4,            # 高
  units = "in",          # 单位
  dpi = 100)


ggsave(
  filename = './fig/figure2/fibro.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 10,             # 宽
  height = 5,            # 高
  units = "in",          # 单位
  dpi = 100)


embryo = colMeans(score)

n = cbind(pbmc,gc,embryo,fibro)
n = as.data.frame(n)
n
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
library(pheatmap)

data <- tibble(
  pbmc = c(0.4462366, 0.4462366, 0.4516129, 0.4677419, 0.5107527),
  gc = c(0.5416667, 0.6388889, 0.6805556, 0.6805556, 0.7361111),
  embryo = c(0.1060606, 0.2727273, 0.1818182, 0.1666667, 0.1818182),
  fibro = c(0.2916667, 0.2500000, 0.2500000, 0.2500000, 0.2708333),
  model = c("ERNIE_score", "Gemini_score", "ChatGPT_score", "Llama_score", "Claude_score")
)
long_data <- data %>%
  pivot_longer(
    cols = -model,   # 指明不需要转换的列，这里是model列
    names_to = "tissue",  # 新列的名称，存放原来的列名
    values_to = "score"  # 新列的名称，存放原来的值
  )
long_data$model = factor(long_data$model, levels = c('ChatGPT_score',
                                               'Claude_score','Gemini_score','Llama_score','ERNIE_score'))
heatmap_plot <- ggplot(long_data, aes(x = tissue, y = model, fill = score)) +
  geom_tile() +  # 使用geom_tile创建热图
  geom_text(aes(label = sprintf("%.2f", score)), color = "black", angle = 0) +  # 添加文字标签
  scale_fill_gradientn(colors = c("#FEF0EA","#EA8568","#DD4633","#B41F23")) +  # 设置填充渐变色
  theme_bw() +  # 使用简洁主题
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +  # 旋转x轴标签以改善显示
  labs(x = "Model", y = "Data Type", fill = "Score")  # 添加标签
heatmap_plot
ggsave(
  filename = './fig/figure2/heatmap.pdf', # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 4,             # 宽
  height = 5,            # 高
  units = "in",          # 单位
  dpi = 100)
