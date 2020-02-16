# ##################### 代码整理 ##################### #
# 原地址：https://mp.weixin.qq.com/s?__biz=MjM5NTk0Mzg2Nw==&mid=2247484467&idx=1&sn=7a8c6744742249fdb25e83875597077d&chksm=a6f180e7918609f1ed7ba9736a0f226f59fb0d3895a73f5b1be24ad4eed53d5a6dd05eb10410&mpshare=1&scene=1&srcid=0202yITyJNw66kjtxYlWJDOg&sharer_sharetime=1580645669642&sharer_shareid=1400662730fdc372bf886bd03ea99487&key=de4a431374d95f8de65f730d8d922a2081d1c4ea8f0e15a3b3a0b7817a384b306f06058cad9ce186bd642a4ce3a8c7ab96c1f3474087834d4398c21417b5bdfc0a3dd189b1ae14d7225186e867f35e04&ascene=1&uin=MTIyNjEyMjQ2OA%3D%3D&devicetype=Windows+10&version=62070158&lang=zh_CN&exportkey=A2d8UqZiWszZXmLX%2Fd5%2BhEU%3D&pass_ticket=MS4VYevCquwI4TgK8UMhfYPuRTvz32R6%2BR55VUxGxhpIF8t%2B3MP0xuch%2Ba2NPAuC
# 作用：GEO数据的下载和分析，主要包括：PCA，火山图，MA图，热图，富集分析，相关性热图
# 处理：将平铺式代码整理为函数
# 整理者：向杨
# --------------------------------------------------
# 有用的函数使用说明：
# 1.下载数据：download_GSE_1(GSE_name); download_GSE_2(GSE_name);
# 2.获取探针-基因转换列表：get_probe_gene_map_1(GPL_name); get_probe_gene_map_2(GPL_name);get_probe_gene_map_3(GPL_name); get_probe_gene_map_4(GPL_name);
# 3.下载数据和探针转换的完整函数：download_and_annotion(GSE_name). 会在本地保存表达数据和分组数据的RData，命名分别为dat和group_list
# 4.PCA：pca(exp_RData_path)
# 5.热图：heatmap_top1000_sd(exp_RData_path)
# 6.相关性热图：cor_top500_mad(exp_RData_path)
# 7.火山图：volcano_MA(exp_RData_path)
# 8.GO和KEGG富集分析：anno_Go_Kegg(Deg_path)
# 9.GSEA：gsea(anno_Deg_path, gmt_path)
# 10.GSVA：gsva(exp_RData_path, gmt_path)
# 11.KEGG扩展图：kegg_expand_plots(anno_Deg_path)
# --------------------------------------------------
# 上述函数参数说明
# GSE_name：GSE名字
# GPL_name：GPL名字
# exp_RData_path：储存在本地的表达矩阵和分组信息的RData地址，数据可以由函数download_and_annotion(GSE_name)获得。
# Deg_path：logFC, t检验等数据的Rdata，可以由volcano_MA(exp_RData_path)获得
# anno_Deg_path：与Deg_path相似的Rdata，数据中多加了symbol等3列数据，可以由anno_Go_kegg(Deg_path)获得
# gmt_path：gmt文件路径，注意，这里是路径！！！！！
# -------------------------------------------------
# 过滤时的默认参数：pca(), cor_top500_mad(), heatmap_top1000_sd(), anno_Go_Kegg()
# exp_min：保留每一个基因在>5个人中表达量>exp_min
# exp_max：去除每一个基因在>5个人中表达量>exp_min
# num：提取的基因个数
# logFC：过滤的倍数值
# -------------------------------------------------
# -------------------------------------------------
## 简易使用, 下面代码直接出全部结果，只需要改GSE名字和gmt路径：
## source("代码整理_函数.R")                            # 可修改
## GSE_name = "GSE64634"                                # 可修改
## gmt_path = ""                                        # 可修改（gmt文件，可以自己制作或者gsea官网下载：msigdb.v7.0.symbols.gmt：https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/7.0/msigdb.v7.0.symbols.gmt； msigdb.v7.0.entrez.gmt：https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/7.0/msigdb.v7.0.entrez.gmt； 所有打包gmt：https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/7.0/msigdb_v7.0_files_to_download_locally.zip）
## exp_RData_path = "./data/step1-output.Rdata"
## Deg_path = "./data/deg.Rdata"
## anno_Deg_path = "./data/anno_DEG.Rdata"
## 这几个参数都是默认的保存名
## download_and_annotion(GSE_name)
## pca(exp_RData_path)
## heatmap_top1000_sd(exp_RData_path)
## cor_top500_mad(exp_RData_path)
## volcano_MA(exp_RData_path)
## anno_Go_Kegg(Deg_path)
## gsea(anno_Deg_path, gmt_path)                        # 网络不好，没下载到gmt文件，所以没试过是否可以运行
## gsva(exp_RData_path, gmt_path)                       # 网络不好，没下载到gmt文件，所以没试过是否可以运行
## kegg_expand_plots(anno_Deg_path)
# --------------------------------------------------
# --------------------------------------------------
# ################## 下面是函数 ################## #
# --------------------------------------------------
# --------------------------------------------------
# 下载GSE数据的2种方法：
# download_GSE_1(GSE_name); download_GSE_2(GSE_name);
download_GSE_1 = function(GSE_name){
  source('http://raw.githubusercontent.com/jmzeng1314/GEOmirror/master/R/geoChina.R') # 加载R代码
  gse = geoChina(gsename) # 下载gse数据，保存为Rdata，默认到当前工作目录
}

download_GSE_2 = function(GSE_name){
  library(GEOquery)
  gse = getGEO(GSE_name, destdir=".", getGPL=F)  # 数据默认到当前工作目录
}

# --------------------------------------------------
# 获取某一平台探针-基因的4种方法：
# get_probe_gene_map_1(GPL_name); get_probe_gene_map_2(GPL_name);
# get_probe_gene_map_3(GPL_name); get_probe_gene_map_4(GPL_name);
get_probe_gene_map_1 = function(GPL_name){
  install.packages("devtools")  # 包含下载github上的R包的函数
  library(devtools)  
  install_github("jmzeng1314/idmap1")  # 安装idmap1
  library(idmap1)
  ids = getIDs(GPL_name)
  return(ids)
}

get_probe_gene_map_2 = function(GPL_name){
  library(devtools)
  install_github("jmzeng1314/idmap2")
  library(idmap2)
  library(idmap2)
  ids=get_soft_IDs(GPL_name)
  return(ids)
}

get_probe_gene_map_3 = function(GPL_name){
  library(devtools)
  install_github("jmzeng1314/idmap3")
  library(idmap3)
  library(idmap3)
  ids=idmap3::get_pipe_IDs(GPL_name)
  return(ids)
}

get_probe_gene_map_4 = function(GPL_name){
  library(devtools)
  install_github("jmzeng1314/AnnoProbe")
  library(AnnoProbe)
  probe2gene=idmap(GPL_name, type = 'pipe')
  return(probe2gene)
}

# --------------------------------------------------
# GSE下载和探针转换的整体函数
# download_and_annotion(GSE_name)
download_and_annotion = function(GSE_name){
  # 下载
  fdata=paste0(GSE_name,"_eSet.Rdata")
  fpath=paste0("data/",fdata)
  if(!file.exists(fpath)){
    gset <- getGEO(GSE_name, destdir="data/",
                   AnnotGPL = F,     ## 注释文件
                   getGPL = F)       ## 平台文件
    gset=gset[[1]]
    save(gset,file=fpath)   ## 保存到本地
  }
  load(fpath)
  # ---------------------------------------
  # 分组
  pd=pData(gset)# 根据disease state:ch1一列得知分组信息
  classes = table(pd[, "disease state:ch1"])
  group_list=c(rep(names(classes)[1], classes[1]), rep(names(classes)[2], classes[2]))
  table(group_list)
  # ---------------------------------------
  # 对数据进行normalization
  dat=exprs(gset)
  dim(dat)
  
  dat[1:4,1:4]
  boxplot(dat,las=2)
  dat=dat[apply(dat,1,sd)>0,]# 去除都是0的探针
  dat[dat<0]=1
  boxplot(dat,las=2)
  
  #dat=log2(dat+1)
  #boxplot(dat,las=2)
  library(limma)
  dat=normalizeBetweenArrays(dat)
  boxplot(dat,las=2)
  # ----------------------------------------
  # 探针注释
  library(hgu133plus2.db)
  ids=toTable(hgu133plus2SYMBOL)
  head(ids)
  ids=ids[ids$symbol != '',]
  ids=ids[ids$probe_id %in%  rownames(dat),]# 过滤没法注释的探针
  
  dat[1:4,1:4]   
  dat=dat[ids$probe_id,]# 调整顺序，让dat的顺序和ids中的一致
  
  ids$median=apply(dat,1,median)
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]# 按照基因名、中位数大小排序
  ids=ids[!duplicated(ids$symbol),]# 只保留相同symbol中中位数最大的探针
  dat=dat[ids$probe_id,]# 调整顺序，让dat的顺序和ids中的一致
  rownames(dat)=ids$symbol# id转换
  dat[1:4,1:4]
  # -----------------------------------------
  # 保存
  save(dat,group_list,file = "./data/step1-output.Rdata")
}

# --------------------------------------------------
# PCA作图:
# pca(exp_RData_path)
pca = function(exp_RData_path, min_exp=1, max_exp=12){
  # 载入数据，数据为RData格式，包含表达数据和分组信息，分别为dat和group_list
  load(file=exp_RData_path)
  dat[1:4,1:4]
  # ---------------------------------------------
  # 数据过滤
  exprSet = dat
  dim(exprSet)
  # 过滤标准需要修改,目前是保留每一个基因在>5个人中表达量>1
  exprSet=exprSet[apply(exprSet,1, function(x) sum(x>min_exp) > 5),]
  boxplot(apply(exprSet,1 , mean))
  dim(exprSet)
  # 过滤标准需要修改,目前是去除每一个基因在>5个人中表达量>12
  exprSet=exprSet[!apply(exprSet,1, function(x) sum(x>max_exp) > 5),]
  dim(exprSet)
  # ---------------------------------------------
  # 作图
  exprSet=log(edgeR::cpm(exprSet)+1)
  dat=exprSet
  dat=as.data.frame(t(dat)) # 画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换。格式要求data.frame
  library("FactoMineR")# 计算PCA
  library("factoextra")# 画图展示
  
  dat.pca <- PCA(dat, graph = F)
  # fviz_pca_ind按样本  fviz_pca_var按基因
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # c("point", "text)2选1
               col.ind = group_list, # color by groups
               # palette = c("#00AFBB", "#E7B800"),# 自定义颜色
               addEllipses = T, # 加圆圈
               legend.title = "Groups"# 图例名称
  )
  if (!file.exists("./pic")){
    dir.create("./pic")
  }
  ggsave('pic/all_samples_PCA.png')
}

# --------------------------------------------------
# 挑选1000个SD最大的基因画表达量热图
# heatmap_top1000_sd(exp_RData_path)
heatmap_top1000_sd = function(exp_RData_path, num=1000){
  # 载入数据
  load(file = exp_RData_path)
  dat[1:4,1:4]
  # -----------------------------------------------
  # 作图
  cg=names(tail(sort(apply(dat,1,sd)),num))# 找到SD最大的1000个基因
  library(pheatmap)
  headmap.dat=dat[cg,]
  # 把每个基因在不同处理和重复中的数据转换为平均值为0，
  # 方差为1的数据，所以在这里也需要先转置(针对基因)。
  headmap.dat=t(scale(t(headmap.dat)))
  headmap.dat[headmap.dat>2]=2 
  headmap.dat[headmap.dat< -2]= -2
  
  # 分组信息设置
  ac=data.frame(group=group_list)
  rownames(ac)=colnames(headmap.dat)
  
  ## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
  pheatmap(headmap.dat,show_colnames =T,show_rownames = F,
           annotation_col=ac,
           filename = 'pic/heatmap_top1000_sd.png')
}

# --------------------------------------------------
# top_500_MAD基因和所有基因相关性热图
# cor_top500_mad(exp_RData_path)
cor_top500_mad = function(exp_RData_path, num=500, min_exp=1, max_exp=12){
  # 载入数据
  load(file=exp_RData_path) 
  dat[1:4,1:4] 
  exprSet=dat
  # -------------------------------------------------------
  # 所有基因相关性图
  ac=data.frame(group=group_list)
  rownames(ac)=colnames(exprSet)
  pheatmap::pheatmap(cor(exprSet),
                     annotation_col = ac,
                     show_rownames = F,
                     filename = 'pic/cor_all_genes.png')
  # -------------------------------------------------------
  # 数据过滤
  dim(exprSet)
  # 过滤标准需要修改,目前是保留每一个基因在>5个人中表达量>1
  exprSet=exprSet[apply(exprSet,1, function(x) sum(x>min_exp) > 5),]
  boxplot(apply(exprSet,1 , mean))
  dim(exprSet)
  # 过滤标准需要修改,目前是去除每一个基因在>5个人中表达量>12
  exprSet=exprSet[!apply(exprSet,1, function(x) sum(x>max_exp) > 5),]
  dim(exprSet)
  # -------------------------------------------------------
  # 数据normalization处理
  # 去除文库大小差异normalization，同时利用log做scale
  exprSet=log(edgeR::cpm(exprSet)+1)
  # 取top500 mad 的基因画图
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:num]),]
  # -------------------------------------------------------
  # top_500作图
  M=cor(exprSet) 
  pheatmap::pheatmap(M,
                     show_rownames = F,
                     annotation_col = ac,
                     filename = 'pic/cor_top500_mad.png')
}

# --------------------------------------------------
# 火山图 和MA 图
# volcano_MA(exp_RData_path)
volcano_MA = function(exp_RData_path){
  # 载入数据
  options(stringsAsFactors = F)
  library(ggpubr)
  library(limma)
  load(file = exp_RData_path)
  # 每次都要检测数据
  dat[1:4,1:4] 
  table(group_list)
  # ----------------------------------------------
  #按照group_list分组画箱线图
  boxplot(dat[1,]~group_list) 
  
  # boxplot的美化版
  bplot=function(g){
    df=data.frame(gene=g,stage=group_list)
    p <- ggboxplot(df, x = "stage", y = "gene",
                   color = "stage", palette = "jco",
                   add = "jitter")
    #  Add p-value
    p + stat_compare_means()
  }
  
  # ----------------------------------------------
  # 制作比较矩阵
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  head(design)
  exprSet=dat
  rownames(design)=colnames(exprSet)
  design
  
  # 比较矩阵
  # 这个矩阵声明，我们要把 npc 组跟 Normal 进行差异分析比较
  contrast.matrix<-makeContrasts("npc-normal",
                                 levels = design)
  contrast.matrix
  
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix)
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    
    ##step3
    # 有了比较矩阵后，coef=1，而number=Inf是把所有结果都打印出来
    tempOutput = topTable(fit2, coef=1, number =Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  deg = deg(exprSet,design,contrast.matrix)
  save(deg,file = './data/deg.Rdata')
  
  # ------------------------------------------
  # 作图
  load('./data/deg.Rdata')
  head(deg)
  bplot(dat[rownames(deg)[1],])
  ## for volcano and MA plot
  # 结果存放在pic/volcano.png和pic/MA.png
  if(T){
    nrDEG=deg
    head(nrDEG)
    attach(nrDEG)
    # 原始版火山图
    plot(logFC,-log10(P.Value))
    library(ggpubr)
    df=nrDEG
    df$y= -log10(P.Value)
    ggscatter(df, x = "logFC", y = "y",size=0.5)
    # 定义logFC=2为阈值
    df$state=ifelse(df$P.Value>0.01,'stable',
                    ifelse( df$logFC >2,'up',
                            ifelse( df$logFC < -2,'down','stable') )
    )
    table(df$state)
    df$name=rownames(df)
    head(df)
    ggscatter(df, x = "logFC", y = "y",size=0.5,color = 'state')
    ggscatter(df, x = "logFC", y = "y", color = "state",size = 0.5,
              label = "name", repel = T,
              #label.select = rownames(df)[df$state != 'stable'] ,
              label.select = c('TTC9', 'AQP3', 'CXCL11','PTGS2'), #挑选一些基因在图中显示出来
              palette = c("#00AFBB", "#E7B800", "#FC4E07"))
    ggsave('pic/volcano.png')
    
    # MA图
    ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
    df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                    ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
    table(df$p_c )
    ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
              palette = c("green", "red", "black") )
    ggsave('pic/MA.png')
  }
  
  ## for heatmap 
  if(T){ 
    load(file = 'data/step1-output.Rdata')
    # 每次都要检测数据
    dat[1:4,1:4]
    table(group_list)
    x=deg$logFC
    names(x)=rownames(deg)
    
    # cg中存放着变化上升和下降的前100个基因名
    cg=c(names(head(sort(x),100)),
         names(tail(sort(x),100)))
    library(pheatmap)
    n=t(scale(t(dat[cg,])))
    
    n[n>2]=2
    n[n< -2]= -2
    n[1:4,1:4]
    pheatmap(n,show_colnames =F,show_rownames = F)
    ac=data.frame(group=group_list)
    rownames(ac)=colnames(n) #将ac的行名也就分组信息（是‘no TNBC’还是‘TNBC’）给到n的列名，即热图中位于上方的分组信息
    pheatmap(n,show_colnames =F,
             show_rownames = F,
             cluster_cols = F, 
             annotation_col=ac,filename = 'pic/heatmap_top200_DEG.png') #列名注释信息为ac即分组信息
    
    
  }
  
  write.csv(deg,file = 'data/deg.csv')
}

# --------------------------------------------------
# KEGG pathway analysis
# 具体结果数据在data/下，图在pic/下
run_kegg <- function(gene_up,gene_down,geneList=F,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  ###   over-representation test
  # 下面把3个基因集分开做超几何分布检验
  # 首先是上调基因集。
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk=kk.up
  dotplot(kk)
  barplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  head(kk)
  write.csv(kk@result,paste0('data/',pro,'_kk.up.csv'))
  
  # 首先是下调基因集。
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  kk=kk.down
  dotplot(kk)
  barplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0('data/',pro,'_kk.down.csv'))
  
  # 最后是上下调合并后的基因集。
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  kk=kk.diff
  dotplot(kk)
  barplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0('data/',pro,'_kk.diff.csv'))
  
  
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.01,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.01,];up_kegg$group=1
  #画图设置, 这个图很丑，大家可以自行修改。
  g_kegg = kegg_plot(up_kegg,down_kegg)
  
  ggsave(g_kegg,filename = paste0('pic/',pro,'_kegg_up_down.png') )
  
  ###  GSEA 
  if(geneList){
    
    ## GSEA算法跟上面的使用差异基因集做超几何分布检验不一样。
    kk_gse <- gseKEGG(geneList     = geneList,
                      organism     = 'hsa',
                      nPerm        = 1000,
                      minGSSize    = 20,
                      pvalueCutoff = 0.9,
                      verbose      = FALSE)
    head(kk_gse)[,1:6]
    gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
    gseaplot(kk_gse, 'hsa04110',title = 'Cell cycle') 
    kk=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    tmp=kk@result
    write.csv(kk@result,paste0(pro,'_kegg.gsea.csv'))
    
    
    # 这里找不到显著下调的通路，可以选择调整阈值，或者其它。
    down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
    up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
    
    g_kegg=kegg_plot(up_kegg,down_kegg)
    print(g_kegg)
    ggsave(g_kegg,filename = paste0(pro,'_kegg_gsea.png'))
  }
}

# GO database analysis 
# 具体结果数据在data/下，图在pic/下
run_go <- function(gene_up,gene_down,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)
  # 因为go数据库非常多基因集，所以运行速度会很慢。
  if(T){
    go_enrich_results <- lapply(g_list, function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',names(gene),ont ))
        ego <- enrichGO(gene          = gene,
                        #universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file =paste0('data/',pro, '_go_enrich_results.Rdata'))
    
  }
  # 只有第一次需要运行，就保存结果哈，下次需要探索结果，就载入即可。
  # load(file=paste0('data/',pro, '_go_enrich_results.Rdata'))
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0("pic/",pro, '_dotplot_',n1[i],'_',n2[j],'.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }
}

# 合并up和down的基因KEGG结果，放在一幅图里展示
kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment")
  return(g_kegg)
}


# Go和kegg注释
# anno_Go_kegg(Deg_path)
anno_Go_Kegg = function(Deg_path, logFC=1.5){
  # 载入数据
  load(file=Deg_path)
  head(deg)
  # ------------------------------------------------
  logFC_t = logFC
  # 预处理1
  deg$g=ifelse(deg$P.Value>0.05,'stable',
               ifelse( deg$logFC > logFC_t,'UP',
                       ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
  )
  table(deg$g)
  head(deg)
  deg$symbol=rownames(deg)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  head(df)
  DEG=deg
  head(DEG)
  
  DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
  head(DEG)
  save(DEG,file = 'data/anno_DEG.Rdata')
  
  # -------------------------------------------------
  # 预处理2
  gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
  gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
  gene_diff=c(gene_up,gene_down)
  gene_all=as.character(DEG[ ,'ENTREZID'] )
  data(geneList, package="DOSE") 
  head(geneList)
  boxplot(geneList)
  boxplot(DEG$logFC)
  
  geneList=DEG$logFC
  names(geneList)=DEG$ENTREZID
  geneList=sort(geneList,decreasing = T)
  # -----------------------------------------------
  # 画图
  # source('kegg_and_go_up_and_down.R')
  run_kegg(gene_up,gene_down,pro='npc_VS_normal')
  # 需要多go数据库的3个条目进行3次富集分析，非常耗时。
  run_go(gene_up,gene_down,pro='npc_VS_normal')
  # ----------------------------------------------
  # 综合显示
  go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
  library(ggplot2)
  library(stringr)
  barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
  barplot(go, split="ONTOLOGY",font.size =10)+ 
    facet_grid(ONTOLOGY~., scale="free") + 
    scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
    ggsave('pic/gene_up_GO_all_barplot.png') 
  
  go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
  barplot(go, split="ONTOLOGY",font.size =10)+ 
    facet_grid(ONTOLOGY~., scale="free") + 
    scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
    ggsave('pic/gene_down_GO_all_barplot.png')
}

# --------------------------------------------------
# GSEA
# gsea(anno_Deg_path, gmt_path)
gsea = function(anno_Deg_path, gmt_path){
  # 载入数据和R包
  options(stringsAsFactors = F)
  load(file = anno_Deg_path)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  # ------------------------------------------------
  # GSEA
  # 对 MigDB中的全部基因集 做GSEA分析。
  # 按照FC的值对差异基因进行排序
  # http://www.bio-info-trainee.com/2105.html
  # http://www.bio-info-trainee.com/2102.html 
  # 自行修改存放gmt文件路径d
  # GSEA每个gene set的具体结果保存在gsea_results这个list中
  # 而最终结果保存在gsea_results_df数据框中
  # d='data/GSEA_gmt/msigdb_v7.0_files_to_download_locally/msigdb_v7.0_GMTs/symbols/'
  if(T){
    geneList=DEG$logFC
    names(geneList)=DEG$symbol
    geneList=sort(geneList,decreasing = T)
    #选择gmt文件（MigDB中的全部基因集）
    
    gmts=list.files(gmt_path, pattern = 'all')
    gmts
    
    #GSEA分析
    library(GSEABase) # BiocManager::install('GSEABase')
    ## 下面使用lapply循环读取每个gmt文件，并且进行GSEA分析
    ## 如果存在之前分析后保存的结果文件，就不需要重复进行GSEA分析。
    f='data/gsea_results.Rdata'
    if(!file.exists(f)){
      gsea_results <- lapply(gmts, function(gmtfile){
        # gmtfile=gmts[2]
        filepath=paste0(d,gmtfile)
        geneset <- read.gmt(filepath) 
        print(paste0('Now process the ',gmtfile))
        egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
        head(egmt)
        # gseaplot(egmt, geneSetID = rownames(egmt[1,]))
        return(egmt)
      })
      # 上面的代码耗时，所以保存结果到本地文件
      save(gsea_results,file = f)
    }
    load(file = f)
    # 提取gsea结果，熟悉这个对象
    gsea_results_list<- lapply(gsea_results, function(x){
      cat(paste(dim(x@result)),'\n')
      x@result
    })
  }
}

# --------------------------------------------------
# gsva
# gsva(exp_RData_path, gmt_path)
gsva = function(exp_RData_path, gmt_path){
  # 载入数据
  options(stringsAsFactors = F)
  
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  load(file=exp_RData_path)
  # -----------------------------------------------
  # gsva
  # GSVA分析
  # 存放gene set的文件路径需要具体修改
  d = gmt_path
  X = dat
  table(group_list)
  ## Molecular Signatures Database (MSigDb) 
  gmts = list.files(d,pattern = 'all')
  gmts
  library(GSVA) # BiocManager::install('GSVA')
  if(!file.exists('data/gsva_msigdb.Rdata')){
    es_max <- lapply(gmts, function(gmtfile){ 
      # gmtfile=gmts[8];gmtfile
      geneset <- read.gmt(file.path(d,gmtfile))  
      es.max <- gsva(X, geneset, 
                     mx.diff=FALSE, verbose=FALSE, 
                     parallel.sz=1)
      return(es.max)
    })
    adjPvalueCutoff <- 0.001
    logFCcutoff <- log2(2)
    es_deg <- lapply(es_max, function(es.max){
      # es.max=es_max[[1]]
      table(group_list)
      dim(es.max)
      design <- model.matrix(~0+factor(group_list))
      colnames(design)=levels(factor(group_list))
      rownames(design)=colnames(es.max)
      design
      library(limma)
      # contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
      contrast.matrix<-makeContrasts("npc-normal",
                                     levels = design)
      
      contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
      
      deg = function(es.max,design,contrast.matrix){
        ##step1
        fit <- lmFit(es.max,design)
        ##step2
        fit2 <- contrasts.fit(fit, contrast.matrix) 
        ##这一步很重要，大家可以自行看看效果
        
        fit2 <- eBayes(fit2)  ## default no trend !!!
        ##eBayes() with trend=TRUE
        ##step3
        res <- decideTests(fit2, p.value=adjPvalueCutoff)
        summary(res)
        tempOutput = topTable(fit2, coef=1, n=Inf)
        nrDEG = na.omit(tempOutput) 
        #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
        head(nrDEG)
        return(nrDEG)
      }
      
      re = deg(es.max,design,contrast.matrix)
      nrDEG=re
      head(nrDEG) 
      return(nrDEG)
    })
    gmts
    save(es_max,es_deg,file='data/gsva_msigdb.Rdata')
  }
  # --------------------------------------------------
  # 画图展示，结果存放在pic/下
  
  load(file='data/gsva_msigdb.Rdata')
  library(pheatmap)
  lapply(1:length(es_deg), function(i){
    # i=8
    print(i)
    dat=es_max[[i]]
    df=es_deg[[i]]
    df=df[df$P.Value<0.01 & abs(df$logFC) > 0.3,]
    print(dim(df))
    if(nrow(df)>5){
      n=rownames(df)
      dat=dat[match(n,rownames(dat)),]
      ac=data.frame(g=group_list)
      rownames(ac)=colnames(dat)
      rownames(dat)=substring(rownames(dat),1,50)
      pheatmap::pheatmap(dat, 
                         fontsize_row = 8,height = 11,
                         annotation_col = ac,show_colnames = F,
                         filename = paste0('[pic/gsva_',strsplit(gmts[i],'[.]')[[1]][1],'.pdf'))
      
    }
  })
  
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  df=do.call(rbind ,es_deg)
  es_matrix=do.call(rbind ,es_max)
  df=df[df$P.Value<0.01 & abs(df$logFC) > 0.5,]
  write.csv(df,file = 'data/GSVA_DEG.csv')

}

# --------------------------------------------------
# kegg 扩展
# kegg_expand_plots(anno_Deg_path)
kegg_expand_plots = function(anno_Deg_path){
  ### 主要是关于KEGG方面的扩展图
  ### 主要是关于KEGG方面的扩展图
  
  # 载入数据
  if (T) {
    rm(list = ls()) 
    options(stringsAsFactors = F)
    
    library(ggplot2)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    load(file = anno_Deg_path)  
    
    head(DEG)
  }
  
  # pre-process data
  if (T) {
    gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
    gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
    gene_diff=c(gene_up,gene_down)
    gene_all=as.character(DEG[ ,'ENTREZID'] )
    # 制作差异基因list L
    boxplot(DEG$logFC)
    
    geneList=DEG$logFC
    names(geneList)=DEG$ENTREZID
    geneList=sort(geneList,decreasing = T)
  }
  
  
  # KEGG富集分析得到结果
  if (T) {
    if (!file.exists("data/enrichkk.rdata")) {
      gene_down
      gene_up
      enrichKK <- enrichKEGG(gene         =  gene_up,
                             organism     = 'hsa',
                             #universe     = gene_all,
                             pvalueCutoff = 0.1,
                             qvalueCutoff =0.1)
      save(enrichKK,file = "data/enrichkk.rdata")
    }
    load(file = "data/enrichkk.rdata")
    # 查看KEGG结果
    head(enrichKK)[,1:6] 
    # 打开网页看相关KEGG通路图
    browseKEGG(enrichKK, 'hsa05150')
    
    # 将数据中的entrz-id变成symbol
    # 更为易读
    enrichKK=DOSE::setReadable(enrichKK, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    enrichKK 
  }
  
  
  ## 可视化
  #条带图
  if (T) {
    # par(mfrow=c(2,1))
    barplot(enrichKK,showCategory=20)
    ggsave("pic/barplot.png")
  }
  
  #气泡图
  if (T) {
    dotplot(enrichKK)
    ggsave("pic/dotplot.png")
  }
  
  #下面的图需要映射颜色，设置和示例数据一样的geneList
  
  # 展示top5通路的共同基因，要放大看。
  #Gene-Concept Network
  if (T) {
    cnetplot(enrichKK, foldChange=geneList,colorEdge = TRUE, circular = F)
    ggsave("pic/cnetplot.png")
    cnetplot(enrichKK, foldChange=geneList, colorEdge = TRUE, circular = T)
    ggsave("pic/cnetplot_circular.png")
  }
  
  
  #Enrichment Map
  if (T) {
    emapplot(enrichKK)
    ggsave("pic/Enrichment_Map.png")
  }
  
  #(4)展示通路关系,仅仅是针对于GO数据库结果。
  # goplot(enrichKK)
  #(5)Heatmap-like functional classification
  if (T) {
    heatplot(enrichKK,foldChange = geneList)
    ggsave("pic/Enrichment_Heatmap.png")
  }
}









