#recovery
UMAP_cluster <- function(x,  pcs = 20, label = NULL, res = 1) {
  
  x.seurat <- CreateSeuratObject(counts = x)
  x.seurat <- NormalizeData(x.seurat)
  x.seurat <- ScaleData(x.seurat)
  x.seurat <- FindVariableFeatures(x.seurat, do.plot = FALSE)
  x.seurat <- RunPCA(x.seurat,  pc.genes = x.seurat@var.genes, pcs.compute = pcs, do.print = FALSE)
  x.seurat <-RunUMAP(x.seurat, dims = 1:10)
  if(is.null(label)){
    x.seurat <- FindNeighbors(x.seurat, dims = 1:pcs)
    x.seurat <- FindClusters(x.seurat, dims.use = 1:pcs, resolution = res,
                             print.output = FALSE, save.SNN = TRUE,graph.name = NULL)
  }else{
    x.seurat@ident =  as.factor(label)
  }
  
  # if (is.null(x.seurat[['tsne']])) {
  #   x.seurat <- RunTSNE(x.seurat, dims.use = 1:pcs, check_duplicates = FALSE,
  #                       do.fast = TRUE)
  # }
  
  x.seurat
}
# x,  pcs = 20, label = NULL, res = 1
#   
# x.seurat <- CreateSeuratObject(counts = x)
# x.seurat <- NormalizeData(x.seurat)
# x.seurat <- ScaleData(x.seurat)
# x.seurat <- FindVariableFeatures(x.seurat, do.plot = FALSE)
# x.seurat <- RunPCA(x.seurat,  pc.genes = x.seurat@var.genes, pcs.compute = pcs, do.print = FALSE)
# 
# if(is.null(label)){
#   x.seurat <- FindNeighbors(x.seurat, dims = 1:pcs)
#   x.seurat <- FindClusters(x.seurat, dims.use = 1:pcs, resolution = res,
#                            print.output = FALSE, save.SNN = TRUE,graph.name = NULL)
# }else{
#   x.seurat@ident =  as.factor(label)
# }

# if (is.null(x.seurat[['tsne']])) {
#   x.seurat <- RunTSNE(x.seurat, dims.use = 1:pcs, check_duplicates = FALSE,
#                       do.fast = TRUE)
# }
# enbaron <- readRDS("D:\\系统默认\\桌面\\后续实验\\revise_exp\\enbaron.rds")
# chu_EC=enbaron[["count.imputed.individual"]][,,1]
# baron <- readRDS("D:\\系统默认\\桌面\\后续实验\\需要插补的数据集\\baron.rds")
# chu_EC=baron[["count.samp"]]
# chu_EC <- readRDS("D:\\系统默认\\桌面\\后续实验\\需要插补的数据集\\chu_EC.rds")
# chu_EC <- read.csv("D:\\系统默认\\桌面\\aaa.csv",header=TRUE,row.names = 1,check.names = FALSE)
# deng <- readRDS("D:\\系统默认\\桌面\\后续实验\\需要插补的数据集\\deng_count_pre.rds")
# count.tsne = tsne_cluster(chu_EC)
# p1=DimPlot(count.tsne, reduction = "umap")
# p2=p1
# p1+p2
# library(Seurat)
# plot(
#       count.tsne@reductions[["pca"]]@cell.embeddings,
#       col = as.numeric(count.tsne@active.ident),
#       pch = 19,
#       axes = FALSE,
#       frame.plot = TRUE,
#       cex = 0.6,
#       main = "Before EnImpute"
#     )


#cluster
library(SingleCellExperiment)
library(SC3)
library(mclust)
library(aricode)
library(Seurat)
seurat_create=function(result){
  x.seurat <- CreateSeuratObject(result)
  x.seurat <- NormalizeData(x.seurat)
  x.seurat <- ScaleData(x.seurat)
  x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)
  
  x.seurat <- RunPCA(x.seurat, features = VariableFeatures(object = x.seurat))
  x.seurat <- JackStraw(x.seurat, num.replicate = 100)
  x.seurat <- ScoreJackStraw(x.seurat, dims = 1:20)
  x.seurat <- FindNeighbors(x.seurat, dims = 1:10)
  x.seurat <- FindClusters(x.seurat, resolution = 0.5)
  x.seurat
}
seurat_cal=function(label,x.seurat){
  ARI <- adjustedRandIndex(as.factor(label), x.seurat$seurat_clusters)
  NMI <- NMI(x.seurat$seurat_clusters, as.factor(label))
  score=data.frame(criteria=c("ARI","NMI"),value=c(ARI,NMI))
  score
}
clusterscore=function(count){
  count.seurat=seurat_create(count)
  label = colnames(count)
  seurat_count_score=seurat_cal(label,count.seurat)
  seurat_count_score
}
# count = read.csv(
#   file = "D:\\系统默认\\桌面\\后续实验\\需要插补的数据集\\pbmc_sw.csv",
#   header = TRUE,
#   row.names = 1, check.names=FALSE
# )
# count = as.matrix(count)
# pbmc_sw <- readRDS("D:\\系统默认\\桌面\\后续实验\\需要插补的数据集\\pbmc_sw.rds")
# a=clusterscore(count)
# label = colnames(count)
# library(ggplot2)
# ggplot(a,aes(x=criteria,y=value,fill=criteria))+
#   geom_bar(stat = "identity")+
#   theme(legend.position ='none')+xlab('Before EnImpute2')+ylab('')
# library(ggplot2)
# k=data.frame(criteria=c("ARI","NMI"),value=c(0.4,0.5))
# ggplot(k,aes(x=criteria,y=value,fill=criteria))+
#   geom_bar(stat = "identity")+
#   theme(legend.position ='none')+xlab('aaa')+theme_classic()


#differential
# library(edgeR)
library(ggplot2)
# my.edgeR <- function(count, group){
#   intcol<-vector("logical")
#   for(i in 1:ncol(count)){
#     intcol<-c(intcol,is.integer(count[,i]))
#   }
#   if (!all(TRUE==intcol)) {
#     warning("WARNING! Non-integer expression levels. Data rounded")
#     count<-round(count)
#   }
#   count.sel = count
#   d = DGEList(count.sel, group=group)
# 
#   d = calcNormFactors(d)
# 
#   d = estimateCommonDisp(d)
#   d = estimateTagwiseDisp(d)
# 
#   dest = exactTest(d)
#   
#   dest.sort = topTags(dest, n = dim(dest)[1])
#   res.edgeR <- data.frame(rownames(dest$table),round(dest$table$logFC,digits=2),signif(p.adjust(dest$table$PValue,method='BH'),digits=3),signif(dest$table$PValue,digits=3))
#   res.edgeR.sort <- data.frame(rownames(dest.sort$table),round(dest.sort$table$logFC,digits=2),signif(p.adjust(dest.sort$table$PValue,method='BH'),digits=3),signif(dest.sort$table$PValue,digits=3))
#   names(res.edgeR) <- c('GeneID','log2FC','p.adjust','pvalue')
#   names(res.edgeR.sort) <- c('GeneID','log2FC','p.adjust','pvalue')
#   out = list(res.edgeR=res.edgeR, res.edgeR.sort = res.edgeR.sort)
#   out
# }
# count = read.csv(
#   file = 'volcano_plot_example_data.csv',
#   header = TRUE, check.names=FALSE,
#   stringsAsFactors = F
# )
# count = as.matrix(count)
volcano_plot=function(deg){#传入的是有pvalue的那种数据！！！
  deg$group<-ifelse(deg$log2FC>=2&deg$pvalue<=0.05,"Up",
                   ifelse(deg$log2FC<=-2&deg$pvalue<=0.05,
                          "Down","Not sig"))
  p=ggplot(deg,aes(x=log2FC,y=-log10(pvalue)))+
    geom_point(aes(color=group))+
    scale_color_manual(values=c("dodgerblue","gray","firebrick"))
  p
}
# volcano_plot(count)


#trajectory
library(SingleCellExperiment)
library(monocle)
library(scater)
library(ggplot2)
library(igraph)
my.monocle <- function(count, cellLabels){
  colnames(count) <- 1:ncol(count)
  geneNames <- rownames(count)
  rownames(count) <- 1:nrow(count)
  pd <- data.frame(timepoint = cellLabels)
  pd <- new("AnnotatedDataFrame", data=pd)
  fd <- data.frame(gene_short_name = geneNames)
  fd <- new("AnnotatedDataFrame", data=fd)
  
  dCellData <- newCellDataSet(count, phenoData = pd, featureData = fd, expressionFamily = uninormal())
  
  dCellData <- detectGenes(dCellData , min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(dCellData),
                                      num_cells_expressed >= 50))
  
  diff_test_res <- differentialGeneTest(dCellData[expressed_genes,],
                                        fullModelFormulaStr = "~timepoint",
                                        cores = 3)
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  
  dCellData <- setOrderingFilter(dCellData, ordering_genes)
  
  dCellData <- reduceDimension(dCellData, max_components = 2,
                               method = 'DDRTree', norm_method='none')  #auto_param_selection = F
  dCellData <- orderCells(dCellData)
  
  cor.kendall = cor(dCellData@phenoData@data$Pseudotime, as.numeric(dCellData@phenoData@data$timepoint), 
                    method = "kendall", use = "complete.obs")
  
  lpsorder2 = data.frame(sample_name = colnames(count), State= dCellData@phenoData@data$State, 
                         Pseudotime = dCellData@phenoData@data$Pseudotime, rank = rank(dCellData@phenoData@data$Pseudotime))
  
  lpsorder_rank = dplyr::arrange(lpsorder2, rank)
  
  lpsorder_rank$Pseudotime = lpsorder_rank$rank
  
  lpsorder_rank = lpsorder_rank[-4]
  
  lpsorder_rank[1] <- lapply(lpsorder_rank[1], as.character)
  
  subpopulation <- data.frame(cell = colnames(count), sub = as.numeric(cellLabels)-1)
  
  POS <- TSCAN::orderscore(subpopulation, lpsorder_rank)[1]
  out=list(score=data.frame(criteria=c("cor.kendall","POS"),value=c(cor.kendall,POS)),dCellData=dCellData)
  out
}
trascore=function(count){
  cellLabels=factor(colnames(count))
  out=my.monocle(count,cellLabels)
  out
}
# deng <- readRDS("D:\\系统默认\\桌面\\后续实验\\需要插补的数据集\\deng_count_pre.rds")
# count = read.csv(
#   file = "D:\\系统默认\\桌面\\后续实验\\需要插补的数据集\\deng.csv",
#   header = TRUE,
#   row.names = 1, check.names=FALSE
# )
# count = as.matrix(count)
# dengscore=trascore(deng)
# ggplot(dengscore$score,aes(x=criteria,y=value,fill=criteria))+
#     geom_bar(stat = "identity")+
#     theme(legend.position ='none')+xlab('aaa')
# plot_cell_trajectory(dengscore$dCellData, color_by = "timepoint")

