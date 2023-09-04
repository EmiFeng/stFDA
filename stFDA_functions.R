
#Function: merge spatial Seurat objects.
#obj_list: list of st seurat objects. example: obj_list = lsit(rds1,rds2,rds3,rds4)
#sample_list: sample names. example: sample_list = c('S1','S2','S3','S4')
#group_list: group names. example: group_list = c('G1','G1','G2','G2')
#strds <- stFDA_mergeST(lsit(rds1,rds2,rds3,rds4), sample_list = c('S1','S2','S3','S4'), group_list = c('G1','G1','G2','G2'))
#return merged spatial seurat object, sample was stored in strds$sample, group was stored in strds$group
stFDA_mergeST <- function(obj_list, sample_list, group_list){
  num <- length(obj_list)
  merge_obj <- readRDS(obj_list[[1]])
  merge_obj <-  RenameCells(merge_obj, add.cell.id = sample_list[1])
  names(merge_obj@images) <- sample_list[1]
  merge_obj$sample <- sample_list[1]
  merge_obj$group <- group_list[1]
  
  for (i in 2:num){
    obj <- readRDS(obj_list[[i]])
    obj <-  RenameCells(obj, add.cell.id = sample_list[i])
    names(obj@images) <- sample_list[i]
    obj$sample <- sample_list[i]
    obj$group <- group_list[i]
    merge_obj <- merge(merge_obj, obj)
  }
  merge_obj$sample <- factor(merge_obj$sample,level = sample_list)
  merge_obj$group <- factor(merge_obj$group,level = unique(group_list))
return(merge_obj)
}

#strds <- stFDA_mergeST(obj_list, sample_list = sample, group_list = group)

#Function: get cell type signature genesets from scrna data.
#scrds: single-cell seurat object with cell type annotation
#number: maximum number of genes selected for a signature
#group.by: whether to change active.ident, default not to change; if change, give the colname stored in meta.data, for example, 'new_ident'
#downsample: whether to downsample cells for idents, default not to; suggest to use when cell number if too large, for exsample, downsample = 2000
#return a list of signatures, can be used for function stFDA_Scoring()
#sc_signature <- stFDA_scSignatures(scrds, downsample = 2000, number = 50)
stFDA_scSignatures <- function(scrds, number= 50 , group.by = 'ident', downsample = NULL){
  if (group.by != 'ident'){Idents(scrds) <- group.by}
  print(levels(scrds))
  if (!is.null(downsample)){scrds <- subset(scrds, downsample = as.numeric(downsample))}

  geneset1 <- list()
  
  for (i in levels(scrds)){
      diff <- FindMarkers(scrds,ident.1 = i)
      diff <- diff[diff$p_val < 0.05 & diff$avg_log2FC > 0,]
      diff <- diff[order(diff$avg_log2FC,decreasing = T),]
      mt.genes <- grep('^MT-', rownames(diff), value = T, ignore.case = T)  #detect mt genes
      rib.genes <- c(grep('^RPL', rownames(diff), value = T, ignore.case = T),grep('^RPS', rownames(diff), value = T, ignore.case = T))   #detect ribosomal genes
      diff <- diff[setdiff(rownames(diff),c(mt.genes, rib.genes)),] #remove mt and ribosomal genes
      nr <- nrow(diff)
      genes <- rownames(diff)[1:min(nr, number)]
	    geneset1[[i ]] = genes

  }
#print(geneset1)
return(geneset1)
}

#Function: spatial feature plot
stFDA_featureplot <- function(obj, matrix, prefix, outdir,image_alpha = 1){
dir.create(outdir,showWarnings = FALSE)
celltype <- colnames(matrix)
ratio = length(obj@images)

for (ct in celltype){
name <- paste(ct, 'Score')
obj@meta.data[[name]] <- matrix[,ct]

#minscore <- min(rds@meta.data[[name]])
#rds@meta.data[[name]] <- rds@meta.data[[name]] - minscore
#maxscore <- max(rds@meta.data[[name]])
#rds@meta.data[[name]] <- rds@meta.data[[name]]/maxscore
p <- SpatialPlot(obj, features = name,image.alpha = image_alpha,stroke = 0.1,pt.size.factor = 1.3, crop = F)
pdf(paste0(outdir,'/',prefix,'.',ct,'.SpatialFeature.pdf'),width = 6*ratio, height = 6)
print(p);dev.off()
png(paste0(outdir,'/',prefix,'.',ct,'.SpatialFeature.png'),units = 'in',res = 300, width = 6*ratio, height = 6)
print(p);dev.off()

}
}

#Function: integrated feature scoring
#obj: spatial seurat object, if you have multiple samples, you can use the result from stFDA_mergeST.
#features: functional genesets, can be feature list (return from stFDA_scSignatures()) or gmt file. 
#method: choose one scoring method from "AUCell", "UCell", "singscore", "ssgsea", default = "ssgsea".
#assay: which assay to use, default = 'Spatial', check Seurat object for more information.
#assay: which data slot to use, default = 'data', check Seurat object for more information.
#custom: Default TRUE. Set it to TRUE when input own genesets; set it to FALSE when use msigdb genesets.
#msigdb: Default FALSE. Set it to TRUE when custom = FALSE.
#category: Choose a msigdb category for analysis, default='H' . Use msigdbr::msigdbr_collections to view all available collections gene sets.
#species: Default Homo sapiens. Use msigdbr::msigdbr_show_species() to view all available species. The parameter works if msigdb is True.
#outdir: output directory path.
#prefix: prefix of output files.


stFDA_Scoring <- function(obj, features,gmt = FALSE, method = 'ssgsea',custom = TRUE, assay = 'Spatial',slot = 'data',msigdb = FALSE, category = 'H', species = 'Homo sapiens', outdir = './',prefix = 'spatial'){
if (custom){
  if (gmt){
	  gsets = read.gmt(features)
	  geneset1 = list()
	  pathway = unique(gsets$term)
      for(i in pathway){
        gset = subset(gsets, term == i)
	      rawgene = as.character(gset$gene)
	      geneset1[[i ]] = rawgene
    }
  }else {geneset1 = features}

  rds <- irGSEA.score(object = obj, assay = assay,
                             slot = slot, seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = custom , geneset = geneset1, msigdb = msigdb,
                             subcategory = NULL, geneid = "symbol",
                             method = method,category = category, species = species,
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                             kcdf = 'Gaussian')
}else{
  rds <- irGSEA.score(object = obj, assay = assay,
                             slot = slot, seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = FALSE , msigdb = TRUE,
                             subcategory = NULL, geneid = "symbol",
                             method = method,category = category, species = species,
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                             kcdf = 'Gaussian')
}

matrix<-rds@assays[[method]]@counts
matrix <- t(as.matrix(matrix))
stFDA_featureplot(obj, matrix, prefix , outdir, image_alpha = 1)
return(matrix)
}



#Plot assigned feature plots.
stFDA_assignedPlot <- function(assign, obj, prefix, outdir,image_alpha = 1){
celltypes <- colnames(assign)
ratio = length(obj@images)
assign = as.data.frame(assign)
for (ct in celltypes){
subassign <- assign[assign[[ct]] == 1,]
cells <- rownames(subassign)
p <- SpatialPlot(obj, cells.highlight = cells,image.alpha = image_alpha,stroke = 0.1,pt.size.factor = 1.3, crop = F, cols.highlight = c("#DE2D26", "lightgray"))
pdf(paste0(outdir,'/',prefix,'.',ct,'.assigned.pdf'),width = 6*ratio, height = 6)
print(p);dev.off()
png(paste0(outdir,'/',prefix,'.',ct,'.assigned.png'),units = 'in',res = 300, width = 6*ratio, height = 6)
print(p);dev.off()

}

}

#Function: Based on the table of feature scoring from stFDA_Scoring(), stFDA_threshold() and stFDA_enrichedSpots() could detect spatially enriched spots. Return assigned table as well as the fraction of enriched spots of every sample.
#data: the table of feature scoring from stFDA_Scoring()
#obj: spatial seurat object, if you have multiple samples, you can use the result from stFDA_mergeST.
#outlier: to define a high-threshold as the maximum score after excluding top X percent outliers with highest scores. default = 0.01
#percent: to determine a low-threshold cross which was regarded as a enriched spot. default = 0.7
#outdir: output directory path.
#prefix: prefix of output files.
stFDA_threshold <- function(data = data, obj, outlier = 0.01, percent = 0.7,outdir = '.',prefix = 'spatial'){
dir.create(outdir)
dir.create(paste0(outdir,'/threshold/'))
assign <- data

geneset <- colnames(data)

for (g in geneset){

sub <- data.frame(Score = data[,g])
rownames(sub) <- rownames(data)

thred = percent*quantile(sub$Score,1-outlier)
#thred2 = percent*max(ref$TCells)
p <- ggplot(sub,aes(Score))+geom_density()+geom_vline(xintercept = thred,linetype = 'dotted')+theme_bw()+ggtitle(paste(g,'threshold:',thred))
pdf(paste0(outdir,'/threshold/',g,'.density_threshold.pdf'));print(p);dev.off()
assign[,g] <- 0
assign[sub$Score > thred,g] <- 1
}

write.table(assign, paste0(outdir,'/assign_table.xls'),sep = '\t', quote = F, row.names = T, col.names = NA)
stFDA_assignedPlot(assign = assign, obj = obj, prefix , outdir)
return(as.data.frame(assign))
}

#ssgsea_result <- stFDA_threshold(data = ssgsea_table, outlier = 0.01, percent = 0.7, outdir = '/SGRNJ03/PiplineTest01/Pipline_test/fengshijing/Scoring_exploreThr/ssgsea/')
#auc_result <- stFDA_threshold(data = auc_table, outlier = 0.01, percent = 0.7, outdir = '/SGRNJ03/PiplineTest01/Pipline_test/fengshijing/Scoring_exploreThr/aucell/')


stFDA_enrichedSpots <- function(assign, strds, assay = 'Spatial', prefix = 'spatial',outdir = '.'){

DefaultAssay(strds) <- assay
cells.use <- intersect(rownames(assign),Cells(strds))
assign <- assign[cells.use,]
strds <- subset(strds,cells = cells.use)
celltypes <- colnames(assign)
sam_group <- unique(strds@meta.data[,c('sample','group')])

tbdata <- data.frame(sample = sam_group$sample, group = sam_group$group,percent = 0,celltype = colnames(assign)[1])
for (i in 2:ncol(assign)){
tmp <- data.frame(sample = sam_group$sample, group = sam_group$group,percent = 0,celltype = colnames(assign)[i])
tbdata <- rbind(tbdata,tmp)
}


for (sp in unique(sam_group$sample)){
cells.sub  <- Cells(subset(strds, sample == sp))
sub <- assign[cells.sub,]
for (celltype in celltypes){tbdata[tbdata$sample == sp & tbdata$celltype == celltype,'percent'] <- nrow(sub[sub[[celltype]] == 1,])/nrow(sub)}
}

write.table(tbdata, paste0(outdir,'/',prefix, '.enrichedSpots_pct.xls'), sep = '\t', quote = F, row.names = F)
return(tbdata)

}

#ssgsea_enriched <-  stFDA_enrichedSpots(assign = ssgsea_result, strds = strds, outdir = '/SGRNJ03/PiplineTest01/Pipline_test/fengshijing/Scoring_exploreThr/ssgsea/', prefix = 'ssgsea')


stFDA_diffplot <- function(table, colors = NULL, prefix = 'spatial',outdir = '.', stat_test = FALSE, test = 'wilcox.test'){
table$sample <- factor(table$sample, level = unique(table$sample))
table$group <- factor(table$group, level = unique(table$group))
p <- ggplot(table, aes(x = celltype, y = percent, fill = sample))+ geom_bar(position = "dodge", width = 0.6, stat = "identity")+
      theme_bw() + xlab('')+ylab('Fraction of Enriched Spots')+theme(axis.text.x = element_text(angle = 45, hjust = 1))

if (!is.null(colors)){p <- p +scale_fill_manual(values = colors)}

pdf(paste0(outdir,'/',prefix,'.sample_barplot.pdf'),width = 6*16/9, height = 6)
print(p);dev.off()
png(paste0(outdir,'/',prefix,'.sample_barplot.png'),units = 'in',res = 300, width = 6*16/9, height = 6)
print(p);dev.off()

p <- ggplot(table, aes(x = celltype, y = percent, fill = group, color = group))+
    geom_boxplot(outlier.size = 0, alpha = 0.5)+theme_bw() + xlab('')+ylab('Fraction of Enriched Spots')+theme(axis.text.x = element_text(angle = 45, hjust = 1))

if (!is.null(colors)){p <- p +scale_fill_manual(values = colors)+scale_color_manual(values = colors) }

pdf(paste0(outdir,'/',prefix,'.group_boxplot.pdf'),width = 6*16/9, height = 6)
print(p);dev.off()
png(paste0(outdir,'/',prefix,'.group_boxplot.png'),units = 'in',res = 300, width = 6*16/9, height = 6)
print(p);dev.off()

if (stat_test){
#library(ggpubr, lib.loc= '/SGRNJ01/Public/Software/conda_env/CeleScanpz/lib/R/library/' )
stat_test<- compare_means(percent~group,data = table,group.by = "celltype",method = test)
write.table(stat_test, paste0(outdir,'/',prefix,'.enrichedSpots_sig.xls'), sep = '\t', quote = F, row.names = F)
}

}



