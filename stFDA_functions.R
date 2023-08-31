#Function: merge spatial Seurat objects.
#obj_list: list of st seurat objects. example: obj_list = lsit(rds1,rds2,rds3,rds4)
#sample_list: sample names. example: sample_list = c('S1','S2','S3','S4')
#group_list: group names. example: group_list = c('G1','G1','G2','G2')
#strds <- stFDA_mergeST(lsit(rds1,rds2,rds3,rds4), sample_list = c('S1','S2','S3','S4'), group_list = c('G1','G1','G2','G2'))
#return merged spatial seurat object, sample was stored in strds$sample, group was stored in strds$group
stFDA_mergeST <- function(obj_list, sample_list, group_list){
  num <- length(obj_list)
  merge_obj <- obj_list[[1]]
  merge_obj <-  RenameCells(merge_obj, add.cell.id = sample_list[1])
  names(merge_obj@images) <- sample_list[1]
  merge_obj$sample <- sample_list[1]
  merge_obj$group <- group_list[1]
  
  for (i in 2:num){
    obj <- obj_list[[i]]
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
