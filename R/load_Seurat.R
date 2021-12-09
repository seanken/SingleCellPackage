#'Make Seurat object
#'
#'Makes a Seurat object given a count matrix, or does downstream analysis on already existing Seurat object
#'
#' @param seur A seurat object, either this or dat required but not both
#' @param dat A count matrix, gene by cell, either this or seur required
#' @param minGenes Minimium number of genes per cell
#' @param regress List of parameters to regress out, pass c() if nothing to regress out
#' @param species Species (not currently used)
#' @return seur A seurat object
#' @export
dir10X<-function(seur=NULL,dat=NULL,minGenes=500,regress=c("nCount_RNA"),species="human")
{
if(is.null(seur))
{
print(paste("Dims: ",toString(dim(dat))))

print("Make object!")
seur<-CreateSeuratObject(dat,"Seurat",min.features=minGenes)#,normalization.method="LogNormalize",scale.factor=1000000)
}

seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=1000000)

print("Get variable genes!")
seur<-FindVariableFeatures(seur)


print("Regress out!")
if(length(regress)>0)
{
seur<-ScaleData(seur,features=seur@assays$RNA@var.features,vars.to.regress=regress)

}
else{
seur<-ScaleData(seur,features=seur@assays$RNA@var.features)
}
print("Run PCA!")
seur<-RunPCA(seur,npcs=60)

print("Run UMAP and Clustering")
seur=RunUMAP(seur,dims=1:20)
seur=FindNeighbors(seur,dims=1:20)
seur=FindClusters(seur)

print("Return!")

return(seur)

}


