library(Seurat)



##
##Gets shared variable genes between all batches
##
#' Shared Variable Genes
#'
#' Calculaates variable genes shared between batches
#'
#' @param seur The Seurat object to run Variable Genes on
#' @param x.low.cutoff Parameter for variable genes
#' @param x.high.cutoff Parameter for variable genes
#' @param minNum Minimium number of batches/samples a variable gene must be found in to be used
#' @param batch Column in meta.data containing batch or sample information
#' @param minCells The min number of cells for a batch to be used
#' @return A seurat object with new variable genes and active.ident set to batch
#' @export
SharedVariable<-function(seur,x.low.cutoff=1,x.high.cutoff=5,minNum=3,batch="orig.ident",minCells=10)
{
batchs<-unique(seur@meta.data[,batch])

seur<-SetIdent(seur,value=batch)

vars<-c()

for(bat in batchs)
{
print(bat)
temp<-SubsetData(seur,cells=WhichCells(seur,idents=bat))
if(length(temp@active.ident)>minCells)
{
temp<-FindVariableFeatures(temp,selection.method="mean.var.plot",mean.cutoff=c(x.low.cutoff,x.high.cutoff))
vars<-c(vars,temp@assays$RNA@var.features)
}
print(" ")
}

print("Combine!")
vars<-table(vars)

genes<-names(vars)[vars>minNum]

print(length(genes))

seur@assays$RNA@var.features<-genes

return(seur)

}


