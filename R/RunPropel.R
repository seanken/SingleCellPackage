##example form: ~0+condition+batch
##example contrasts: conditionko-conditionwt
#' Run Propel
#'
#' Run Propeller to test for chances in cluster abundance
#' @param seur A seurat object to test
#' @param form A formula to use (for example ~0+condition+batch)
#' @param contrasts The contrast to use, for example conditionko-conditionwt
#' @param samples The column in meta.data with sample of origin information
#' @param clusters The column in meta.data with cluster/cell type information
#' @return out Results of differential analysis
#' @export
runPropel<-function(seur,form,contrasts,samples,clusters)
{
print("Get Prop")
props <- getTransformedProps(seur@meta.data[,clusters],seur@meta.data[,samples], transform="asin")
print("Get metadata")
tokeep=trimws(strsplit(as.character(form)[[2]],"+",fixed=T)[[1]])
tokeep=intersect(tokeep,colnames(seur@meta.data))
print(tokeep)
meta=seur@meta.data[,c(samples,tokeep)]
print(head(meta))
tab<-dplyr::distinct(seur@meta.data[,c(samples,tokeep)])
rownames(tab)=tab[,samples]
tab=tab[colnames(props$TransformedProps),]
print("Set up")
design <- model.matrix(form,tab)
print(head(design))
print(contrasts)
mycontr <- makeContrasts(contrasts=contrasts,levels=design)
print("Run!")
out=propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE,sort=TRUE)
return(out)
}
