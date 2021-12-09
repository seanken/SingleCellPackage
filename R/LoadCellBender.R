#' Load CellBender
#'
#' Loads the output of CellBender
#'
#' @param filename File output by cellbender
#' @return sparse.mat A sparse Matrix with expression counts
#' @export
loadCellBender=function(filename)
{
infile <- hdf5r::H5File$new(filename = filename, mode = "r")
genome="matrix"
counts <- infile[[paste0(genome, "/data")]]
indices <- infile[[paste0(genome, "/indices")]]
indptr <- infile[[paste0(genome, "/indptr")]]

barcodes <- infile[[paste0(genome, "/barcodes")]]
gene_names=infile[[paste0(genome, "/features/name")]][]
gen_id=infile[[paste0(genome, "/features/id")]][]
gene_names[duplicated(gene_names)]=gen_id[duplicated(gene_names)]



sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[],x = as.numeric(x = counts[]), giveCsparse = FALSE,dims=c(length(gen_id),length(barcodes[])))

colnames(x = sparse.mat) <- barcodes[]
rownames(x = sparse.mat) <- gene_names
return(sparse.mat)

}



#' Load CellBender For Multiple Files
#'
#' Loads the output of CellBender for multiple files
#'
#' @param filenames List of files output by cellbender
#' @return dat A sparse Matrix with expression counts
#' @export
loadCellBender_multi=function(filenames)
{
nams=names(filenames)

out=lapply(filenames,function(x){dat=loadCellBender(x)})
for(i in 1:length(out)){colnames(out[[i]])=sub("^",paste(nams[i],"_",sep=""),colnames(out[[i]]))}
lapply(out,function(x){print(dim(x))})
dat=do.call(cbind,out)
return(dat)
}


