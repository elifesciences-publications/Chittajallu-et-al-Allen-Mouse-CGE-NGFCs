#
# Collection of functions to convert data to cellranger-like format
#

library(Matrix)
library(hdf5r)

#
# Input: TSV or CSV file from GEO, .tome (HDF5) file from Allen institute
# Output: cellranger < v3 files: matrix.mtx, genes.tsv, barcodes.tsv
#

# function to write sparse matrix to cellranger-like files
# 
write_10x <- function(sparse.mat, out.dir, gene.names=NULL, cell.names=NULL){

    if(!dir.exists(out.dir)){
        dir.create(out.dir, showWarnings=FALSE, recursive=TRUE)
    }

    writeMM(obj=sparse.mat, file=file.path(out.dir, 'matrix.mtx'))

    # gene file, barcode file
    genes <- file.path(out.dir, 'genes.tsv')
    cells <- file.path(out.dir, 'barcodes.tsv')

    #
    # if sparse matrix has row/colnames output as genes.tsv, barcodes.tsv
    # otherwise, output specified gene.names, cell.names
    #
    if(!is.null(rownames(sparse.mat))){
        write(x = rownames(sparse.mat), file=genes)
    } else {
        if(is.null(gene.names)){
            stop('gene names not specified')
        }
        write(x = gene.names, file=genes)
    }

    if(!is.null(colnames(sparse.mat))){
        write(x = colnames(sparse.mat), file=cells)
    } else {
        if(is.null(cell.names)){
            stop('cell names not specified')
        }
        write(x = cell.names, file=cells)
    }
}

# function to convert HDF5 file from Allen institute (.tome)
# to sparse matrix
#
# tome = .tome file (HDF5 format) from Allen institute
# out.dir = output directory, will create if not existing
#
tome_to_10x <- function(tome, out.dir){
    h5file <- H5File$new(tome, mode='r')

    # 
    # Relevant data is in
    # counts: data/t_exon, gene x cell matrix
    # gene names: gene_names
    # sample names: sample_names
    #
    counts <- h5file[['data/t_exon']]
    gene_names <- h5file[['gene_names']]
    sample_names <- h5file[['sample_names']]

    # create sparse matrix
    sparse.mat <- sparseMatrix(i = counts[['i']][],
                               p = counts[['p']][],
                               x = counts[['x']][],
                               dims = counts[['dims']][],
                               dimnames = list(gene_names[], sample_names[]),
                               index1 = FALSE)

    write_10x(sparse.mat, out.dir)
}

# function to convert counts TSV/CSV to sparse matrix
#
# count = TSV or CSV file name
# out.dir = output directory, will create if not existing
# sep = '\t' for TSV, ',' for CSV
#
count_to_10x <- function(count, out.dir, sep='\t'){
    
    mat <- read.table(count, sep=sep, header=TRUE, row.names=1)

    # sparse matrix
    sparse.mat <- Matrix(as.matrix(mat), sparse=TRUE)

    write_10x(sparse.mat, out.dir)
}


#
# this function parses the TSV/CSV from fluidigm C1 (GSE109796)
# and creates cellranger-like files: matrix.mtx, genes.tsv, barcodes.tsv
#
# 1. strip trailing bit from column names
# 2. parse cell barcode
# 3. parse gene names (expected format: ID|symbol)
#
fluidigm_to_10x <- function(count, out.dir,
                            sample.pattern,
                            strip.trailing=NULL,
                            gene.sep='\\|',
                            sep='\t'){
    library(Matrix)

    mat <- read.table(count, sep=sep, header=TRUE, row.names=1)

    # sparse matrix
    sparse.mat <- Matrix(as.matrix(mat), sparse=TRUE)

    # strip trailing string: 
    # e.g. strip.trailing = '_L.+', matches '_L008_R1_all'
    if(!is.null(strip.trailing)){
        new.colnames <- sub(strip.trailing,'',colnames(sparse.mat))
    } else {
        new.colnames <- colnames(sparse.mat)
    }

    # get sample id
    # e.g. sample.pattern='C\\d.+[A-Z]\\d+' matches C1.101.A10
    p <- regexpr(sample.pattern, new.colnames)
    sample.ids <- substr(new.colnames,
                         p, attr(p, 'match.length'))
    
    # get cell umis
    cell.umis <- substr(new.colnames,
                        p + attr(p, 'match.length') + 1,
                        nchar(new.colnames))
    # remove any non-alphabetic characters from UMIs
    #cell.umis <- sub('[^A-Z]+','',cell.umis)

    colnames(sparse.mat) <- paste0(sample.ids, '_', cell.umis)

    # replace gene id and symbol separators with tab
    rownames(sparse.mat) <- sub(gene.sep,'\t',rownames(sparse.mat))

    write_10x(sparse.mat, out.dir)
}
