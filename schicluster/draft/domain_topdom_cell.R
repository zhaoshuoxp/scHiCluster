library(rhdf5)
library(Matrix)
library(data.table)

# === 1. 动态寻找 TopDom.R ===
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", grep(file.arg.name, initial.options, value = TRUE))
script.dir <- dirname(script.name)
# 尝试向上寻找
topdom.path <- file.path(script.dir, "../domain/TopDom.R")
if (!file.exists(topdom.path)) topdom.path <- file.path(script.dir, "TopDom.R")

if (file.exists(topdom.path)) {
    source(topdom.path)
} else {
    stop("CRITICAL ERROR: Could not find TopDom.R")
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) stop("Not enough arguments.")

cell = args[1]
if (substr(args[2], 1, 3) == 'chr'){chrom = args[2]}else{chrom = paste('chr', args[2], sep='')}
mode = args[3]
ws = as.integer(args[4])
indir = args[5]
outdir = args[6]

fmat = paste(indir, chrom, '/', cell, '_', chrom, '_', mode, '.hdf5', sep = '')
fbin = paste(indir, 'bins/', chrom, '.bed', sep='')

fout = paste(outdir, chrom, '/', cell, '_', chrom, '_', mode, '.w', ws, '.domain.bed', sep = '')

tryCatch({
    if (!file.exists(fmat)) stop(paste("Matrix missing:", fmat))
    if (!file.exists(fbin)) stop(paste("Bins missing:", fbin))

    bins = fread(fbin)
    bins = data.frame(bins)
    colnames(bins) <- c("chr", "from.coord", "to.coord")
    n_bins = nrow(bins)

    # 读取 HDF5
    indices = as.numeric(h5read(fmat, 'Matrix/indices'))
    indptr = as.numeric(h5read(fmat, 'Matrix/indptr'))
    count = as.numeric(h5read(fmat, 'Matrix/data'))
    
    if (length(count) == 0) stop("Matrix data is empty")

    # === [核心修复] ===
    # 1. 使用 i (行索引) 而不是 j (列索引)
    # 2. 这里的逻辑是：Python CSR (row_ptr, col_ind) -> R CSC (col_ptr, row_ind) = Transpose
    #    因为 Hi-C 矩阵是对称的，所以 Transpose 后数据依然正确
    matrix.data = as.matrix(sparseMatrix(i=indices+1, p=indptr, x=count, dims=c(n_bins, n_bins)))
    
    # === [参数优化] ===
    # statFilter=FALSE: 关闭 P-value 过滤，这对稀疏的单细胞数据至关重要
    tad = TopDom(matrix.data, bins, window.size=ws, statFilter=FALSE)
    
    # 确保目录存在并写入
    dir.create(dirname(fout), showWarnings = FALSE, recursive = TRUE)
    write.table(tad$bed, file=fout, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    
    # 打印成功日志 (会被 Python 捕获)
    cat(paste("Success:", cell, "Domains:", sum(tad$bed$tag=='domain'), "\n"))

}, error = function(e) {
    # 打印错误信息
    message(paste("Error processing cell", cell, ":", e$message))
})
