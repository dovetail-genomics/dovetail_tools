#!/usr/bin/env Rscript

library("HiCcompare")
library("argparse")
library("parallel")

parser <- ArgumentParser(description='Run HiCCompare on case and control')
parser$add_argument('--hic1', help='First .mcool file')
parser$add_argument('--hic2', help='Second .mcool file')
parser$add_argument("--res", type="integer", help="Resolution to call differential interactions")
parser$add_argument("--out", help="Output file to write interactions")

args <- parser$parse_args()

hic1 <- args$hic1
hic2 <- args$hic2
res <- args$res

chr.t <- read.table(pipe(paste("cooler dump -t chroms ",paste0(hic1,'::/resolutions/',res))))

getCool <- function(mcool,chr,res){
    cool.f <- paste0(mcool,'::/resolutions/',res)
    region <- paste0(chr)
    cmd <- pipe(paste("cooler dump -t pixels -r",region,"-r2",region,"--header --join",cool.f))
    interaction <- read.table(cmd,header=TRUE,fill=TRUE,stringsAsFactors=FALSE)
}


ans <- do.call(rbind,mclapply(chr.t$V1, function(chr) {
    ##first_hic <- paste(tmpdir,"/first.txt",sep="")
    ##system(paste("straw NONE",hic1,chr,chr,"BP", res, ">", first_hic, sep= " "))
    ##second_hic <-  paste(tmpdir,"/second.txt",sep="")
    ##system(paste("straw NONE",hic2,chr,chr,"BP", res, ">", second_hic, sep= " "))
    ##first_hic_table <- read.table(first_hic, header=FALSE, stringsAsFactors = FALSE)
    ##second_hic_table <- read.table(second_hic, header=FALSE, stringsAsFactors = FALSE)
    first_hic_table <- getCool(hic1,chr,res)
    second_hic_table <- getCool(hic2,chr,res)
    hic_table <- create.hic.table(first_hic_table, second_hic_table, chr=paste("chr",chr,sep=""))
    hic.table <- hic_loess(hic_table, Plot = FALSE, Plot.smooth = FALSE,  parallel=TRUE)
    hic.table <- hic_compare(hic.table, A.min = 15, adjust.dist = TRUE, p.method = 'fdr', Plot = FALSE,  parallel=TRUE)
    return (hic.table)
}, mc.cores = as.integer(system("nproc", intern=TRUE))/2, mc.preschedule = FALSE))


write.table(as.data.frame(ans), file=args$out, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE)

