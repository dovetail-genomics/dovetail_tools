suppressPackageStartupMessages(library(Sushi))
suppressPackageStartupMessages(library(argparse))




parser <- ArgumentParser(description='Plotting HiChIP Interactions')
parser$add_argument('--bg', help='Coverage BedGraph File')
parser$add_argument('--cool', help='Cooler matrix dump text file')
parser$add_argument('--chrom', help="Chromosome name")
parser$add_argument("--start", type="integer", help="Start of the region")
parser$add_argument("--end", type="integer", help="End of the region")
parser$add_argument('--prefix', help='Output Prefix')

args <- parser$parse_args(c("--bg", "--cool", "--prefix"))


cov <- read.table(args$bg)
arc <- read.table(args$cool)

arc <- arc[which(arc$V1 == arc$V4),]
arc$V8 <- abs(arc$V2 - arc$V6)  

#this is arbitrary filtering
arc <- arc[which(arc$V8 > 20000 & arc$V7 > 1),]
chrom = args$chrom
chromstart = args$start
chromend = args$end

plotBedgraph(cov,chrom,chromstart,chromend)
labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)

plotBedpe(arc,chrom,chromstart,chromend,heights = arc$V8,plottype="loops", flip=TRUE)
labelgenome(chrom, chromstart,chromend,side=3, n=3,scale="Mb")
axis(side=2,las=2,tcl=.2) 
mtext("distance",side=2,line=1.75,cex=.75,font=2)

pdfname <- paste(args$prefix,".hichip.cov.arcs.pdf")
makepdf = TRUE
if(makepdf==TRUE)
{
	pdf(pdfname , height=10, width=12)
}


