#!/usr/bin/env Rscript

# only applicable when running on rackham cluster
# module load R/4.0.0
# module load R_packages/4.0.0
## .libPaths("/sw/apps/R_packages/4.0.0/rackham")

bam_path <- "/crex/course_data/ngsintro/rnaseq/main_full/3_mapping/"

# For running as command line script and parsing options from command line
args = commandArgs(trailingOnly=TRUE)
# test if correct number of arguments are given as input: if not, return an error
if (length(args)<3) {
    stop("Please supply chromosome number as well as start and stop position for the plot\n Example: Rscript gene.R 2 234243 238555",
         call.=FALSE)
}

library(Gviz)
library(EnsDb.Mmusculus.v79)
edb <- EnsDb.Mmusculus.v79
gen <- "mm10"
chr <- args[1] # parse the first argument after the r-script to chr
start <- args[2] # parse the second argument after the r-script to start
stop <- args[3] # parse the second argument after the r-script to stop

cat("Reading BAM files ...\n")

# Collect bam filenames from bam folder
bams <- list.files(bam_path, pattern = "*.bam$",full.names = TRUE, recursive = TRUE)

# Allow for using different chromosome names than ucsc
options(ucscChromosomeNames=FALSE)

# create track see Gviz manual for details
gat <- GenomeAxisTrack()
gr <- getGeneRegionTrackForGviz(edb, chromosome = args[1], start = start, end = stop)
genome(gr) <- "mm10"

tracklist <- lapply(seq_along(bams), function(x) AlignmentsTrack(bams[[x]],chromosome = chr,genome = gen,name = basename(bams[x])))
# Merge all tracks
toplot <- append(list(gat, GeneRegionTrack(gr)),tracklist)

width = 14
height = 14

# Plot to pdf
cat("Exporting plots ...\n")
pdf("Coverage.pdf", width = width, height = height)  
plotTracks(toplot, type = c("coverage"))
dev.off()
pdf("Sashimi.pdf", width = width, height = height)  
plotTracks(toplot, type = c("sashimi"))
dev.off()
