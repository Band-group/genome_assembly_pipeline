library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(argparse)

# Create parser object
parser <- ArgumentParser(description = " Filter BAM files based on regions")

# Define arguments
parser$add_argument("--outdir",
                    type = "character",
                    help = "Output directory for filtered BAMs",
                    required = TRUE)

parser$add_argument("--regions_bed",
                    type = "character",
                    help = "BED file containing regions to filter on",
                    required = TRUE)

parser$add_argument("--bams",
                    type = "character",
                    nargs = "+",
                    help = "Input BAM files to filter",
                    required = TRUE)

# Parse arguments
args <- parser$parse_args()

# Rest of script using args$outdir, args$regions_bed, args$bams
regions <- import(args$regions_bed, format = "BED")
names <- file.path(args$outdir, sub("\\..*$", "_regionRestricted", basename(args$bams)))

for (i in seq_along(args$bams)) {
    temp_bam <- tempfile(pattern = "temp", tmpdir = args$outdir, fileext = ".bam")
    reads_gal <- readGAlignments(args$bams[i])
    filtered_reads <- subsetByOverlaps(reads_gal, regions, type = "within")
    export(filtered_reads, BamFile(temp_bam))
    sortBam(temp_bam, destination = names[i])
    indexBam(paste0(names[i], ".bam"))
    unlink(c(temp_bam, paste0(temp_bam, ".bai")))
}