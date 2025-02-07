#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(argparse)

# Create parser object
parser <- ArgumentParser(description = " Plot breadth of coverage from bedtools genomecov output.")

# Define arguments
parser$add_argument("--outdir",
                    type = "character",
                    help = "Output directory for coverage plot",
                    required = TRUE)

parser$add_argument("--sample_names",
                    type = "character",
                    nargs = "+",
                    help = "Sample names corresponding to each coverage file",
                    required = TRUE)

parser$add_argument("--coverage_files",
                    type = "character",
                    nargs = "+",
                    help = "Coverage data produced by bedtools genomecov",
                    required = TRUE)

parser$add_argument("--regions_bed",
                    type = "character",
                    help = "Optional bed file specifying regions that should be included in the coverage calculation",
                    required = FALSE)

parser$add_argument("--ref_length",
                    type = "integer",
                    help = "Length of the reference sequence. Required if --regions_bed is provided",
                    required = FALSE)

# Parse arguments
args <- parser$parse_args()

# Function to process coverage file
process_coverage <- function(file, sample_name, regions_bed = NULL, ref_length = NULL) {
    cov <- read.table(file)
    colnames(cov) <- c("chromosome", "depth", "bases", "size", "fraction")
    if (!is.null(regions_bed)) {
        regions <- read.table(regions_bed)
        diff <- ref_length - max(cumsum(regions$V3 - regions$V2))
        df_coverage <- cov %>%
            filter(chromosome == "genome") %>%
            mutate(bases = ifelse(row_number() == 1, bases - diff, bases)) %>%
            mutate(fraction = bases / size) %>%
            reframe(depth = depth + 1, fraction = 1 - cumsum(fraction)) %>%
            mutate(fraction = pmax(fraction, 0), sample = sample_name)
    } else {
        df_coverage <- cov %>%
            filter(chromosome == "genome") %>%
            reframe(depth = depth + 1, fraction = 1 - cumsum(fraction)) %>%
            mutate(fraction = pmax(fraction, 0), sample = sample_name)
    }
    return(df_coverage)
}

# Clean sample names
unclean_sample_names <- args$sample_names
clean_sample_names <- sub("\\[|\\]", "", unclean_sample_names)

# Process all coverage files
coverage_list <- mapply(process_coverage,
                        file = args$coverage_files,
                        sample_name = clean_sample_names,
                        regions_bed = args$regions_bed,
                        ref_length = args$ref_length,
                        SIMPLIFY = FALSE)
all_coverage <- bind_rows(coverage_list)

# Create single plot with multiple lines
p <- ggplot(all_coverage, aes(x = depth, y = fraction, colour = sample)) +
    geom_line(linewidth = 1) +
    scale_x_log10(
        limits = c(1, 525),
        breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500),
        expand = c(0, 0),
        labels = c("1", "2", "5", "10", "20", "50", "100", "200", "500")
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
        expand = c(0, 0),
        labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")
    ) +
    labs(
        x = "Minimum Read Depth (log scale)",
        y = "Proportion of Genome Covered"
    ) +
    theme_bw() +
    theme(
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        panel.border = element_rect(linewidth = 1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

ggsave(file.path(args$outdir, "coverage_plot.pdf"), p, width = 12, height = 9, units = "in")
