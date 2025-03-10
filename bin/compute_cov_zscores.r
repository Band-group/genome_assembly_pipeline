#!/usr/bin/env Rscript

library(dplyr)
library(argparse)

parser <- ArgumentParser(description = " Plot breadth of coverage from bedtools genomecov output.")

# Define arguments
parser$add_argument("--coverage_file",
                    type = "character",
                    help = "Per base coverage data produced by bedtools genomecov",
                    required = TRUE)

parser$add_argument("--outdir",
                    type = "character",
                    help = "Output directory for zscore regions file",
                    required = TRUE)

# Parse arguments
args <- parser$parse_args()

# Read coverage file into table 
data <- read.table(args$coverage_file)

# Calculate mean and var coverage for each contig
contig_stats <- data %>% 
  group_by(V1) %>% 
  summarise(mean_coverage = mean(V3),
            var_coverage = var(V3))

# Function to calculate stats for windows along each contig
calc_sliding_windows <- function(per_base_coverage_stats, contig_stats, contig_name, window_size = 1000) {
  contig_data <- per_base_coverage_stats %>% 
    filter(V1 == contig_name)
  min_pos <- min(contig_data$V2)
  max_pos <- max(contig_data$V2)
  window_starts <- seq(min_pos, max_pos, by = window_size)
  
  # Calculate mean for each window
  results <- data.frame()
  for (start in window_starts) {
    end <- start + window_size - 1
    window_data <- contig_data %>% filter(V2 >= start, V2 <= end)
    
    if (nrow(window_data) > 0) {
      window_mean <- mean(window_data$V3)
      contig_mean <- contig_stats$mean_coverage[contig_stats$V1 == contig_name]
      contig_sd <- sqrt(contig_stats$var_coverage[contig_stats$V1 == contig_name])
      z_score <- (window_mean - contig_mean) / contig_sd
      
      new_row <- data.frame(
        contig = contig_name,
        window_start = start,
        window_end = end,
        window_mean = window_mean,
        z_score = z_score
      )
      results <- rbind(results, new_row)
    }
  }
  return(results)
}

contig_names <- unique(data$V1)
all_windows <- data.frame()
for (name in contig_names) {
  cat("Processing contig:", contig, "\n")
  contig_windows <- calc_sliding_windows(per_base_coverage_stats = data, contig_stats = contig_stats, contig_name = name)
  all_windows <- rbind(all_windows, contig_windows)
}

# Filter windows with absolute z-score greater than threshold
options(scipen = 999) # avoid writing bed coordinates in standard form
zscore_threshold <- args$zscore_threshold
high_zscore_windows <- all_windows %>%
  filter(abs(z_score) > zscore_threshold) %>%
  arrange(contig, window_start)

# Merge contiguous windows into regions
merged_regions <- data.frame(
  contig = character(),
  start = numeric(),
  end = numeric(),
  mean_zscore = numeric(),
  stringsAsFactors = FALSE
)

# Process each contig separately
for (curr_contig in unique(high_zscore_windows$contig)) {
  contig_windows <- high_zscore_windows %>% filter(contig == curr_contig) %>% arrange(window_start)
  
  if (nrow(contig_windows) > 0) {
    # Initialize with first window
    current_start <- contig_windows$window_start[1]
    current_end <- contig_windows$window_end[1]
    current_zscores <- contig_windows$z_score[1]
    
    # Process remaining windows
    for (i in 2:nrow(contig_windows)) {
      # Check if this window is contiguous with current region
      if (contig_windows$window_start[i] <= current_end + 1) {
        # Extend current region
        current_end <- contig_windows$window_end[i]
        current_zscores <- c(current_zscores, contig_windows$z_score[i])
      } else {
        # Save current region and start a new one
        merged_regions <- rbind(merged_regions, data.frame(
          contig = curr_contig,
          start = current_start,
          end = current_end,
          mean_zscore = mean(current_zscores),
          stringsAsFactors = FALSE
        ))
        
        # Start new region
        current_start <- contig_windows$window_start[i]
        current_end <- contig_windows$window_end[i]
        current_zscores <- contig_windows$z_score[i]
      }
    }
    
    # Add the last region
    merged_regions <- rbind(merged_regions, data.frame(
      contig = curr_contig,
      start = current_start,
      end = current_end,
      mean_zscore = mean(current_zscores),
      stringsAsFactors = FALSE
    ))
  }
}

# Convert to BED format (0-based start coordinate)
bed_regions <- merged_regions %>%
  mutate(start = start - 1) %>%  # Convert to 0-based coordinates for BED format
  select(contig, start, end, mean_zscore)

# Write BED file
output_file <- file.path(args$outdir, "high_zscore_regions.bed")
write.table(bed_regions, file = output_file, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat("High z-score regions written to:", output_file, "\n")