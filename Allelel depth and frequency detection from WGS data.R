if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")

library(VariantAnnotation)

# Function to read VCF and extract allele depth information along with chromosome and sample information
extract_allele_depth_chrom_and_samples <- function(vcf_path) {
  vcf <- readVcf(vcf_path, "hg19")
  ad <- geno(vcf)$AD
  chrom <- seqnames(rowRanges(vcf))
  samples <- colnames(vcf)
  return(list(ad = ad, chrom = chrom, samples = samples))
}

# Function to calculate allele depth ratios for each sample and chromosome
calculate_ad_ratios_per_sample_chrom <- function(ad, chrom, samples) {
  ad_ratios_list <- lapply(1:ncol(ad), function(i) {
    sample_ad_ratios <- data.frame(sample = samples[i], chrom = chrom, 
                                   ratio = apply(ad[,,i], 1, function(x) x[2] / sum(x)))
    return(sample_ad_ratios)
  })
  ad_ratios <- do.call(rbind, ad_ratios_list)
  return(ad_ratios)
}

# Analyzing allele depth ratios per sample and chromosome for aneuploidy detection
detect_aneuploidy_per_sample_chrom <- function(ad_ratios) {
  # Assuming a simple threshold for demonstration purposes
  threshold <- 0.6
  ad_ratios$potential_aneuploidy <- ad_ratios$ratio < threshold | ad_ratios$ratio > (1 - threshold)
  
  # Summarizing potential aneuploidy by sample and chromosome
  aneuploidy_summary <- aggregate(potential_aneuploidy ~ sample + chrom, data = ad_ratios, FUN = sum)
  names(aneuploidy_summary)[3] <- "num_variants_deviant"
  
  return(aneuploidy_summary)
}

# Main function to run the analysis
run_aneuploidy_detection_per_sample_chrom <- function(vcf_path) {
  data <- extract_allele_depth_chrom_and_samples(vcf_path)
  ad_ratios <- calculate_ad_ratios_per_sample_chrom(data$ad, data$chrom, data$samples)
  aneuploidy_summary <- detect_aneuploidy_per_sample_chrom(ad_ratios)
  return(aneuploidy_summary)
}

# Example usage
vcf_path <- "path/to/your/vcf_file.vcf"
aneuploidy_summary <- run_aneuploidy_detection_per_sample_chrom(vcf_path)
print(aneuploidy_summary)