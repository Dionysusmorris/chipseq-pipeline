# Generating the diffbind metadata csv 

#### GENERATE META DATA CSV #####

# Populate parameters and then
# execute this script to load bam and peaks into R env and generate diffbind metadata csv
bam_dir <- "results/bam/sorted"
peak_dir <- "results/callpeaks"

chip_bams <- list.files(
  bam_dir,
  pattern = "_chip\\.sorted\\.bam$",
  full.names = TRUE
)

sample_ids <- sub("_chip\\.sorted\\.bam$", "", basename(chip_bams))

metadata <- data.frame(
  SampleID = sample_ids,
  
  # Same tissue/cell type for all samples
  Tissue = "PRDM16",
  
  # ChIP target
  Factor = "H3K27ac",
  
  # WT vs cKO
  Condition = ifelse(grepl("^WT_", sample_ids), "WT", "cKO"),
  
  # REP1 -> 1, REP2 -> 2, etc.
  Replicate = sub(".*REP([0-9]+).*", "\\1", sample_ids),
  
  # ChIP BAM
  bamReads = chip_bams,
  
  # Matching input sample name
  ControlID = paste0(sample_ids, "_input"),
  
  # Input BAM
  bamControl = file.path(
    bam_dir,
    paste0(sample_ids, "_input.sorted.bam")
  ),
  
  # MACS3 peak file
  Peaks = file.path(
    peak_dir,
    paste0(sample_ids, "_peaks.narrowPeak")
  ),
  
  # DiffBind expects "macs" for MACS2/3
  PeakCaller = "macs",
  
  stringsAsFactors = FALSE
)

# Validate before writing
stopifnot(all(file.exists(metadata$bamReads)))
stopifnot(all(file.exists(metadata$bamControl)))
stopifnot(all(file.exists(metadata$Peaks)))

# Write to csv
write.csv(metadata,
          "meta/diffbind_metadata.csv",
          row.names = FALSE)
