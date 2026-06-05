# -----------------------------------------------------------------------------
# FUTURE OPTIMIZATION IDEA (NOT IMPLEMENTED)
#
# Goal:
#   Download SRA runs, merge technical replicates, and compress FASTQ output
#   without storing large intermediate FASTQ files on disk.
#
# Possible workflow:
#
#   SRR1 -----\
#              \             
#   SRR2 -------------> merged FASTQ stream --> pigz --> sample.fastq.gz
#              /
#   SRR3 -----/
#
# Challenges:
#   - fasterq-dump does not naturally support arbitrary sample renaming.
#   - Technical replicate metadata must be mapped from SRR IDs to sample names.
#   - Failure recovery is more difficult than the current workflow.
#   - Need to verify compatibility with paired-end data.
#   - Need to benchmark against the current download->FASTQ->compress workflow.
#
# Current production workflow:
#   1. Download each SRR independently.
#   2. Generate FASTQ files.
#   3. Rename/merge using sample metadata.
#   4. Compress with pigz.
#
# Example sketch only (NOT TESTED):
#
# parallel ... prefetch ... fasterq-dump ... pigz ...
#
# -----------------------------------------------------------------------------
