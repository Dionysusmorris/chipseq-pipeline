# Rechecking input/output files for all R scripts 
# to ensure smooth transitions across scripts
#
# -----------------------------------------------------------------------------
#
# Have to use multiple software tools to generate a QC report analagous to multiqc
# or nf-core's chip-seq pipeline
# not sure whether to recreate the qc report or use only multiqc report
#
# -----------------------------------------------------------------------------
# git tag v0.1 
# -----------------------------------------------------------------------------
#
# This hbc_chipseq project reflects an updated project_startup.sh :
# 
# mkdir -p data/{raw,processed,reference} results/{qc,bam/sorted,callpeaks/idr,bed,counts,de, \
# figures,tables,tmp} meta scripts envs docs logs
# touch README.md .gitignore git_startup.txt docs/notes.md docs/workflow.md
#
# -----------------------------------------------------------------------------
#
# The pipeline currently requires two conda environments.
#
# chipseq_env:
# - alignment
# - QC
# - peak calling
#
# idr_env:
# - IDR reproducibility analysis
#
# The environments were separated because of dependency
# conflicts between IDR and the primary ChIP-seq software stack.
# 
# Packages not yet installed on chipseq_env but are required
# - homer
# - Picard
# - phantompeakqualtools
# - bedGraphToBigWig 
# 
# R environment files are not up to date on the repo currently
#
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
# 
# Large datasets, reference genomes, and analysis outputs are intentionally excluded from Git 
# tracking and remain local to the compute environment.