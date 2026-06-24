#!/bin/bash
# Not yet functioning

# # Run picard CollectAlignmentSummaryMetrics for a sample
# java -jar picard.jar CollectAlignmentSummaryMetrics \
#   --INPUT $INPUT_BAM \
#   --REFERENCE_SEQUENCE $REFERENCE \
#   --OUTPUT $OUTPUT_METRICS_FILE

# # Run picard CollectAlignmentSummaryMetrics for a sample
# java -jar picard.jar CollectAlignmentSummaryMetrics \
#   --INPUT $INPUT_BAM \
#   --REFERENCE_SEQUENCE $REFERENCE \
#   --OUTPUT $OUTPUT_METRICS_FILE

# # Run phantompeakqualtools for a sample
# R_LIBS_USER=/n/groups/hbctraining/phantompeakqualtools/ \
#  Rscript --no-environ /n/groups/hbctraining/phantompeakqualtools/run_spp.R \
#  -c="$COORDINATE_SORTED_BAM" \
#  -savp="$OUTPUT_PDF" \
#  -savd="$OUTPUT_RDATA" \
#  -out="$OUTPUT_FILE" \
#  -p=$CORES

# load("SAMPLE_PHANTOMPEAKQUALTOOLS_OUTPUT.Rdata")
# crosscorr$phantom.coeff

# load("SAMPLE_PHANTOMPEAKQUALTOOLS_OUTPUT.Rdata")
# crosscorr$rel.phantom.coeff

# # Calculate the fraction of reads in peaks
# # Requires bedtools to run
# sh calculate_frip.sh \
#  $DIRECTORY_WITH_BAM_FILES \
#  $DIRECTORY_WITH_PEAK_FILES \
#  $OUTPUT_FILE

# # Mark duplicates and create metrics file
# java -jar picard.jar MarkDuplicates \
#   --INPUT <SORTED_BAM_FILE> \
#   --OUTPUT <REMOVE_DUPLICATES_BAM_FILE> \
#   --METRICS_FILE <METRICS_FILE>