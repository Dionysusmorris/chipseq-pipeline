#!/bin/bash




# Create directories
cd ./results/callpeaks/
mkdir -p ./idr/


# Concordant True Peaks
  # Set IDR threshold to 0.05 or 0.01 depending on ChIP or TF-ChIP
  # and set MACS peak caller p-val to 1e-3

    #   macs3 callpeak \
    # -t "$chip" \
    # -c "$input" \
    # -f BAM -g mm \
    # -p 1e-3 \
    # -n "$name" \
    # --outdir ./callpeaks/ \
    # 2> ./callpeaks/${name}.log


parallel -j 3 '
  set -euo pipefail

  rep1="{}"
  base=$(basename "$rep1")
  name=${base%_sample1_peaks.narrowPeak}

  rep2="./${name}_sample2_peaks.narrowPeak"

  idr \
    --samples "$rep1" "$rep2" \
    --input-file-type narrowPeak \
    --rank signal.value \
    --idr-threshold 0.05 \ 
    --output-file ./idr/${name}_samples_idr.txt \
    --plot \ ${name}_samples_idr
    --log-output-file ../../logs/${name}_samples_idr.log
  echo "${name} sample IDR complete"
' ::: ./*_sample1_peaks.narrowPeak














'''
# Concordant Pooled Pseudo Peaks

# peak calling on IDRs should have pval 1e-3

date 

inputFile1=`basename $1`
treatFile1=`basename $2`
inputFile2=`basename $3`
treatFile2=`basename $4`
EXPT=$5

NAME1=`basename $treatFile1 _full.bam`
NAME2=`basename $treatFile2 _full.bam`

# Make Directories
mkdir -p /n/scratch2/mm573/idr_chipseq/macs
mkdir -p /n/scratch2/mm573/idr_chipseq/pooled_pseudoreps
mkdir -p /n/scratch2/mm573/idr_chipseq/tmp

# Set paths
baseDir=/n/groups/hbctraining/ngs-data-analysis-longcourse/chipseq/bowtie2
macsDir=/n/scratch2/mm573/idr_chipseq/macs
outputDir=/n/scratch2/mm573/idr_chipseq/pooled_pseudoreps
tmpDir=/n/scratch2/mm573/idr_chipseq/tmp

#Merge treatment BAMS
echo "Merging BAM files for pseudoreplicates..."
samtools merge -u ${tmpDir}/${NAME1}_${NAME2}_merged.bam $baseDir/${treatFile1} $baseDir/${treatFile2}
samtools view -H ${tmpDir}/${NAME1}_${NAME2}_merged.bam > ${tmpDir}/${EXPT}_header.sam

#Split merged treatments
nlines=$(samtools view ${tmpDir}/${NAME1}_${NAME2}_merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${tmpDir}/${NAME1}_${NAME2}_merged.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${EXPT}" # This will shuffle the lines in the file and split it
 into two SAM files
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}00 | samtools view -bS - > ${outputDir}/${EXPT}00.bam
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}01 | samtools view -bS - > ${outputDir}/${EXPT}01.bam

#Merge input BAMS
echo "Merging input BAM files for pseudoreplicates..."
samtools merge -u ${tmpDir}/${NAME1}input_${NAME2}input_merged.bam $baseDir/${inputFile1} $baseDir/${inputFile2}

#Split merged treatment BAM
nlines=$(samtools view ${tmpDir}/${NAME1}input_${NAME2}input_merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${tmpDir}/${NAME1}input_${NAME2}input_merged.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${EXPT}_input" # This will shuffle the lines in the file and split in two 
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}_input00 | samtools view -bS - > ${outputDir}/${EXPT}_input00.bam
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}_input01 | samtools view -bS - > ${outputDir}/${EXPT}_input01.bam


#Peak calling on pseudoreplicates
echo "Calling peaks for pseudoreplicate1 "
macs2 callpeak -t ${outputDir}/${EXPT}00.bam -c ${outputDir}/${EXPT}_input00.bam -f BAM -g hs -n $macsDir/${NAME1}_pr -B -p 1e-3  2> $macsDir/${NAME1}_pr_macs2.log

echo "Calling peaks for pseudoreplicate2"
macs2 callpeak -t ${outputDir}/${EXPT}01.bam -c ${outputDir}/${EXPT}_input01.bam -f BAM -g hs -n $macsDir/${NAME2}_pr -B -p 1e-3  2> $macsDir/${NAME2}_pr_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${NAME1}_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${NAME1}_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${NAME2}_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${NAME2}_pr_sorted.narrowPeak

#Independent replicate IDR
echo "Running IDR on pseudoreplicates..."
idr --samples $macsDir/${NAME1}_pr_sorted.narrowPeak $macsDir/${NAME2}_pr_sorted.narrowPeak --input-file-type narrowPeak --output-file ${EXPT}_pseudorep-idr --rank p.value --plot


# Remove the tmp directory
rm -r $tmpDir

'''


# Concordant Pseudo subreplicate peaks




'''
If starting with < 100K pre-IDR peaks for large genomes (human/mouse): 
For true replicates and self-consistency replicates an IDR threshold of 0.05 is more appropriate

Use a tighter threshold for pooled-consistency since pooling and subsampling equalizes
the pseudo-replicates in terms of data quality. 
Err on the side of caution and use more stringent IDR threshold of 0.01
'''





# Evaluating Self-Consistency Ratio 

  # Pseudo Subreplicate 1 / Subreplicate 2

  # Pooled Pseudo replicates / True Peakset



# peak IDR ranks are field $5. Found by int(-125log2(0.05))
# awk '{if($5 >= 540) print $0}' outputfile.txt | wc -l # 5% IDR peak counts
# awk '{if($5 >= 830) print $0}' outputfile.txt | wc -l # 1% IDR peak counts


# Cleanup

cd ../