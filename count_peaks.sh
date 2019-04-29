#!/bin/bash

# merge broad peaks from each samples
bedtools merge -i <(cat *_broad_peaks.bed | sort -k1,1 -k2,2n) -sorted | awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3}' | cut -f1,2,3 > Tn5_peaks_ALL.bed

# convert bed to gtf format
python /exec5/GROUP/barreiro/barreiro/barreiro_group/ALAIN/bed2gtf.py -i Tn5_peaks_ALL.bed -o Tn5_peaks.gtf
sed -i -e 's/peak_id/gene_id/' -e 's/peak/exon/' Tn5_peaks.gtf

# count number of reads using subread's featureCounts program
suffix='bam'
for f in *.$suffix; do f=${f/\.$suffix/}; echo "module load subread/1.4.6
featureCounts -p -P -d 19 -D 1000 -a Tn5_peaks.gtf $f.$suffix -o $f.counts.txt" > $f.sh ; chmod u+x $f.sh; qsub -A nnb-306-ad -lwalltime=1:0:0,nodes=1 -d `pwd` $f.sh; done

# tabulate read counts
paste -d "\t" $(ls -1v *.counts.txt) | sed '/^#/d' | cut -f 1,`seq --separator="," 7 7 126` | sed -e 's/.bam//g' -e 's/Geneid/peakID/g' > allSamples_Tn5.counts
