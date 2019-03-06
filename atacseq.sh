#!/bin/batch

dir0=$1

cd $dir0

n=1
while [$n -le $(ls fastq_files |grep -e "fastq.gz" | wc -l)]; do

  ### STEP : TRIMMOMATIC

  cd $dir0
  mkdir -p trim/trimlog

  export trim_job_${n}=$(echo "#!/bin/bash
    #SBATCH --time=12:00:00
    #SBATCH --job-name=Trimmomatic_$i
    #SBATCH --output=%x-%j.out
    #SBATCH --error %x-%j.err
    #SBATCH --ntasks=12
    #SBATCH --mem-per-cpu=2000
    #SBATCH --mail-user=$JOB_MAIL
    #SBATCH --mail-type=END, FAIL
    #SBATCH --A=$SLURM_ACCOUNT

    SAMPLE=$i

    FASTQ_DIR=${dir0}/QD303
    trimdir=${dir0}/trim
    trimlogdir=${trimdir}/trimlog

    cd ${trimdir}
    mkdir ${SAMPLE}

    module load  java/1.8.0_121 trimmomatic/0.36 && \
            java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 \
            -trimlog ${trimlogdir}/${SAMPLE}_trim.txt\
            ${FASTQ_DIR}/${SAMPLE}/*.fastq.gz\
            ${trimdir}/${SAMPLE}_forward_paired.fastq.gz\
            ${trimdir}/${SAMPLE}_forward_unpaired.fastq.gz\
            ${trimdir}/${SAMPLE}_reverse_paired.fastq.gz\
            ${trimdir}/${SAMPLE}_reverse_unpaired.fastq.gz\
            ILLUMINACLIP:${dir0}/adapters.txt:2:30:15:8:true\
            TRAILING:30\
            MINLEN:32\" |
    sbatch | grep "[0-9]" | cut -d\ -f4)
          


 ### STEP: BOWTIE
 cd $dir0
 mkdir -p bowtie/bowtie_index
 cd bowtie/bowtie_index
 
 wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip
 unzip *.zip
 
 cd $dir0 
 
 export bowtie_job_${n}=$(echo "#!/bin/bash
  #SBATCH --time=12:00:00
  #SBATCH --job-name=bowtie_align$i
  #SBATCH --output=%x-%j.out
  #SBATCH --error %x-%j.err
  #SBATCH --ntasks=16
  #SBATCH --mem=80000
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=$JOB_MAIL
  #SBATCH --A=$SLURM_ACCOUNT
  
  SAMPLE=$i
  aligndir=${dir0}/bowtie
  indexdir=${dir0}/bowtie_index
  trimdir=${dir0}/trim

  cd ${aligndir}

  module load bowtie2

  bowtie2 -x ${indexdir} -q --local --phred33 --mm --threads 16 \
          -1 ${trimdir}/${SAMPLE}_forward_paired.fastq.gz \
          -2 ${trimdir}/${SAMPLE}_reverse_paired.fastq.gz \
    -S ${SAMPLE}_aligned.sam \ |
  sbatch --depend=afterok:trim_job_${n} | grep "[0-9]" | cut -d\ -f4)


  ### STEP : SAM_SORT
  
  cd $dir0
  mkdir -p sorted_bam/coordinate
  mkdir -p sorted_bam/queryname
  
  export sort_job_${n}=$(echo "#!/bin/bash
    #SBATCH --time=5:00:00
    #SBATCH --job-name=sort_bam_$i
    #SBATCH --output=%x_%j.out
    #SBATCH --error=%x_%j.err
    #SBATCH --ntasks=16
    #SBATCH --mem=50000
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=$JOB_MAIL

    SAMPLE=$i
    sort_dir=${dir0}/sorted_bam
    bam_dir=${dir0}/bowtie

    cd ${sort_dir}

    module load java/1.8.0_121a picard/2.10.7 samtools/1.5 &&\

      samtools view -S -b ${bam_dir}/${SAMPLE}_aligned.sam > ${bam_dir}/${SAMPLE}_aligned.bam
      rm ${bam_dir}/${SAMPLE}_aligned.sam

      java -jar $EBROOTPICARD/picard.jar SortSam\
      I=${bam_dir}/${SAMPLE}_aligned.bam\
      O=${sort_dir}/coordinate/${SAMPLE}_sorted_coordinate.bam\
      SORT_ORDER=coordinate" |
  sbatch --depend=afterok:bowtie_job${$n} | grep "[0-9]" | cut -d\ -f4)

  ### STEP : MARK DUPLICATES
  
  cd $dir0
  mkdir -p mark_dup/metrics
  
  export dup_job_${n}=$(echo "#!/bin/bash
  #SBATCH --time=5:00:00
  #SBATCH --job-name=mark_dup_${i}
  #SBATCH --output=%x_%j.out
  #SBATCH --error=%x_%j.err
  #SBATCH --ntasks=16
  #SBATCH --mem=50000
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=$JOB_MAIL
  
  SAMPLE=$i
  dup_dir=${dir0}/mark_dup
  sort_dir=${dir0}/sorted_bam
  metrics_dir=${dup_dir}/metrics

  cd ${dup_dir}
  mkdir ${metrics_dir}

  module load java/1.8.0_121a picard/2.10.7 &&\
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates\
    I=${sort_dir}/coordinate/${SAMPLE}_sorted_coordinate.bam\
    O=${SAMPLE}_marked_duplicates_coordinate.bam\
    M=${metrics_dir}/marked_duplicates_metrics_coordinate_${SAMPLE}.txt\
    REMOVE_DUPLICATES=false" |
    sbatch --depend=afterok:sort_job_${n} | grep "[0-9]" | cut -d\ -f4)
    
  ### STEP: REMOVE FLAGS
  
  cd $dir0
  mkdir -p final_bam/stats
  
  export flag_job_${n}=$(echo "#!/bin/bash
  #SBATCH --time=12:00:00
  #SBATCH --job-name=remove_flags_${i}
  #SBATCH --output=%x_%j.out
  #SBATCH --error=%x_%j.err
  #SBATCH --ntasks=12
  #SBATCH --mem=20000
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=$JOB_MAIL

  SAMPLE=$i
  dup_dir=${dir0}/mark_dup/${SAMPLE}
  final_bam_dir=${dir0}/final_bam
  stats_dir=${final_bam_dir}/stats

  cd ${final_bam_dir}

  module load samtools/1.5\
    samtools view -F 1804\
          ${dup_dir}/marked_duplicates_coordinate_${SAMPLE}.bam>${SAMPLE}_final.bam

    wait 1

    samtools index\
           ${SAMPLE}_final.bam\
           ${SAMPLE}_final.bam.bai

    wait 1

    samtools flagstat\
      ${SAMPLE}_final.bam>${stats_dir}/${SAMPLE}_final_bam_stats.txt" |
    sbatch --depend=afterok:dup_job_${n} | grep "[0-9]" | cut -d\ -f4)
    
   ### STEP : PEAK CALLING 
   
  cd $dir0
  mkdir -p peaks
  
  export peak_job_${n}=$(echo "#!/bin/bash
  #SBATCH --time=6:00:00
  #SBATCH --job-name=peak_calling_${i}
  #SBATCH --output=%x_%j.out
  #SBATCH --error=%x_%j.err
  #SBATCH --ntasks=6
  #SBATCH --mem=20000
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=$JOB_MAIL

  SAMPLE=$i
  
  final_bam_dir=${dir0}/final_bam
  peak_dir=${dir0}/peaks
  
  cd peak_dir
  
  module load mugqic/MACS2
  
  macs2 callpeak -t ${final_bam_dir}/${SAMPLE}_final.bam --keep-dup all --broad --nomodel --extsize 200 --nolambda -n $SAMPLE -f BAMPE " |
  sbatch --depend=afterok:flag_job_${n} | grep "[0-9]" | cut -d\ -f4)
  
#STEP: QC

cd $dir0
mkkdir stats

export stat_job_${n}=$(echo "#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --job-name=stats_atac_test_$(i)
#SBATCH --output=%x-%j.out
#SBATCH --error %x-%j.err
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-user=$JOB_MAIL
#SBATCH --mail-type=END, FAIL
#SBATCH --A=$SLURM_ACCOUNT


bam_dir=${dir0}/final_bam
cd ${dir0}/stats

module load samtools

for f in $(ls $bam_dir | grep ".sh"); do
  file=${bam_dir}/$f
  samtools flagstat $file > ${f}_stats.txt
  CHROMOSOMES=$(samtools view -H $file | grep '^@SQ' | cut -f 2 | grep -v -e _ -e MT -e 'VN:' -e "GL" | sed 's/SN://' | xargs echo)
  samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 10 -@ 12 $file \$CHROMOSOMES > $f.mapped.bam
  samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 10 -@ 12 $file MT > $f.chrM.bam

done

for f in $(ls *.mapped.bam); do
        samtools view -@ 12 $f | cut -f1 | sort -u | wc -l
done > processed_reads.txt

for f in $(ls *.chrM.bam); do
        samtools view -@ 12 $f | cut -f1 | sort -u | wc -l
done > processed_reads_MT.txt

  
(( n++ ))

done
