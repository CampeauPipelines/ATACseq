#!/bin/batch
## bash atacseq.sh <working_dir> <fastq_dir> <adapters>
## Directories should be called by their full addresses 
## Aligns by defaults to human hg19, modify the script for other index


dir0=$1
fastq_dir=$2
adapters=$3

cd $dir0
mkdir -p trim/trimlog
mkdir -p bowtie/bowtie_index
mkdir -p sorted_bam/
mkdir -p mark_dup/metrics
mkdir -p final_bam/stats
mkdir peaks
mkdir final_stats

trimdir=${dir0}/trim
trimlogdir=${trimdir}/trimlog
aligndir=${dir0}/bowtie
indexdir=${aligndir}/bowtie_index
sortdir=${dir0}/sorted_bam
dupdir=${dir0}/mark_dup
metricsdir=${dupdir}/metrics
finalbamdir=${dir0}/final_bam
statsdir=${finalbamdir}/stats
peakdir=${dir0}/peaks
finalstatsdir=${dir0}/final_stats


cd bowtie/bowtie_index 
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip ### Modify here for a different alignment index
unzip *.zip
 
cd $dir0

n=1

for i in $(ls $fastq_dir | grep -e "fastq.gz"); do

  SAMPLE=$(echo $i | awk -F ".fastq.gz" '{print $1}') ##Get name of sample

  ### STEP : TRIMMOMATIC

  export trim_job_${n}=$(echo "#!/bin/bash
    #SBATCH --time=12:00:00
    #SBATCH --job-name=Trimmomatic_$SAMPLE
    #SBATCH --output=%x-%j.out
    #SBATCH --error %x-%j.err
    #SBATCH --ntasks=12
    #SBATCH --mem-per-cpu=2000
    #SBATCH --mail-user=$JOB_MAIL
    #SBATCH --mail-type=END, FAIL
    #SBATCH --A=$SLURM_ACCOUNT    
    
    cd ${fastq_dir}
    fastqs=$(ls | grep -e "${SAMPLE}")  ## Get forward and reverse fastqs

    module load  java trimmomatic && \
            java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 \
            -trimlog ${trimlogdir}/${SAMPLE}_trim.txt\
            $fastqs\
            ${trimdir}/${SAMPLE}_forward_paired.fastq.gz\
            ${trimdir}/${SAMPLE}_forward_unpaired.fastq.gz\
            ${trimdir}/${SAMPLE}_reverse_paired.fastq.gz\
            ${trimdir}/${SAMPLE}_reverse_unpaired.fastq.gz\
            ILLUMINACLIP:${adapters}:2:30:15:8:true\
            TRAILING:30\
            MINLEN:32\ " |
    sbatch | grep "[0-9]" | cut -d\ -f4)
          


      ### STEP: BOWTIE
 
 cd $dir0 
 
 export bowtie_job_${n}=$(echo "#!/bin/bash
  #SBATCH --time=12:00:00
  #SBATCH --job-name=bowtie_align${i}
  #SBATCH --output=%x-%j.out
  #SBATCH --error %x-%j.err
  #SBATCH --ntasks=16
  #SBATCH --mem=80000
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=$JOB_MAIL
  #SBATCH --A=$SLURM_ACCOUNT

  cd ${aligndir}

  module load bowtie2

  bowtie2 -x ${indexdir} -q --local --phred33 --mm --threads 16 \
          -1 ${trimdir}/${SAMPLE}_forward_paired.fastq.gz \
          -2 ${trimdir}/${SAMPLE}_reverse_paired.fastq.gz \
    -S ${SAMPLE}_aligned.sam \ |
  sbatch --depend=afterok:trim_job_${n} | grep "[0-9]" | cut -d\ -f4)


  ### STEP : SAM_SORT
  
  cd $dir0
  
  export sort_job_${n}=$(echo "#!/bin/bash
    #SBATCH --time=5:00:00
    #SBATCH --job-name=sort_bam_${i}
    #SBATCH --output=%x_%j.out
    #SBATCH --error=%x_%j.err
    #SBATCH --ntasks=16
    #SBATCH --mem=50000
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=$JOB_MAIL


    cd ${sortdir}

    module load java picard samtools &&\

      samtools view -S -b ${aligndir}/${SAMPLE}_aligned.sam > ${aligndir}/${SAMPLE}_aligned.bam
      rm ${aligndir}/${SAMPLE}_aligned.sam

      java -jar $EBROOTPICARD/picard.jar SortSam\
      I=${aligndir}/${SAMPLE}_aligned.bam\
      O=${sortdir}/${SAMPLE}_sorted.bam\
      SORT_ORDER=coordinate" |
  sbatch --depend=afterok:bowtie_job${$n} | grep "[0-9]" | cut -d\ -f4)

  ### STEP : MARK DUPLICATES
  
  cd $dir0
  
  export dup_job_${n}=$(echo "#!/bin/bash
  #SBATCH --time=5:00:00
  #SBATCH --job-name=mark_dup_${i}
  #SBATCH --output=%x_%j.out
  #SBATCH --error=%x_%j.err
  #SBATCH --ntasks=16
  #SBATCH --mem=50000
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=$JOB_MAIL
  
  
  cd ${dupdir}
  
  module load java picard &&\
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates\
    I=${sort_dir}/${SAMPLE}_sorted.bam\
    O=${SAMPLE}_markdup.bam\
    M=${metrics_dir}/markdup_metrics_${SAMPLE}.txt\
    REMOVE_DUPLICATES=false" |
    sbatch --depend=afterok:sort_job_${n} | grep "[0-9]" | cut -d\ -f4)
    
  ### STEP: REMOVE FLAGS
  
  cd $dir0
  
  export flag_job_${n}=$(echo "#!/bin/bash
  #SBATCH --time=12:00:00
  #SBATCH --job-name=remove_flags_${i}
  #SBATCH --output=%x_%j.out
  #SBATCH --error=%x_%j.err
  #SBATCH --ntasks=12
  #SBATCH --mem=20000
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=$JOB_MAIL

  cd ${finalbamdir}

  module load samtools\
    samtools view -F 1804\
       ${dup_dir}/${SAMPLE}_markdup.bam>${SAMPLE}_final.bam

    wait 1

    samtools index\
           ${SAMPLE}_final.bam\
           ${SAMPLE}_final.bam.bai

    wait 1

    samtools flagstat\
      ${SAMPLE}_final.bam>${statsdir}/${SAMPLE}_final_bam_stats.txt" |
    sbatch --depend=afterok:dup_job_${n} | grep "[0-9]" | cut -d\ -f4)
    
   ### STEP : PEAK CALLING 
   
  cd $dir0
 
  export peak_job_${n}=$(echo "#!/bin/bash
  #SBATCH --time=6:00:00
  #SBATCH --job-name=peak_calling_${i}
  #SBATCH --output=%x_%j.out
  #SBATCH --error=%x_%j.err
  #SBATCH --ntasks=6
  #SBATCH --mem=20000
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=$JOB_MAIL
  
  
  cd peakdir
  
  module load mugqic/MACS2
  
  macs2 callpeak -t ${finalbamdir}/${SAMPLE}_final.bam --keep-dup all \
  --broad --nomodel --extsize 200 --nolambda -n $SAMPLE -f BAMPE " |
  sbatch --depend=afterok:flag_job_${n} | grep "[0-9]" | cut -d\ -f4)
  
#STEP: QC

cd $dir0

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


  cd ${finalstatsdir}

  module load samtools

  for f in $(ls ${finalbamdir} | grep ".bam"); do
    file=${finalbamdir}/$f
    samtools flagstat $file > ${SAMPLE}_final_stats.txt
    CHROMOSOMES=$(samtools view -H $file | grep '^@SQ' | cut -f 2 | grep -v -e _ -e MT -e 'VN:' -e "GL" | sed 's/SN://' | xargs echo)
    samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 10 -@ 12 $file $CHROMOSOMES > ${SAMPLE}.mapped.bam
    samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 10 -@ 12 $file MT > ${SAMPLE}.chrM.bam

  done

  for f in $(ls *.mapped.bam); do
          samtools view -@ 12 $f | cut -f1 | sort -u | wc -l
  done > processed_reads.txt

  for f in $(ls *.chrM.bam); do
          samtools view -@ 12 $f | cut -f1 | sort -u | wc -l
  done > processed_reads_MT.txt " |
 sbatch --depend=afterok:flag_job_${n} | grep "[0-9]" | cut -d\ -f4)

n=n+1

done
