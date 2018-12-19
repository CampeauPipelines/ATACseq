#!/bin/batch

dir0=$1

cd $dir0

n=1
while [$n -le $(ls fastq_files |grep -e "fastq.gz" | wc -l)]; do

  ### STEP : TRIMMOMATIC

  cd $dir0
  mkdir -p trim/trimlog

  export trim_job_$n=$(echo " #!/bin/bash
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
            ${trimdir}/${SAMPLE}/${SAMPLE}_forward_paired.fastq.gz\
            ${trimdir}/${SAMPLE}/${SAMPLE}_forward_unpaired.fastq.gz\
            ${trimdir}/trim/${SAMPLE}/${SAMPLE}_reverse_paired.fastq.gz\
            ${trimdir}/${SAMPLE}/${SAMPLE}_reverse_unpaired.fastq.gz\
            ILLUMINACLIP:${dir0}/adapters.fasta:2:30:15:8:true\
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
 
 export bowtie_job_$n=$(echo " #!/bin/bash
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
  indexdir=${dir0}/bowtie_index/
  trimdir=${dir0}/trim/${SAMPLE}

  cd ${aligndir}
  mkdir ${SAMPLE} 
  cd ${SAMPLE}

  module load bowtie2

  bowtie2 -x ${indexdir} -q --local --phred33 --mm --threads 16 \
          -1 ${trimdir}/forward_paired.fastq.gz \
          -2 ${trimdir}/reverse_paired.fastq.gz \
    -S ${SAMPLE}_aligned.sam \ |
  sbatch --depend=afterok:trim_job_${n} | grep "[0-9]" | cut -d\ -f4)


  ### STEP : SAM_SORT
  
  cd $dir0
  mkdir -p sorted_bam/coordinate
  mkdir -p sorted_bam/queryname
  
  export sort_job_$n=$( echo "#!/bin/bash
    #SBATCH --time=5:00:00
    #SBATCH --job-name=sort_bam_$i
    #SBATCH --output=%x_%j.out
    #SBATCH --error=%x_%j.err
    #SBATCH --ntasks=16
    #SBATCH --mem=50000
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=sophie.ehresmann@gmail.com

    SAMPLE=$1
    sort_dir=${dir0}/sorted_bam/
    bam_dir=${dir0}/bowtie/${SAMPLE}

    cd ${sort_dir}

    module load java/1.8.0_121a picard/2.10.7 samtools/1.5 &&\

      samtools view -S -b ${bam_dir}/${SAMPLE}_aligned.sam > ${bam_dir}/${SAMPLE}_aligned.bam
      rm ${bam_dir}/${SAMPLE}_aligned.sam

      java -jar $EBROOTPICARD/picard.jar SortSam\
      I=${bam_dir}/${SAMPLE}_aligned.bam\
      O=${sort_dir}/coordinate/${SAMPLE}_sorted_$2.bam\
      SORT_ORDER=coordinate
      
      java -jar $EBROOTPICARD/picard.jar SortSam\
      I=${bam_dir}/${SAMPLE}_aligned.bam\
      O=${sort_dir}/queryname/${SAMPLE}_sorted_$2.bam\
      SORT_ORDER=queryname      
      
  
(( n++ ))
