#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=compute               #REQUIRED: which partition to use
#SBATCH --cpus-per-task=1                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB
#SBATCH --output=/projects/bgmp/apwat/bioinfo/Bi622/Demultiplex/Assignment-the-first/logs/base_qual_out_%j.log
#SBATCH --error=/projects/bgmp/apwat/bioinfo/Bi622/Demultiplex/Assignment-the-first/logs/base_qual_err_%j.log

conda activate bgmp-demux

DIR="/projects/bgmp/apwat/bioinfo/Bi622/Demultiplex"

/usr/bin/time -v $DIR/Assignment-the-first/qual_dist.py \
    -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
    -o $DIR/Assignment-the-first/histograms/read1 \
    -l read1

/usr/bin/time -v $DIR/Assignment-the-first/qual_dist.py \
    -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
    -o $DIR/Assignment-the-first/histograms/index1 \
    -l index1

/usr/bin/time -v $DIR/Assignment-the-first/qual_dist.py \
    -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
    -o $DIR/Assignment-the-first/histograms/index2 \
    -l index2

/usr/bin/time -v $DIR/Assignment-the-first/qual_dist.py \
    -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
    -o $DIR/Assignment-the-first/histograms/read2 \
    -l read2
