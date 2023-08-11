#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=compute               #REQUIRED: which partition to use
#SBATCH --cpus-per-task=4                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB
#SBATCH --output=/projects/bgmp/apwat/bioinfo/Bi622/Demultiplex/Assignment-the-third/logs/demuxing_out_%j.log
#SBATCH --error=/projects/bgmp/apwat/bioinfo/Bi622/Demultiplex/Assignment-the-third/logs/demuxing_err_%j.log

conda activate bgmp-demux

DIR="/projects/bgmp/apwat/bioinfo/Bi622/Demultiplex"

/usr/bin/time -v $DIR/Assignment-the-third/demuxer.py \
    -i /projects/bgmp/shared/2017_sequencing \
    -o $DIR/Assignment-the-third/FULL_OUTPUT \
    -x $DIR/Assignment-the-third/idx_list.txt \
    -q 26

#zip files
gzip $DIR/Assignment-the-third/FULL_OUTPUT/*.fastq