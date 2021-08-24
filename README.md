# Single Cell Genomic


To make sure it has installed:

```bash
#module load Anaconda3
#conda env update --file singlecell.yml --prefix /data/wraycompute/alejo/aleconda/singlecell
#conda activate /data/wraycompute/alejo/aleconda/singlecell

module load cellranger
module load bcl2fastq2/v2.20.0.422-gcb01

```

Download data

```bash
module load ddsclient
ddsclient download -p Massri_7208 input 

```

Demultiplex NovaSeq Results:

```bash
do_bcl2fastq.sh 
#! /bin/bash -l

#SBATCH -J bcl2fastq
#SBATCH -o bcl2fastq_demux.log
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=64
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16


bcl2fastq \
    --runfolder-dir=$PWD/ \
    --output-dir=$PWD/fastqs \
    --loading-threads 4 \
    --processing-threads 8 \
    --writing-threads 4


```



