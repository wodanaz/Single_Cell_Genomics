# Single Cell Genomic


Create the environment

```bash
nano singlecell.yml
channels:
  - bih-cubi
dependencies:
  - bcl2fastq2
```

load Anaconda3 and create environment


```bash
module load Anaconda3
conda env update --file singlecell.yml --prefix /data/wraycompute/alejo/aleconda/singlecell
conda activate /data/wraycompute/alejo/aleconda/singlecell

```

Download data

```bash
module load ddsclient
ddsclient download -p Massri_7208 input 

```

Demultiplex NovaSeq Results:

```bash
module load Anaconda3
conda activate /data/wraycompute/alejo/aleconda/singlecell


cat do_novaseq2fastq.sh
#! /bin/bash -l
#SBATCH -J cellranger2fastq
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem 15G
cellranger mkfastq --run=/data/wraycompute/alejo/singlecell/input/210820_A00201R_0483_BHHV7YDRXY --csv=/data/wraycompute/alejo/singlecell/input/Lv_micro_samplesheet.csv --id=Lv_fastq_micro --output-dir=/data/wraycompute/alejo/singlecell/input/


```


```bash
sbatch do_novaseq2fastq.sh
```





