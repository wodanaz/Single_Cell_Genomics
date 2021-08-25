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


nano do_novaseq2fastq.sh
#! /bin/bash -l
#SBATCH -J cellranger2fastq
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem 15G
cellranger mkfastq --run=/data/wraycompute/alejo/singlecell/input/210820_A00201R_0483_BHHV7YDRXY \
                   --csv=/data/wraycompute/alejo/singlecell/input/Lv_micro_samplesheet.csv \
                   --id=Lv_fastq_micro \
                   --output-dir=/data/wraycompute/alejo/singlecell/input/


```


```bash
sbatch do_novaseq2fastq.sh
```

Create a Genome Reference


```bash
nano do_ref_genome.sh
#! /bin/bash -l
#SBATCH -J genomeref_cr
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem 15G
cellranger mkref --genome=L_var_3.0 \
	               --fasta=/data/wraycompute/alejo/singlecell/input/genome/Lvar_scaffolds.fasta \
                 --genes=/data/wraycompute/alejo/singlecell/input/genome/Lvar.final.gtf \
                 --nthreads=8 --memgb=32 
```
```bash
sbatch do_ref_genome.sh
```



Map into the Lv genome

```bash
nano do_fastq2counts.sh
#! /bin/bash -l
#SBATCH -J fastq2counts
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem 15G
cellranger count --id=Lv_7hpf_micro \
                 --transcriptome=/data/wraycompute/alejo/singlecell/input/genome/L_var_3.0 \
                 --fastqs=/data/wraycompute/alejo/singlecell/input/HHV7YDRXY/Lv_7hpf_micro \
                 --sample=Lv_7hpf_micro --expect-cells=3000
cellranger count --id=Lv-9hpf_micro \
                 --transcriptome=/data/wraycompute/alejo/singlecell/input/genome/L_var_3.0 \
                 --fastqs=/data/wraycompute/alejo/singlecell/input/HHV7YDRXY/Lv_9hpf_micro \
                 --sample=Lv_9hpf_micro  --expect-cells=3000 
cellranger count --id=Lv_11hpf_micro \
                 --transcriptome=/data/wraycompute/alejo/singlecell/input/genome/L_var_3.0 \
                 --fastqs=/data/wraycompute/alejo/singlecell/input/HHV7YDRXY/Lv_11hpf_micro \
                 --sample=Lv_11hpf_micro --expect-cells=3000 
cellranger count --id=Lv_13hpf_micro \
                 --transcriptome=/data/wraycompute/alejo/singlecell/input/genome/L_var_3.0 \
                 --fastqs=/data/wraycompute/alejo/singlecell/input/HHV7YDRXY/Lv_13hpf_micro \
                 --sample=Lv_13hpf_micro --expect-cells=3000 
cellranger count --id=Lv_15hpf_micro \
                 --transcriptome=/data/wraycompute/alejo/singlecell/input/genome/L_var_3.0 \
                 --fastqs=/data/wraycompute/alejo/singlecell/input/HHV7YDRXY/Lv_15hpf_micro \
                 --sample=Lv_15hpf_micro --expect-cells=3000 
cellranger count --id=Lv_17hpf_micro \
                 --transcriptome=/data/wraycompute/alejo/singlecell/input/genome/L_var_3.0 \
                 --fastqs=/data/wraycompute/alejo/singlecell/input/HHV7YDRXY/Lv_17hpf_micro \
                 --sample=Lv_17hpf_micro --expect-cells=3000 
                   
```



```bash
sbatch do_fastq2counts.sh
```



