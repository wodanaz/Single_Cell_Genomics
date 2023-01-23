# Single Cell Genomics in Sea Urchins


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
conda env create --file singlecell.yml --prefix /data/wraycompute/alejo/aleconda/singlecell
conda activate /data/wraycompute/alejo/aleconda/singlecell

```

Download data

```bash
module load ddsclient
ddsclient download -p Massri_7208 input 

```

Demultiplex NovaSeq Results:

```bash
nano Lv_micro_samplesheet.csv
Lane,Sample,Index
*,Lv-6hpf,SI-TT-H6
*,Lv-8hpf,SI-TT-A7
*,Lv-10hpf,SI-TT-B7
*,Lv-12hpf,SI-TT-C7
*,Lv-14hpf,SI-TT-D7
*,Lv-16hpf,SI-TT-E7
*,Lv-18hpf,SI-TT-F7
*,Lv-20hpf,SI-TT-G7
```


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
cellranger mkref --genome=L_var_3_3 \
	               --fasta=/data/wraycompute/alejo/singlecell/input/genome/Lvar_scaffolds.fasta \
                 --genes=/data/wraycompute/alejo/singlecell/input/genome/Lvar.final.gtf \
                 --nthreads=8 --memgb=32 
```
```bash
sbatch do_ref_genome.sh
```



Map into the Lv genome. 

But first, we need to create a list of libraries to loop through

```bash
nano sample.list
Lv_7hpf_micro
Lv_9hpf_micro
Lv_11hpf_micro
Lv_13hpf_micro
Lv_15hpf_micro
Lv_17hpf_micro
```

CONTROL O + ENTER + CONTROL X 
	
Now, lets create multiple batch jobs to run in parallel	
	
```bash		 
for i in `cat sample.list`; do
echo '#!/usr/bin/env bash' > $i.fastq2counts.sh;
echo "#SBATCH -N 1" >> $i.fastq2counts.sh;
echo "#SBATCH -J fq2count.$i" >> $i.fastq2counts.sh;
echo "#SBATCH --mail-user=alebesc@gmail.com" >> $i.fastq2counts.sh;
echo "#SBATCH --mail-type=END,FAIL"  >> $i.fastq2counts.sh;
echo "#SBATCH --mem 15G" >> $i.fastq2counts.sh;
echo "cellranger count --id=${i} --transcriptome=/gpfs/fs1/data/covid19lab/L_var_3_3 --fastqs=/data/wraycompute/alejo/singlecell/input/HHV7YDRXY --sample=${i} --expect-cells=3000  star_parameters=/"--outFilterMatchNminOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNoverLmax=0.05 --outFilterMultimapNmax 0/"" >> $i.fastq2counts.sh;
done


for file in *fastq2counts.sh ; do sbatch $file ; done
                 
```




