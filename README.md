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



