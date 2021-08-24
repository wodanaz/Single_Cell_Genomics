# Single Cell Genomic


To make sure it has installed:

```bash
#module load Anaconda3
#conda env update --file singlecell.yml --prefix /data/wraycompute/alejo/aleconda/singlecell
#conda activate /data/wraycompute/alejo/aleconda/singlecell

module load cellranger

```

Download data

```bash
module load ddsclient
ddsclient download -p Massri_7208 input 

```



