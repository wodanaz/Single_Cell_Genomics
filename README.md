# Single Cell Genomics

First we need to install the envinronment for single cell genomics


```bash
nano singlecell.yml


```
To make sure it has installed:

```bash
module load Anaconda3
conda env update --file singlecell.yml --prefix /data/wraycompute/alejo/aleconda/singlecell
conda activate /data/wraycompute/alejo/aleconda/singlecell

cellranger

```

The output should be similiar to the following:

```bash
  /xxx/xxx/user.xxx/xxx/xxx/cellranger-3.1.0/cellranger-cs/3.1.0/bin
   cellranger  (3.1.0)
   Copyright (c) 2019 10x Genomics, Inc.  All rights reserved.
```

# Perform a Sitecheck
The purpose of sitecheck is to check your system to make sure it meets the system requirements for running the cellranger pipeline. Run the command and use the > symbol to direct the output to a file.

```bash
cellranger sitecheck > sitecheck.txt
```

Now use the less command to take a look at this file. Use the up and down arrow keys to scroll through the file. Or spacebar to scroll down by page. Press the q key on the keyboard to quit out of the less program.
less sitecheck.txt

We will take a look at the following sections of the sitecheck file to see if they meet the minimum requirements for running Cell Ranger:


