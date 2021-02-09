# Notes on getting BAIT to run

### Install BAIT
Install tarbell and unzip (https://sourceforge.net/p/bait/wiki/Home/)

`tar -jcvf BAIT_v1.0.tar.bz /usr/local/bin/BAIT_v1.4/`


### Package versions
* samtools version 1.4.1 (https://sourceforge.net/projects/samtools/files/samtools/1.4.1/)
* bedtools version 2.27.1 (https://github.com/arq5x/bedtools2/releases/tag/v2.27.1)

### Paths
Set PATH to bait executable

`export PATH=$PATH:"/path/to/executable"`

Set BAITPATH to bait.r file

`export BAITPATH="/usr/local/bin/BAIT_v1.4/"`

### Running

Run:

`BAIT -abcdr -i BAM_files -o OutputName -7 "Homo_sapiens" -8 "hg18"`
