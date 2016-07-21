---
layout: post2
permalink: /install_tools_2016/
title: Tool Installation 2016 Student Page
header1: Tool Installation 2016
image: bioinformatics_LOGO.jpg
---

**The tools installation instruction is based on ubuntu server we used on Amazon cloud.**

Tools for HT-seq, RNA-seq, Cancer Genomics workshops
----------------------------------------------------

### Openjdk-7-jre-headless

```
sudo apt-get install openjdk-7-jre-headless
```

### Python 2

```
sudo apt-get install python2.7 python2.7-dev
sudo apt-get install python-numpy
sudo apt-get install python-scipy
```

### R, bioconductor and packages

```
sudo apt-get install r-base r-base-dev
source("`[`http://bioconductor.org/biocLite.R`](http://bioconductor.org/biocLite.R)`")
biocLite()
biocLite("gplots")
biocLite("ggplot2")
biocLite("cummeRbund")
biocLite("edgeR")
```

### annovar

Download from annovar website: <http://annovar.openbioinformatics.org/en/latest/user-guide/download/> (need registration)

```
tar -zxvf annovar.latest.tar.gz
```

Download databases

```
cd annovar
./annotate_variation.pl -buildver hg19 -webfrom annovar -downdb `<db_name>` ./humandb
```

### APOLLOH

#### MATLAB Runtime

```
wget `[`ftp://ftp.bcgsc.ca/public/shahlab/Apolloh/MCRInstaller.bin`](ftp://ftp.bcgsc.ca/public/shahlab/Apolloh/MCRInstaller.bin)
./MCRInstaller.bin -console
```

#### APOLLOH

```
unzip APOLLOH_0.1.0.zip
APOLLOH_0.1.1/run_apolloh.sh `<MCR location>
```

### bam-readcount

```
git clone --recursive `[`git://github.com/genome/bam-readcount.git`](git://github.com/genome/bam-readcount.git)
cmake bam-readcount/
make deps
make
#sudo make install
```

### bamtools

```
git clone `[`git://github.com/pezmaster31/bamtools.git`](git://github.com/pezmaster31/bamtools.git)
cd bamtools
mkdir build
cd build
cmake ..
make
cd ..
#sudo cp bin/bamtools /usr/local/bin/
```

### bedtools

```
git clone `[`https://github.com/arq5x/bedtools2.git`](https://github.com/arq5x/bedtools2.git)
cd bedtools2/
make
#sudo make install
```

### bowtie

```
wget `[`http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip`](http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip)
unzip bowtie-0.12.7-linux-x86_64.zip
#sudo mv bowtie2-2.2.4 /usr/local
#sudo ln -s /usr/local/bowtie2-2.2.4/bowtie2* /usr/local/bin
```

### bowtie2

```
wget `[`http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/bowtie2-2.2.4-linux-x86_64.zip`](http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/bowtie2-2.2.4-linux-x86_64.zip)
unzip bowtie2-2.2.4-linux-x86_64.zip
#sudo mv bowtie-0.12.7/ /usr/local/
#sudo ln -s /usr/local/bowtie-0.12.7/bowtie* /usr/local/bin/
```

### bwa

```
wget `[`http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2`](http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2)
tar -jxvf bwa-0.7.12.tar.bz2
cd bwa-0.7.12/
make
#sudo mv bwa /usr/local/bin
```

### chimerascan

**dependencies**

-   pysam
-   bowtie 0.12.7
-   cython

```
wget `[`https://chimerascan.googlecode.com/files/chimerascan-0.4.5a.tar.gz`](https://chimerascan.googlecode.com/files/chimerascan-0.4.5a.tar.gz)
tar -zxf chimerascan-0.4.5a.tar.gz
cd chimerascan-0.4.5/
sudo python setup.py install
```

### cufflinks

```
wget `[`http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz`](http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz)
tar -zxvf cufflinks-2.2.1.Linux_x86_64.tar.gz
```

### FastQC

```
wget `[`http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip)
unzip fastqc_v0.11.2.zip
cd FastQC
chmod +x fastqc
#sudo mv FastQC/ /usr/local/
#sudo ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc
```

### Flexbar

```
wget `[`http://downloads.sourceforge.net/project/flexbar/2.4/flexbar_v2.4_linux64.tgz`](http://downloads.sourceforge.net/project/flexbar/2.4/flexbar_v2.4_linux64.tgz)
tar -xzvf flexbar_v2.4_linux64.tgz
cd flexbar_v2.4_linux64/
chmod +r libtbb.so.2
sudo cp libtbb.so.2 /usr/lib
#sudo cp flexbar /usr/local/bin
```

### GATK

```
download from `[`https://www.broadinstitute.org/gatk/download`](https://www.broadinstitute.org/gatk/download)` (need to log in)
tar -jxf GenomeAnalysisTK-3.3-0.tar.bz2
```

### GMAP

```
wget `[`http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-12-23.tar.gz`](http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-12-23.tar.gz)
tar -zxf gmap-gsnap-2014-12-23.tar.gz
cd gmap-2014-12-23/
./configure
make
sudo make install
```

### HMMCopy

```
wget `[`http://compbio.bccrc.ca/files/2013/12/HMMcopy.zip`](http://compbio.bccrc.ca/files/2013/12/HMMcopy.zip)
unzip HMMcopy.zip
cd HMMcopy/
cmake .
make
#sudo mv bin/* /usr/local/bin/
```

### HTSeq

```
sudo pip install HTSeq
```

### Hydra

```
wget `[`https://hydra-sv.googlecode.com/files/Hydra.v0.5.3.tar.gz`](https://hydra-sv.googlecode.com/files/Hydra.v0.5.3.tar.gz)
tar -zxf Hydra.v0.5.3.tar.gz
cd Hydra-Version-0.5.3/
make clean
make all
#sudo cp bin/* /usr/local/bin
```

### lumpy

**dependencies**

-   bamtools
-   yaha

```
wget `[`https://github.com/arq5x/lumpy-sv/releases/download/0.2.9/lumpy-sv-0.2.9.tar.gz`](https://github.com/arq5x/lumpy-sv/releases/download/0.2.9/lumpy-sv-0.2.9.tar.gz)
tar -zxf lumpy-sv-0.2.9.tar.gz
cd lumpy-sv-0.2.9/
make
#sudo cp bin/lumpy /usr/local/bin
```

### MutationSeq

```
wget `[`ftp://ftp.bcgsc.ca/public/shahlab/MutationSeq/museq.zip`](ftp://ftp.bcgsc.ca/public/shahlab/MutationSeq/museq.zip)
unzip museq.zip
cd museq
wget `[`http://downloads.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.zip`](http://downloads.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.zip)
unzip boost_1_57_0.zip
make clean
python ./setup.py build --boost_source=./boost_1_57_0/
mv build/lib.linux-x86_64-2.7/pybam.so .
rm -rf boost*
cd ..
#sudo mv museq /usr/local
```

### OncoSNP-SEQ

**dependency MCR R2013b**

```
wget `[`http://www.mathworks.co.uk/supportfiles/downloads/R2013b/deployment_files/R2013b/installers/glnxa64/MCR_R2013b_glnxa64_installer.zip`](http://www.mathworks.co.uk/supportfiles/downloads/R2013b/deployment_files/R2013b/installers/glnxa64/MCR_R2013b_glnxa64_installer.zip)
unzip and make optionfile:
       destinationFolder=/opt/MCR
       agreeToLicense=yes
       mode=silent
./MCRInstaller.bin -console

git clone `[`https://github.com/cwcyau/oncosnpseq.git`](https://github.com/cwcyau/oncosnpseq.git)
```

### PennCNV

```
wget `[`http://www.openbioinformatics.org/penncnv/download/penncnv.latest.tar.gz`](http://www.openbioinformatics.org/penncnv/download/penncnv.latest.tar.gz)
sudo apt-get install libperl-dev
tar -zxf penncnv.latest.tar.gz
cd penncnv/kext
#sudo mv penncnv /usr/local/
```

### Picard

```
wget `[`https://github.com/broadinstitute/picard/releases/download/1.124/picard-tools-1.124.zip`](https://github.com/broadinstitute/picard/releases/download/1.124/picard-tools-1.124.zip)` -O picard-tools-1.124.zip
unzip picard-tools-1.124.zip
#sudo mv picard-tools-1.124 /usr/local
```

### RSEM

```
wget `[`http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.19.tar.gz`](http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.19.tar.gz)
tar -zxf rsem-1.2.19.tar.gz
make
make ebseq
#sudo mv rsem-1.2.19 /usr/local
```

### SAMStat

```
wget `[`http://downloads.sourceforge.net/project/samstat/samstat-1.5.tar.gz`](http://downloads.sourceforge.net/project/samstat/samstat-1.5.tar.gz)
tar -xzvf samstat-1.5.tar.gz
cd samstat-1.5/
./configure
make
#sudo make install
```

### samtools

```
wget `[`http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2/download`](http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2/download)` -O samtools-1.1.tar.bz2
tar -jxvf samtools-1.1.tar
cd samtools-1.1
make
#make install
```

### snpEff

```
wget `[`http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip`](http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip)
unzip snpEff_latest_core.zip
```

### STAR

```
wget `[`https://github.com/alexdobin/STAR/archive/STAR_2.4.0f1.tar.gz`](https://github.com/alexdobin/STAR/archive/STAR_2.4.0f1.tar.gz)
tar -zxvf STAR_2.4.0f1.tar.gz
cd STAR-STAR_2.4.0f1/source
make
#sudo mv STAR /usr/local/bin
```

### Strelka

```
wget `[`ftp://strelka:%27%27@ftp.illumina.com/v1-branch/v1.0.14/strelka_workflow-1.0.14.tar.gz`](ftp://strelka:%27%27@ftp.illumina.com/v1-branch/v1.0.14/strelka_workflow-1.0.14.tar.gz)
tar -zxf strelka_workflow-1.0.14.tar.gz
cd strelka_workflow-1.0.14/
./configure
make
#sudo make install
```

### TopHat

```
wget `[`http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.13.Linux_x86_64.tar.gz`](http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.13.Linux_x86_64.tar.gz)
tar -zxvf tophat-2.0.13.Linux_x86_64.tar.gz
#sudo mv tophat-2.0.13.Linux_x86_64 /usr/local
```

### YAHA

```
git clone `[`git://github.com/GregoryFaust/yaha.git`](git://github.com/GregoryFaust/yaha.git)
cd yaha
make
#sudo cp bin/yaha /usr/local/bin/
```

Tools for Metagenomics Workshop
-------------------------------

### Openjdk-7-jre-headless

sudo apt-get install openjdk-7-jre-headless

### python2

```
sudo apt-get install python2.7 python2.7-dev
sudo apt-get install python-numpy
sudo apt-get install python-scipy
```

### bioperl

```
 wget `[`http://bioperl.org/DIST/BioPerl-1.6.1.tar.gz`](http://bioperl.org/DIST/BioPerl-1.6.1.tar.gz)
 tar -zxf BioPerl-1.6.1.tar.gz
 cd BioPerl-1.6.1/
 perl Build.PL
 ./Build
 sudo ./Build install
```

### R

```
 add deb `[`http://cran.mtu.edu/bin/linux/ubuntu`](http://cran.mtu.edu/bin/linux/ubuntu)` trusty/ into /etc/apt/sources.list.d/R.list
 sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
 sudo apt-get update
 sudo apt-get install r-base r-base-dev
 install.packages(c('ape', 'vegan', 'getopt'))
```

### blat

```
 wget `[`http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat`](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat)
 chmod +x blat
 sudo mv blat /usr/local/bin
```

### blast+

```
 wget `[`ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.23/ncbi-blast-2.2.23+-x64-linux.tar.gz`](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.23/ncbi-blast-2.2.23+-x64-linux.tar.gz)
 tar -zxvf ncbi-blast-2.2.23+-x64-linux.tar.gz
 sudo mv ncbi-blast-2.2.23+ /usr/local
```

### bowtie

```
 wget `[`http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip`](http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip)
 unzip bowtie-0.12.7-linux-x86_64.zip
 sudo mv bowtie-0.12.7 /usr/local/
 sudo ln -s /usr/local/bowtie-0.12.7/bowtie* /usr/local/bin
```

### bowtie2

```
 wget `[`http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/bowtie2-2.2.4-linux-x86_64.zip`](http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/bowtie2-2.2.4-linux-x86_64.zip)
 unzip bowtie2-2.2.4-linux-x86_64.zip
 sudo mv bowtie2-2.2.4/ /usr/local/
 sudo ln -s /usr/local/bowtie2-2.2.4/bowtie* /usr/local/bin/
```

### bwa

```
 wget `[`http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.5a.tar.bz2`](http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.5a.tar.bz2)
 tar -xjvf bwa-0.7.5a.tar.bz2
 cd bwa-0.7.5a/
 make
 sudo mv bwa /usr/local/bin
```

### cdhit

```
 wget `[`https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz`](https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz)
 make openmp=yes
 mv into /usr/local/cd-hit-v4.6.1
 add /usr/local/cd-hit-v4.6.1 into path
```

### DETECT

```
 wget `[`http://www.compsysbio.org/projects/DETECT/detect_1.0.tar.gz`](http://www.compsysbio.org/projects/DETECT/detect_1.0.tar.gz)
 tar -zxvf detect_1.0.tar.gz
```

### DIAMOND

```
 wget `[`http://github.com/bbuchfink/diamond/releases/download/v0.7.9/diamond-linux64.tar.gz`](http://github.com/bbuchfink/diamond/releases/download/v0.7.9/diamond-linux64.tar.gz)
 tar -xzvf diamond-linux64.tar.gz
 sudo mv diamond /usr/local/bin
```

### EMBOSS

```
 sudo apt-get install libplplot-dev
 wget `[`ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz`](ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz)
 tar -zxvf EMBOSS-6.6.0.tar.gz
 cd EMBOSS-6.6.0/
 ./configure --without-x
 make
 sudo make install
```

### FastQC

```
 wget `[`http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip)
 unzip fastqc_v0.11.3.zip
 cd FastQC/
 chmod +x ./fastqc
 symlink into /usr/local/bin
```

### FASTX

```
 wget `[`https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz`](https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz)
 tar -zxvf libgtextutils-0.7.tar.gz
 cd libgtextutils-0.7/
 ./configure
 make
 sudo make install
 wget `[`https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2`](https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2)
 tar -xjvf fastx_toolkit-0.0.14.tar.bz2
 cd fastx_toolkit-0.0.14/
 ./configure
 make
 sudo make install
```

### FLASH

```
 wget `[`http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.7.tar.gz`](http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.7.tar.gz)
 tar -zxvf FLASH-1.2.7.tar.gz
 cd FLASH-1.2.7/
 make
 sudo mv ./flash /usr/local/bin
```

### HUMAnN

#### SCons (dependency)

```
 wget `[`https://bitbucket.org/biobakery/humann/downloads/humann-v0.99.tar.gz`](https://bitbucket.org/biobakery/humann/downloads/humann-v0.99.tar.gz)
 tar -xzvf humann-v0.99.tar.gz
 sudo mv humann-0.99/ /usr/local
```

### MetaPhlAn2

```
 hg clone `[`https://bitbucket.org/biobakery/metaphlan2`](https://bitbucket.org/biobakery/metaphlan2)
 sudo mv ./metaphlan2 /usr/local
```

### infernal

```
 wget `[`http://selab.janelia.org/software/infernal/infernal-1.1.1-linux-intel-gcc.tar.gz`](http://selab.janelia.org/software/infernal/infernal-1.1.1-linux-intel-gcc.tar.gz)
 tar -zxvf infernal-1.1.1-linux-intel-gcc.tar.gz
 cd infernal-1.1.1-linux-intel-gcc
 sudo mv ./binaries/* /usr/local/bin
```

### Microbiome Helper

#### Requirment

-   Both pipelines
    -   FastQC (v0.11.2) (optional)
    -   PEAR
-   Metagenomics
    -   Bowtie2
    -   Human pre-indexed database
    -   DIAMOND (&gt; v.0.7.0)
    -   MetaPhlAn2
    -   HUMAnN
-   16S
    -   FASTX toolkit (v0.0.14)
    -   QIIME (v1.9)
    -   SortMeRNA
    -   SUMACLUST
    -   PICRUSt
-   Visualization
    -   STAMP

```
 wget `[`https://github.com/mlangill/microbiome_helper/archive/master.zip`](https://github.com/mlangill/microbiome_helper/archive/master.zip)
 unzip master.zip
 sudo mv microbiome_helper-master/ /usr/local
```

### mothur

```
 wget `[`https://github.com/mothur/mothur/releases/download/v1.35.1/Mothur.cen_64.zip`](https://github.com/mothur/mothur/releases/download/v1.35.1/Mothur.cen_64.zip)
 untip and move into /usr/local and add into PATH
```

### PANDAseq

```
 sudo apt-add-repository ppa:neufeldlab/ppa
 sudo apt-get update
 sudo apt-get install pandaseq
```

### PEAR

```
 wget `[`http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.6-bin-64.tar.gz`](http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.6-bin-64.tar.gz)
 tar -xzvf pear-0.9.6-bin-64.tar.gz
 sudo mv ./pear-0.9.6-bin-64 /usr/local/bin
 sudo ln -s /usr/local/bin/pear-0.9.6-bin-64 /usr/local/bin/pear
```

### Phrap/Cross\_match

```
 tar -zxvf phrep.distrib.tar.Z
 make
 make manyreads
```

### PICRUSt

```
 sudo pip install -Iv biom-format==1.3.1
 wget `[`https://github.com/picrust/picrust/releases/download/1.0.0/picrust-1.0.0.tar.gz`](https://github.com/picrust/picrust/releases/download/1.0.0/picrust-1.0.0.tar.gz)
 tar -zxvf picrust-1.0.0.tar.gz
 cd picrust-1.0.0/
 sudo python ./setup.py install
```

### QIIME

```
 sudo apt-get install libfreetype6-dev
 sudo pip install qiime
```

### RDP

```
 wget `[`http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.2.zip`](http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.2.zip)
 unzip and move into /usr/local
```

### samtools

```
 wget `[`http://downloads.sourceforge.net/project/samtools/samtools/1.2/samtools-1.2.tar.bz2`](http://downloads.sourceforge.net/project/samtools/samtools/1.2/samtools-1.2.tar.bz2)
 tar -jxvf samtools-1.2.tar.bz2
 cd samtools-1.2/
 make
 sudo make install
```

### SCons

```
 wget `[`http://prdownloads.sourceforge.net/scons/scons-2.3.4.tar.gz`](http://prdownloads.sourceforge.net/scons/scons-2.3.4.tar.gz)
 tar -zxvf scons-2.3.4.tar.gz
 cd scons-2.3.4
 sudo python ./setup.py install
```

### SortMeRNA

```
 wget `[`http://bioinfo.lifl.fr/RNA/sortmerna/code/sortmerna-2.0-linux-64.tar.gz`](http://bioinfo.lifl.fr/RNA/sortmerna/code/sortmerna-2.0-linux-64.tar.gz)
 tar -zxvf sortmerna-2.0-linux-64.tar.gz
 cd sortmerna-2.0-linux-64
```

### STAMP

```
 sudo apt-get install freetype*
 sudo apt-get install python-qt4
 sudo pip install STAMP
```

### SUMACLUST

```
 svn co `[`http://www.grenoble.prabi.fr/public-svn/LECASofts/sumatra/tags/V_1.0.01`](http://www.grenoble.prabi.fr/public-svn/LECASofts/sumatra/tags/V_1.0.01) suma_package_V_1.0.01
 cd suma_package_V_1.0.01/sumaclust
 make
 sudo mv ./sumaclust /usr/local/bin
```

### Trimmomatic

```
 wget `[`http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip`](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip)
 unzip Trimmomatic-0.33.zip
```

### Trinity

```
 wget `[`https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.0.6.tar.gz`](https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.0.6.tar.gz)
 tar -zxvf v2.0.6.tar.gz
 cd trinityrnaseq-2.0.6/
 make
 make plugins
```

### USEARCH

```
 request from `[`http://drive5.com`](http://drive5.com)
 chmod +x usearch*
 sudo mv usearch* /usr/local/bin
```
