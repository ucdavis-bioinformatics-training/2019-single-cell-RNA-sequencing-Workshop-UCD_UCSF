# Data Setup

Let's set up a project directory for the analysis, and talk a bit about project philosophy..

**1\.** First, create a directory for your user and the example project in the workshop directory:

    cd
    mkdir -p /share/workshop/$USER/scrnaseq_example

---

**2a\.** Next, go into that directory, create a raw data directory (we are going to call this 00-RawData) and cd into that directory. Lets then create symbolic links to the fastq files that contains the raw read data.

    cd /share/workshop/$USER/scrnaseq_example
    mkdir 00-RawData
    cd 00-RawData/
    ln -s /share/biocore/workshops/2019_scRNAseq/2017_10X_mouse_comparative_V3/cellranger-fastqs/HG27NBBXX/* .

This directory now contains a folder for each "sample" (in this case just 1) and the fastq files for each "sample" are in the sample folders.

**2b\.** lets create a sample sheet for the project, store sample names in a file called samples.txt

    ls > ../samples.txt
    cat ../samples.txt

---
**3a\.** Now, take a look at the raw data directory.

    ls /share/workshop/$USER/scrnaseq_example/00-RawData


**3b\.** To see a list of the contents of each directory.

    ls *

**3c\.** Lets get a better look at all the files in all of the directories.

    ls -lah */*

---

**5\.** View the contents of the files using the 'less' command, when gzipped used 'zless' (which is just the 'less' command for gzipped files, q to exit):

    cd 654/
    zless 654_S1_L008_R1_001.fastq.gz
    zless 654_S1_L008_R2_001.fastq.gz

Make sure you can identify which lines correspond to a single read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen. Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this use cat to output the uncompressed file and pipe that to "wc" to count the number of lines:

    zcat 654_S1_L008_R1_001.fastq.gz | wc -l

Divide this number by 4 and you have the number of reads in this file. One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

    zcat 654_S1_L008_R1_001.fastq.gz  | head -2 | tail -1

Note the header lines (1st and 3rd line) and sequence and quality lines (2nd and 4th) in each 4-line fastq block. Then, copy and paste the DNA sequence line into the following command (replace [sequence] with the line):

    echo -n [sequence] | wc -c

This will give you the length of the read.

Also can do the bash one liner:

    echo -n $(zcat 654_S1_L008_R1_001.fastq.gz  | head -2 | tail -1) | wc -c

See if you can figure out how this command works.
