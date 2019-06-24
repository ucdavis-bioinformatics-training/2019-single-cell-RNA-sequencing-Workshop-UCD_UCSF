#!/bin/bash

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"

## Where cellranger and bcl2fastq executables are located
## a) by loading a module
module load cellranger/3.0.2

## b) or, by placing the location of the executables on the path (edit to your location)
# export PATH=/share/pathtosoftware/bin/:$PATH

## c) or if they are already on the path, do nothing

## Set the parameters for the run
transcriptome="/share/genomes/cellranger_genomes/refdata-cellranger-mm10-3.0.0"
basepath="/share/workshop/$USER/scrnaseq_example"
fastqs="${basepath}/00-RawData"

## provide the script the row # of the sample to be run
sample=`sed "$1q;d" samples.txt`

## https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome
## Create the call
call="cellranger count \
  --id=${sample} \
  --sample=${sample} \
  --transcriptome=${transcriptome} \
  --fastqs=${fastqs} \
  --localcores=2 \
  --localmem=64"

## Some other parameters that may be usefull/needed
## --expect-cells=NUM, number of cells expected
## --nosecondary, skip the unnecessary secondary analysis
## --r2-length=NUM, if your R2 qualities are really poor
## --chemistry=CHEM, should it fail chemistry detection

## Echo the call
echo $call
## Evaluate the call
#eval $call

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
