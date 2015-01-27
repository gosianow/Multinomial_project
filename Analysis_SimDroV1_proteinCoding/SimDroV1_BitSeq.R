# create bowtie reference index for the annotation
$ bowtie-build -f --ntoa ensSelect1.fasta ensSelect1-index
# align reads in data-c0b0.fastq against index
$ bowtie -q -v 3 -3 0 -p 4 -a -m 100 --sam ensSelect1-index \
data-c0b0.fastq data-c0b0.sam
