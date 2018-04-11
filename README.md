# fastq2bam 
Automated script to align paired-end FASTQ reads to Bam using BWA. 
Once Bam is created goes through several standard processing before Bam can be analyzed. 

* Picard MarkDuplicates
* GATK Indel Realignment
* GATK Base Recalibration

```
python fastq2bam.py -s test -f1 ./testData/test1.fq -f2 ./testData/test2.fq -o ./bam
```

