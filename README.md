# scafN50

### Contents
scaflens.py : Python script to examine scaffold lengths and calculate N50  
scafN50.sh : wrapper shell script for scaflens.py with specific options to calculate statistics, including N50, after excluding contigs/scaffolds smaller than five length cutoffs (201 nt, 300 nt, 500 nt, 1,000 nt, and 10,000 nt)  

### Usage
$ git clone https://github.com/calacademy-research/scafN50.git  
$ ./scafN50.sh \<assembly_file.fasta\>  
\# Results output to a file named "scafN50s.txt"

### scaflens.py options
$ python scaflens.py  
```
Usage: scaflens.py <scafSeq_file_from_SOAP_assembly> [ [-t|-n [-C] [scaff_len ...]] | -r ]

         Outputs scaffold names and their sizes in a list (unless -t, -n or -r option used).
         You can pipe output to sort -nr -k2,2 to see the scaffolds listed largest to smallest.

         -t option shows total length of scaffolds (includes Ns). Unscaffolded singleton contigs are excluded (unless -C used).
         -n shows the N50 count along with the total.
         To see N50 and/or total count excluding scaffolds of a particular size or shorter
            you can add a number (or list of numbers) after the -t or -n option
         -C without this only scaffolds are included, using this also includes singleton contigs;
            when using -C 200 contigs/scaffolds 200 or smaller are ignored

         -r shows ranges of scaffold sizes and ranges of singleton contig sizes (other options ignored)


         Examples: scaflens.py soap_asm12.scafSefSeq
                   scaflens.py soap_asm12.fasta | sort -nr -k2,2
                   scaflens.py soap_asm12.fasta -n50 9999 99999
                   scaflens.py *.fasta -n -C 500 1000 10000 100000 1000000
                   scaflens.py -r *.fasta
```

### scafN50.sh description
$ ./scafN50.sh 
```

    scafN50.sh computes the scaffold N50 for the scafSeq file of argument 1.
    It computes N50s for all the scaffolds/contigs and also excluding scaffold of several sizes and below.
    The output is in a file named scafN50s.txt

```

### Citing

#### Authorship
Code author: James B. Henderson, jhenderson@calacademy.org  
README.md authors: Zachary R. Hanna, James B. Henderson  
