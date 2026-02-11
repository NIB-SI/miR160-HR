# scripts

## degradome sequencing

```
Mode: 1
Checking Dependencies
        bowtie: PASS: /usr/bin/bowtie
        bowtie-build: PASS: /usr/bin/bowtie-build
        RNAplex: PASS RNAplex 2.4.18
        GSTAr: PASS GSTAr.pl version 1.0
        samtools: PASS Version: 1.9-210-g72d140b (using htslib 1.9-437-g71d8683)


perl CleaveLand4.5.pl -e /DATA/... -u ./sRNAs.fa -n ./Unitato_CDS.fasta -c 3 -t
```

## _in silico_ target prediction

1. open https://www.zhaolab.org/psRNATarget/

2. Import ...\input\potato_sRNAs.fa and input\UniTato_cDNA.fa

3. Select:Schema V2 (2017 release)

  - Of top targets
  ```
  200
  Expectation:
  5
  Penalty for G:U pair:
  0.5
  Penalty for other mismatches:
  1
  Extra weight in seed region:
  1.5
  Seed region:
  2 -  13 NT
  ```

  - Of mismatches allowed in seed region
  ```
  2
  HSP size:
  19
  Allow bulge(gap)
  Penalty for opening gap:
  2
  Penalty for extending gap:
  0.5
  Calculate target accessibility
  Max UPE:
  25
  Flank length:
  17
   / 
  13
   NT in up/downstream
  ```

  - Translation inhibition range:
  ```
  10
   NT - 
  11
   NT
  ```

