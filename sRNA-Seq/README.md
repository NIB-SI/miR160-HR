# sRNA-Seq

sRNA-Seq analysis was performed as described previously (Križnik et al., 2021).

- Reads were preprocessed to remove adapter sequences and low-quality reads (Phred quality score < 20) using [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html) and [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html). 
- Resulting reads were mapped with no mismatches to annotated potato miRNAs from [miRbase (release 22)](https://www.mirbase.org/), and to previously identified potato small RNAs (novel miRNAs, miRNA variants and phasiRNAs (Križnik et al., 2017)). 
- Detected sRNAs were then quantified and subjected to differential expression analysis using limma-voom approach as previously described (Križnik et al., 2021). 

<hr style="width: 40ch; border: 2px solid gray">


Check README files at each layer.
