# RNA-seq


- RNA-Seq, reads were quantified with Salmon using a transcriptome built from the S. tuberosum Phureja DM1-3 v6.1 genome (Pham et al., 2020) and the unified potato annotation (Zagorščak et al., 2024). 
- Differential expression analysis was conducted using [edgeR (v3.42.4)](https://bioconductor.org/packages//release/bioc/html/edgeR.html) and [limma (v3.56.2)](https://bioconductor.org/packages//release/bioc/html/limma.html).
- Salmon-derived counts were first scaled by the transcript-specific overdispersion estimates (function Overdispersion), then TMM-normalized and transformed using the voom function, followed by empirical Bayes moderation.

<hr style="width: 40ch; border: 2px solid gray">


Check README files at each layer.
