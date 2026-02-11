# scripts

## salmon quantification
```
for i in $(ls ./data | sed s/_[12].fastq.gz// | sort -u)	
do salmon quant \
	-i /DKED/salmon_index\
	-l ISR \
	-1 /data-repository/data/${i}_1.fastq.gz \
	-2 /data-repository/data/${i}_2.fastq.gz \
	--validateMappings \
	--numBootstraps 120 \
	-p 100 \
	-o /output/${i}
done
```

## limma voom
