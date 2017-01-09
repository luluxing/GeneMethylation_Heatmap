# GeneMethylation_Heatmap

##'count_meth.py' 

Counting the methylation level from bisulfite-sequencing.

Setting TSS as coordiante zero, expanding to 5' and 3', dividing the upstream and the downstream sequences into BN bin numbers with BS bin size.

The input format is: Methylation file, GFF file, chromosome number, bin size, bin number and output prefix.

The methylation file is the transformation of the output of BRAT-nova 'acgt-count'.

##'matrix2image.py'

Drawing heatmaps based on the above counts.

The example image:

![alt tag](https://github.com/luluxing/GeneMethylation_Heatmap/blob/master/example.png)
