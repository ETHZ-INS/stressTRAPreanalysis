---
title: "Trap Reanalysis"
output:
  html_document:
    keep_md: yes
    theme: united
    toc: yes
  pdf_document:
    toc: yes
---



required r libraries

```r
library(edgeR)
library(SEtools)
library(SummarizedExperiment)
```


## About this reanalysis

Here, we outline the re-analysis of two GEO datasets, GSE100579 and GSE131972, both of which are investigateing the effect of acute stress on the translatome of CA3 pyramidal neurons in the mouse hippocampus.

the reanalysis be performed on the output of two different alginment methods a genome alignment using or a pseudoalignment to the transcriptome using kallisto

first, we are loading both the kallisto and the genome alignment data


```r
kallistodata <- readRDS("data/TRAPReanalysis.SE.rds")
```

## Re-Analysis of Marrocco et al. 2019

first, lets run a reanalysis of Marrocco et al. 2019 checking for ELS dependent changes in the acute stress response

```r
se <- kallistodata
se <- subset(se,select  = se$Set == "GSE131972")
```

In their Publications the authors unfortunately do not upload a list with differentially expressed genes. However, in their discussion they mention a number of genes that they claim are differentially expressed between ELS and non-ELS mice after acute stress.
They claim that acute stress reduces the expression of Grin1, Grin2a, Gabbr2, and Gabra1 in CA3 neurons of non-ELS mice, but not ELS mice:

```r
sehm(se,c("Grin1","Grin2a","Gabbr2","Gabra1"),do.scale = T,anno_columns = c("ELS","Condition"), cluster_rows = T)
```

```
## Using assay logcpm
```

![](README_figs/README-unnamed-chunk-4-1.png)<!-- -->

it becomes apparent, that this finding does not look very solid

Further, they claim that there is restricted list of genes selectively induced by AS in ELS mice (Per1, Npy, Nfkbia, Penk,Dusp1, Cst3, Trib1, Htra1, Sdc4, Plekhf1) but not non-ELS mice.

```r
sehm(se,c("Per1", "Npy", "Nfkbia", "Penk","Dusp1", "Cst3", "Trib1", "Htra1", "Sdc4", "Plekhf1"),do.scale = T,anno_columns = c("ELS","Condition"), cluster_rows = T)
```

```
## Using assay logcpm
```

![](README_figs/README-unnamed-chunk-5-1.png)<!-- -->

it becomes apparent, that this finding also does not look very solid. while these genes might be increased in expression in ELS mice it is hard to see this being a significant finding, especially, since the same effect can be observed for one replicate in the non-ELS group.

The authors claim that there are a number of genes that appear to be induced by AS in both ELS and non ELS mice. these include (Egr1/2/4, Arc, Fos, and Fosb)

```r
sehm(se,c("Egr1", "Egr2", "Egr4", "Arc","Fos", "Fosb"),do.scale = T,anno_columns = c("ELS","Condition"), cluster_rows = T)
```

```
## Using assay logcpm
```

![](README_figs/README-unnamed-chunk-6-1.png)<!-- -->

For these genes it looks as if their finding might be significant.

Let's investigate the whole data with a statistical approach using egeR and a GLM type of analysis for swim effects, early life effects and interactions

```r
#experimental design, interactive model
design <- model.matrix(~se$Condition * se$ELS)

y <- DGEList(counts=assays(se)$counts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)

#filter out genes that are below 10 counts in more than 75% of samples
keep <- rowSums(y$counts>10) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]

Results <- list()
fit <- glmQLFit(y,design)
for(i in colnames(design)[-1]){
  Results[[i]] <- glmQLFTest(fit, i)
}
```

Les's investigate if there are any genes altered by acute stress

```r
topTags(Results$`se$ConditionNone`)
```

```
## Coefficient:  se$ConditionNone 
##            logFC   logCPM        F       PValue       FDR
## Egr4  -1.4403405 5.885444 53.26253 7.366511e-05 0.6397289
## Nr4a1 -1.3308764 6.692671 49.36148 9.681327e-05 0.6397289
## Fos   -2.1307567 5.100641 45.14032 1.330274e-04 0.6397289
## Egr2  -2.6011381 2.423339 36.51426 2.782330e-04 0.8744991
## Junb  -1.0440860 6.806899 35.61133 3.030773e-04 0.8744991
## Fosb  -1.4051175 4.423401 27.72609 6.981470e-04 0.9843760
## Sik1  -1.1741066 3.640998 26.05440 8.536970e-04 0.9843760
## Arc   -1.5191907 7.690720 20.90101 1.703817e-03 0.9843760
## Trib1 -1.0700034 3.599083 19.38456 2.139942e-03 0.9843760
## Dusp5 -0.9521124 5.808549 18.05644 2.641813e-03 0.9843760
```
Indeed, there are two candidate genes that pass the multiple testing correction, Fos and Egr4


Are there any genes altered by early life stress?

```r
topTags(Results$`se$ELSNone`)
```

```
## Coefficient:  se$ELSNone 
##              logFC     logCPM        F       PValue       FDR
## Plekhg3 -1.9967921  2.4615140 32.28299 0.0004222831 0.9996194
## Itgb8   -0.9937680  3.5789149 19.23668 0.0021895669 0.9996194
## Trpc6   -1.7023306  2.7981601 16.34132 0.0035283596 0.9996194
## Igfbp5  -0.9692028  5.7252455 15.51154 0.0040903459 0.9996194
## Sft2d1  -1.3952806  2.0899822 14.30578 0.0051219008 0.9996194
## Galnt4  -1.6777129 -0.2792201 14.24235 0.0051846914 0.9996194
## Grm2    -0.7297558  4.5540417 13.93444 0.0055035518 0.9996194
## Ctdspl2 -1.0900578  3.7103135 13.35386 0.0061743225 0.9996194
## Fgfr3   -1.0389808  3.9430810 12.14425 0.0079339010 0.9996194
## Gabrg1  -1.2811709  1.8892472 11.61273 0.0089039549 0.9996194
```
Unfortunately, there are no genes that pass the multiple testing correction

Let's investigate if there are genes with a significant interaction

```r
topTags(Results$`se$ConditionNone:se$ELSNone`)
```

```
## Coefficient:  se$ConditionNone:se$ELSNone 
##              logFC     logCPM        F      PValue       FDR
## Plekhg3  2.8250275  2.4615140 36.30672 0.000283712 0.9998581
## Piga     2.6242173  1.5385751 22.16459 0.001422287 0.9998581
## Nadsyn1 -8.1411065 -0.1856326 21.32271 0.002252620 0.9998581
## Stard9   5.7100942 -0.4341469 14.37339 0.005056005 0.9998581
## Mill2   -4.5346623 -1.3626955 14.29680 0.005130735 0.9998581
## Ctdspl2  1.3933425  3.7103135 11.61326 0.008902932 0.9998581
## Il10ra   8.2935183 -1.1359543 13.82403 0.009279340 0.9998581
## Gpr17    2.1180496  2.8256340 11.32226 0.009496970 0.9998581
## B3gat2  -0.9225619  4.8564532 11.24795 0.009656549 0.9998581
## Nnt     -1.7460846  3.4884480 10.80953 0.010669564 0.9998581
```
Unfortunately, no genes have a altered acute stress response in ELS vs normal animals



## Meta-Analysis of all data

Let's run an overarching analysis over all data to determine if there are any significant effects for acute stress (=Condition), Genotype, GEOdataset or early life stress

```r
se <- kallistodata

#experimental design, full additive model
design <- model.matrix(~se$Condition + se$Gender + se$Geontype + se$Set + se$ELS)

y <- DGEList(counts=assays(se)$counts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)

#filter out genes that are below 10 counts in more than 75% of samples
keep <- rowSums(y$counts>10) >= 5
y <- y[keep, , keep.lib.sizes=FALSE]

Results <- list()
fit <- glmQLFit(y,design)
for(i in colnames(design)[-1]){
  Results[[i]] <- glmQLFTest(fit, i)
}
```

Les's investigate if there are any genes altered by acute stress

```r
topTags(Results$`se$ConditionNone`)
```

```
## Coefficient:  se$ConditionNone 
##            logFC   logCPM         F       PValue          FDR
## Egr4  -1.4504092 5.700383 214.14618 5.534782e-11 7.942412e-07
## Fos   -2.1211555 4.906998 149.96297 8.624405e-10 6.188011e-06
## Egr2  -2.3426004 2.247529 134.08514 2.007587e-09 9.602956e-06
## Fosb  -1.4963321 4.223948 110.24820 8.589276e-09 3.081403e-05
## Nr4a1 -1.1486927 6.534812  87.33636 4.627058e-08 1.327966e-04
## Sik1  -1.1092486 3.493908  77.98207 1.028066e-07 2.458792e-04
## Junb  -0.9520758 6.615866  67.02760 2.918232e-07 5.982375e-04
## Dusp5 -1.1245332 5.633421  57.96897 7.723457e-07 1.385395e-03
## Arc   -1.3361997 7.613828  49.21993 2.232716e-06 3.559942e-03
## Egr1  -0.8278663 8.097496  34.85452 1.828138e-05 2.623379e-02
```
Indeed, there are multiple candidate genes that are significantly altered by acute stress across other conditions

Are there any genes altered by gender

```r
topTags(Results$`se$Gendermale`)
```

```
## Coefficient:  se$Gendermale 
##              logFC   logCPM          F       PValue          FDR
## Eif2s3y 11.3972639 3.455656 1281.93513 2.172442e-13 3.117454e-09
## Uty     10.3866137 2.653583  432.75005 1.184239e-10 8.496914e-07
## Kdm5d    2.8606312 1.942235  128.78858 2.715137e-09 1.298741e-05
## Ddx3y   11.8469936 4.146537  125.85803 1.215087e-07 4.359125e-04
## mt-Nd5  -3.0704215 5.319974   47.54160 2.782169e-06 7.984825e-03
## Eif2s3x -0.8739168 7.331434   39.41080 8.845076e-06 2.115447e-02
## Pou3f1   1.2923702 6.737009   35.91281 1.535610e-05 3.148001e-02
## Midn     0.6843066 5.651226   33.51025 2.293767e-05 4.114445e-02
## mt-Nd6  -3.3953597 4.653982   32.61893 2.675486e-05 4.265914e-02
## Spry4    0.7043658 4.093567   28.19129 6.014896e-05 8.631376e-02
```
Indeed, there are multiple candidate genes that are significantly altered by gender across other conditions


Are there any genes altered by BDNF Val66Met phenotype?

```r
topTags(Results$`se$GeontypeWild Type`)
```

```
## Coefficient:  se$GeontypeWild Type 
##                   logFC     logCPM        F       PValue       FDR
## Bloc1s6      -1.6845415  3.6033657 31.83990 3.068073e-05 0.2304958
## Cd59a         1.7142445  1.9530726 31.58108 3.212485e-05 0.2304958
## Tcp11l1      -1.4606140  2.8547082 23.37359 1.606079e-04 0.4938836
## Luzp1        -0.6417807  7.1666383 18.10729 5.481632e-04 0.4938836
## Prl          -7.0088172  3.3317522 18.10258 5.488134e-04 0.4938836
## Alg3         -1.3060785  2.0287865 17.24838 6.823469e-04 0.4938836
## RP23-78D19.4  3.6049972 -0.1249779 16.74718 7.775085e-04 0.4938836
## Strc          1.3320817  0.4421829 16.53116 8.230530e-04 0.4938836
## Fyttd1       -0.5620122  6.8522505 16.01712 9.439905e-04 0.4938836
## Reps2        -0.5554824  8.6385483 15.48113 1.091807e-03 0.4938836
```
Unfortunately, there are no genes that pass the multiple testing correction

Are there any genes altered early life stress?

```r
topTags(Results$`se$ELSNone`)
```

```
## Coefficient:  se$ELSNone 
##                   logFC     logCPM         F      PValue       FDR
## RP23-445H7.1 -0.7291894 2.74282804 15.428622 0.001107636 0.9997824
## Scel         -0.8256697 0.06239698 12.354656 0.002700464 0.9997824
## Rrp15         0.7312827 3.30166498 11.562444 0.003456599 0.9997824
## Gm2163       -1.6578356 0.62733604  9.331060 0.007248642 0.9997824
## Map2k5        0.5565325 5.51777906  8.849308 0.008589268 0.9997824
## Psip1        -0.4244342 8.26355564  8.842136 0.008611244 0.9997824
## Serbp1       -0.4102638 8.22281402  8.387176 0.010144829 0.9997824
## Atpaf2        0.4073688 4.44786602  8.224273 0.010767483 0.9997824
## Rsl1d1       -0.3888628 6.32758454  8.186028 0.010919874 0.9997824
## Wipf3         0.4153641 9.72897807  7.749044 0.012846573 0.9997824
```
Unfortunately, there are no genes that pass the multiple testing correction
