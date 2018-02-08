# TD-multiOmics
Tensor Based unsupervised feature extraction was applied to multiomics data which was used for performance evaluation of DAIBLO algorithm implemented in mixOmics pacage

# Introduction

Recently, I poposed multi-view data analysis strategy using tensor based unsupervised feature extraction ([paper](https://www.ncbi.nlm.nih.gov/pubmed/28841719)), and applied it to multi-omics data among them. However, because there is no multi-omics in the title, I feel that I am pretty defeated by the appealing power to the [paper](https://doi.org/10.1371/journal.pcbi.1005752) of the recent [mixOmics package](http://mixomics.org/) which emphasized this point. Therefore, I would like to compare the performance when applying the variable selection multi-view data application by unsupervised learning using tensor decomposition to the dataset whose performance is tried with mixOmics as well as the propagation of own method.

# About mixOmics

As its name suggests, [mixOmcs](http://mixomics.org/) is an R package specialized for analyzing multi-omics data, [published](https://cran.r-project.org/web/packages/mixOmics/index.html) in CRAN, and can be normally installed using the install.packages () command. Various methods are included, but here we compare the performance with the method named [DIABLO](http://mixomics.org/mixdiablo/) that they regard as the most advanced. Figure 1 in [preprint](https://www.biorxiv.org/content/early/2016/08/03/067611) is the easiest to understand in DIABLO, but you can create a variable by multiplying multiple multi-omics data in pair wise and then summing up the variables, and using that variable to set the predictor (Supervised learning) method to construct multiple classes. A proper explanation using mathematical expressions is posted in detail in the Sec. 4 in [text S1](https://doi.org/10.1371/journal.pcbi.1005752.s001) of the mixOmics [paper](https://www.ncbi.nlm.nih.gov/pubmed/29099853). Since it is not the purpose of this paper to enter the detailed explanation of DIABLO method, please refer to the papers linking here for details.


# DIABLO: A case study

The page of the above-mentioned [preprint](https://www.biorxiv.org/content/early/2016/08/03/067611) or [DIABLO execution example](http://mixomics.org/mixmint/stemcells-example/) describes an example of integrated analysis of mRNA, miRNA, proteomics data using TCGA data.

This data is displayed on the page of [DIABLO execution example](http://mixomics.org/mixmint/stemcells-example/) as

>\#\# \$mRNA  
>\#\# [1] 150 200  
>\#\#   
>\#\# \$miRNA  
>\#\# [1] 150 184  
>\#\#   
>\#\# \$proteomics  
>\#\# [1] 150 142  

these say that there are 150 samples with 200 mRNA, 184 miRNA and 142 proteomics measurements. 150 samples consist of three classes,

>\#\# Basal  Her2  LumA   
>\#\#    45    30    75

as shown in  [DIABLO execution example](http://mixomics.org/mixmint/stemcells-example/).
