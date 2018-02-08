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

![図１](https://i.imgur.com/xX7SJXB.png,"図１")

The figure above shows the classification performance by how many components created with DIABLO are used. The vertical axis is the error rate, but sufficient performance is obtained with the first two components. Also, the embedding of all 150 samples in this space

![図２](https://i.imgur.com/rxzWhea.png,"図２")

can be shown in the above.　It is rather obvious that thee sub classes are well separated on this plane. In addition to this, DIABLO has a function that select important features. 

![図３](https://i.imgur.com/X6uT5I6.png,"図３")


The above figure is the heatmap of the variable selected by DIABLO. Although rows are samples and columns are (selected) omics data, it is possible to determine three classes sufficiently by hierarchical clustering alone without having to carry out multivariate analysis each time variable selection can be performed on omics data. 

# Applying multi-view data analysis using tensor decomposition based unsupervised feature extraction to multi-omics data

Now, let us use the application of multi-view data analysis of variable selection by unsupervised learning using tensor decomposition to the above data set. In my [paper] (https://www.ncbi.nlm.nih.gov/pubmed/28841719), it corresponds to Type I Case I. Here, <img src="https://latex.codecogs.com/gif.latex?x_{i_1,j}&space;\in&space;\mathbb{R}^{200&space;\times&space;150}" title="x_{i_1,j} \in \mathbb{R}^{200 \times 150}" /> expresses the expression level in <img src="https://latex.codecogs.com/gif.latex?j" title="j" /> th sample of <img src="https://latex.codecogs.com/gif.latex?i_1" title="i_1" /> th mRNA, <img src="https://latex.codecogs.com/gif.latex?x_{i_1,j}&space;\in&space;\mathbb{R}^{184&space;\times&space;150}" title="x_{i_2,j} \in \mathbb{R}^{184 \times 150}" /> represents the expression level in the <img src="https://latex.codecogs.com/gif.latex?j" title="j" /> th sample of the <img src="https://latex.codecogs.com/gif.latex?i_2" title="i_2" /> th miRNA, <img src="https://latex.codecogs.com/gif.latex?x_{i_3,j}&space;\in&space;\mathbb{R}^{142&space;\times&space;150}" title="x_{i_3,j} \in \mathbb{R}^{142 \times 150}" /> represents the expression level in the <img src="https://latex.codecogs.com/gif.latex?j" title="j" /> th sample of <img src="https://latex.codecogs.com/gif.latex?i_3" title="i_3" /> th proteomics, then the tensor <img src="https://latex.codecogs.com/gif.latex?x_{i_1,i_2,i_3,j}" title="x_{i_1,i_2,i_3,k}" /> is defined as

<img src="https://latex.codecogs.com/gif.latex?x_{i_1i_2i_3j}&space;\equiv&space;x_{i_1,j}&space;\cdot&space;x_{i_2,j}&space;\cdot&space;x_{i_3,j}" title="x_{i_1i_2i_3j} \equiv x_{i_1,j} \cdot x_{i_2,j} \cdot x_{i_3,j}" />

After applying [HOSVD](https://en.wikipedia.org/wiki/Higher-order_singular_value_decomposition) to <img src="https://latex.codecogs.com/gif.latex?x_{i_1,i_2,i_3,j}" title="x_{i_1,i_2,i_3,k}" />, we get

<img src="https://latex.codecogs.com/gif.latex?x_{i_1i_2i_3j}&space;=&space;\sum_{\ell_1,\ell_2,\ell_3,\ell_4}&space;G(\ell_1,\ell_2,\ell_3,\ell_4)&space;x_{\ell_1i_1}x_{\ell_2i_2}x_{\ell_3i_3}x_{\ell_4j}" title="x_{i_1i_2i_3j} = \sum_{\ell_1,\ell_2,\ell_3,\ell_4} G(\ell_1,\ell_2,\ell_3,\ell_4) x_{\ell_1i_1}x_{\ell_2i_2}x_{\ell_3i_3}x_{\ell_4j}" />

![plot.jpg](https://qiita-image-store.s3.amazonaws.com/0/199087/333c0a15-0509-dc34-4016-7168de3cfdca.jpeg)

The above figure shows the embedding of 150 samples coposed of three sub classes into plane spaned by <img src="https://latex.codecogs.com/gif.latex?x_{\ell_1,j}" title="x_{\ell_1,k}" />, <img src="https://latex.codecogs.com/gif.latex?\ell_4=1,4" title="\ell_4=1,4" />. It is also obvious that 150 samples are well separated into three sab classes. Confusion table (rows: prediction, columns: true classes) obtained by linear discriminat analysis using these two components is

|       |Basal |  Her2 |LumA |
|:-----|----:|-----:|----:|
| **Basal** |  **42** |     4 |    0|
| **Her2**   |    2 |  **25**   |  2  |
|  **LumA** |    1 |    1  | **73**  |

The accuracy is as large as 0.94. 
