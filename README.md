# TD-multiOmics
Tensor Based unsupervised feature extraction was applied to multiomics data which was used for performance evaluation of [DIABLO](http://mixomics.org/mixdiablo/) algorithm implemented in [mixOmics package](http://mixomics.org/)

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

can be shown in the above.　It is rather obvious that three sub classes are well separated on this plane. In addition to this, DIABLO has a function that select important features. 

![図３](https://i.imgur.com/X6uT5I6.png,"図３")


The above figure is the heatmap of the variable selected by DIABLO. Although rows are samples and columns are (selected) omics data, it is possible to determine three classes sufficiently by hierarchical clustering alone without having to carry out multivariate analysis each time variable selection can be performed on omics data. 

# Applying multi-view data analysis using tensor decomposition based unsupervised feature extraction to multi-omics data

Now, let us use the application of multi-view data analysis of variable selection by unsupervised learning using tensor decomposition to the above data set. In my [paper](https://www.ncbi.nlm.nih.gov/pubmed/28841719), it corresponds to Type I Case I. Here, <img src="https://latex.codecogs.com/gif.latex?x_{i_1,j}&space;\in&space;\mathbb{R}^{200&space;\times&space;150}" title="x_{i_1,j} \in \mathbb{R}^{200 \times 150}" /> expresses the expression level in <img src="https://latex.codecogs.com/gif.latex?j" title="j" /> th sample of <img src="https://latex.codecogs.com/gif.latex?i_1" title="i_1" /> th mRNA, <img src="https://latex.codecogs.com/gif.latex?x_{i_2,j}&space;\in&space;\mathbb{R}^{184&space;\times&space;150}" title="x_{i_2,j} \in \mathbb{R}^{184 \times 150}" /> represents the expression level in the <img src="https://latex.codecogs.com/gif.latex?j" title="j" /> th sample of the <img src="https://latex.codecogs.com/gif.latex?i_2" title="i_2" /> th miRNA, <img src="https://latex.codecogs.com/gif.latex?x_{i_3,j}&space;\in&space;\mathbb{R}^{142&space;\times&space;150}" title="x_{i_3,j} \in \mathbb{R}^{142 \times 150}" /> represents the expression level in the <img src="https://latex.codecogs.com/gif.latex?j" title="j" /> th sample of <img src="https://latex.codecogs.com/gif.latex?i_3" title="i_3" /> th proteomics, then the tensor <img src="https://latex.codecogs.com/gif.latex?x_{i_1,i_2,i_3,j}" title="x_{i_1,i_2,i_3,k}" /> is defined as

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

The accuracy is as large as 0.94. Feature selection can be done using these two selected components.
 <img src="https://latex.codecogs.com/gif.latex?G(\ell_1,\ell_2,\ell_3,\ell_4)" title="G(\ell_1,\ell_2,\ell_3,\ell_4)" /> are sorted by thier absolute values for <img src="https://latex.codecogs.com/gif.latex?\ell_4=1,4" title="\ell_4=1,4" />, <img src="https://latex.codecogs.com/gif.latex?\ell_1,\ell_2,\ell_3" title="\ell_1,\ell_2,\ell_3" /> with  <img src="https://latex.codecogs.com/gif.latex?G" title="G" /> having larger absolute values are selected.
 
 | rank | <img src="https://latex.codecogs.com/gif.latex?\ell_1" title="\ell_1" /> | <img src="https://latex.codecogs.com/gif.latex?\ell_2" title="\ell_2" />  | <img src="https://latex.codecogs.com/gif.latex?\ell_3" title="\ell_3" /> |<img src="https://latex.codecogs.com/gif.latex?\ell_4" title="\ell_4" />  | <img src="https://latex.codecogs.com/gif.latex?G(\ell_1,\ell_2,\ell_3,\ell_4)" title="G(\ell_1,\ell_2,\ell_3,\ell_4)" /> |
|:---:|:---:|:---:|:---:|:---:|:--:|
|1   |  1 | 1 | 1 | 1 |   -407857.582 |
|2   |  1 | 1 | 4 | 4 |    -209720.615 |
|3   |  2 | 1 | 1  |4 |    -20452.480 |
|4   |  2 | 1 | 3 |1  |  -11677.505 |
|5   |  2 | 1 | 4 |1  |  -10428.742 |
|6   |  2 | 1 | 2 |1  |  10157.467 |
|7   |  1 | 1 | 2 | 1 |  -8973.774 |
|8   |  1 | 2 | 1 | 4 |   8360.976 |
|9   |  2 | 1 | 5 | 4 |   -6628.467 |
|10  |  1 | 1 | 3 | 4 |   6623.046 |

The above table that lists top 10 ranked ones, <img src="https://latex.codecogs.com/gif.latex?\ell_1=1,2,\ell_2=1,2,\ell_3=1,2,3,4" title="\ell_1=1,2,\ell_2=1,2,\ell_3=1,2,3,4" /> turn out to be selected. With assuming that <img src="https://latex.codecogs.com/gif.latex?x_{\ell_1i_1},x_{\ell_2i_2},x_{\ell_3i_3}" tile="x_{\ell_1i_1},x_{\ell_2i_2},x_{\ell_3i_3}"> obey multiple Gaussian distribution as null hypothesis、P-values are attributed to <img src="https://latex.codecogs.com/gif.latex?i_1,i_2,i_3" title="i_1,i_2,i_3">.  The top 10 <img src="https://latex.codecogs.com/gif.latex?i_1,i_2,i_3" title="i_1,i_2,i_3"> are selected respectively. 

![heatmap.jpg](https://qiita-image-store.s3.amazonaws.com/0/199087/e1deb27d-3a3a-4184-eab8-2538aa482e4d.jpeg)

Rows are samples(black:Basel, red:Her2, green:Luma), columns are omics (blue:mRNA,pink:miRNA,cyan:proteomics).
Comapative to DIABLO, feature that allows hierarchical clautsering classify three subclasses are well selected.


# Discussion

[DIABLO](http://mixomics.org/mixdiablo/) is a fairly complicated calculation, you must create a design matrix on how to combine multi-omics data yourself, and use label information There is supervised learning. On the other hand, tensor decomposition is unsupervised learning, and what you are doing is simple. Actually, although I raised the execution code of R as tensor.R in this repository, it is simple enough to beat it. I feel a little unbelievable if I think that the same performance as DIABLO is done with this. However, DIABLO also thinks of the product between multi-omics data, and thinks that if it performs linear discrimination it may mean that the direction is not deviated so much in the meaning.

In the tensor decomposition, the tensor becomes huge (in the present case, there are 150 × 200 × 184 × 142 elements), so there is a drawback that calculation time is required. Actually, although DIABLO can be executed by a lap-top, tensor decomposition can not be performed without a server machine with dozens of giga of memory. Nevertheless, in the future, tensor decomposition will be frequently used for multi-omics data analysis.
