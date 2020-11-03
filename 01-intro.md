Introduction
============

A recurring problem in biomedical research is how to isolate signals of distinct populations (cell types, tissues, and genes) from composite measures obtained by a single analyte or sensor. This deconvolution problem often stems from the prohibitive cost of profiling each population separately [@subramanian2017next; @cleary2017efficient] and has important implications for the analysis of transcriptional data in mixed samples [@shen2010cell; @zhong2012gene; @newman2015robust; @zaitsev2019complete], single-cell data [@deng2019scalable], the study of cell dynamics [@lu2003expression], and imaging data [@preibisch2014efficient]. 

<!-- Deconvolution methods have a long history with many different applications.  -->

In the context of gene expression analysis, available deconvolution approaches offer several advantages but also have various limitations [@shen2010cell; @shen2013computational]. Machine learning approaches present a promising route towards improvement but implementing and benchmarking these methods can be challenging. 
A prominent issue is the scarcity of realistic "ground truth" data for training and validation. Existing methods are often trained, and have their results validated, on synthetic data rather than biologically mixed expression data, which are more difficult to generate and evaluate. But even when ground truth data is available, other difficulties arise that may prevent an effective evaluation of machine learning approaches, such as complex parameter optimization which requires substantial experience in more advanced machine learning techniques.

In this paper, we addressed these challenges using an open innovation competition. The goal of the competition was to solve a gene-expression deconvolution problem with application to the Connectivity Map (CMap), a catalog of over 1.3 million human gene-expression profiles [@subramanian2017next]. To reach such a massive scale, CMap focuses on a reduced representation of the transcriptome, consisting of approximately 1,000 human genes, called "landmarks." These genes were selected to capture a large portion of the cell’s transcriptional state, thus yielding significant cost reductions compared to traditional methods (RNA-sequencing). However, the data-generation technology (a bead-based multiplex assay called L1000)  is limited to 500 different bead colors per sample. Hence, to measure 1,000 genes per sample, CMap coupled each of the available bead types with a gene pair, using a k-means algorithm called "dpeak" to separate the expression of 1,000 genes from the 500 bead measurements. 

The goal of the competition was to develop methods to improve upon the dpeak algorithm. Examples from past competitions have shown that machine learning algorithms developed through open competitions tend to outperform solutions derived through conventional means [@lakhani2013prize; @good2013crowdsourcing; @blasco2019advancing]. Accordingly, we wanted to see if machine learning solutions developed through the contest could yield considerable improvements over the current dpeak algorithm. 
 
We generated a novel dataset for the contest, which is available online (\nameref{s1-availability-and-implementation}). It consists of a collection of gene-expression profiles for 122 different perturbagens, both short hairpin RNA (shRNA) and compound treatments at multiple replicates, for a total of over 2,200 gene expression experiments. The same set of experiments was profiled twice while varying the detection mode for acquiring the data. Detection varied between a dual procedure (DUO), two genes per bead barcode, and a uni procedure (UNI), one gene per bead barcode. The UNI data served as the ground truth against which deconvolution procedures applied to the DUO data were compared.
 
The contest lasted 21 days and was run on Topcoder (Wipro, India), a popular crowdsourcing platform. Competitors had access to both UNI and DUO datasets, which could be used to inform the development of their deconvolution solutions. Performance evaluation was based on pre-specified metrics of accuracy and computational speed.  A prize purse of $23,000 in cash was offered to competitors as an incentive to be divided among the top 9 submissions. To align incentives and prevent problems like overfitting, the cash prizes were awarded based on the performance evaluation on the _holdout_ (not seen by the competitors during the competition) subset of the data.

<!-- 
The contest drew about 300 competitors from 20 different countries. The top solutions included machine-learning methods, such as Random Forests and Convolutional Neural Networks (CNNs), as well as more traditional models such as Gaussian mixtures and k-means algorithms. These approaches performed significantly better than the L1000 benchmark in various measures of accuracy and computational speed, and likely have application beyond gene expression.
 -->

<!-- 
TO READ: 

"DeconvSeq: deconvolution of cell mixture distribution in sequencing data" https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz444/5506629?redirectedFrom=fulltext

"Evaluation of methods to assign cell type labels to cell clusters from single-cell RNA-sequencing data [version 3; peer review: 2 approved, 1 approved with reservations]"
https://f1000research.com/articles/8-296

Systematic comparative analysis of single cell RNA-sequencing methods
https://www.biorxiv.org/content/10.1101/632216v1

The Technology and Biology of Single-Cell RNA Sequencing
https://www.sciencedirect.com/science/article/pii/S1097276515002610

Fractional proliferation: a method to deconvolve cell population dynamics from single-cell data
https://www.nature.com/articles/nmeth.2138

Bulk tissue cell type deconvolution with multi-subject single-cell expression reference
https://www.nature.com/articles/s41467-018-08023-x

A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure
https://www.sciencedirect.com/science/article/pii/S2405471216302666

-->