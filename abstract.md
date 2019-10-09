<!--
 In the context of biomedical research with human cell lines, the ability to extract the expression levels of individual genes from composite samples is of paramount importance.
 -->

Deconvolution problems are ubiquitous and have a long history in many areas of science and engineering. Such methods are of paramount importance for applications in which individual sensors are used to measure multiple distinct signals. The L1000 assay, which uses 500 sensors to measure the expression levels of 1,000 unique genes, offers a concrete example. The validation of existing deconvolution methods is often difficult because of the lack of ground-truth datasets, as well as the costs involved in the development of sound alternatives. 

To address the problem, we used L1000, to generate a novel dataset of the response of 1,000 genes to various perturbagens (shRNA and compounds) where genes were measured both in tandem and individually. We then ran an open competition to explore different approaches. The competition produced several approaches to the problem based on a wide variety of techniques including cutting edge machine-learning algorithms, such as Convolutional Neural Networks and Random Forests, and more traditional approaches, such as gaussian mixtures and k-means. All top solutions achieved notable improvements over the benchmark. 

