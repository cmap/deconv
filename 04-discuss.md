
Discussion 
===========

<!-- Our data contain several replicate samples for each perturbagen (4 for shRNAs and 10 for compounds). Minimizing inter-replicate variation is of crucial practical importance to biologists because it translates to higher reproducibility of results. In this regard, the benchmark k-means solution is likely suboptimal because it does little to mitigate the discrepancy in variability between the genes measured with high and low bead proportion. MOVE TO DISCUSSION -->

We created a novel dataset of L1000 profiles for over 120 shRNA and compound experiments with several replicates for a total of 2,200 gene expression profiles of genes measured independently, and in tandem. This dataset constitutes now a public resource (\nameref{s1-availability-and-implementation}) to all the researchers in this area who are interested in testing their deconvolution approaches. 

Using an open innovation competition, we collected and evaluated multiple and diverse deconvolution methods. 
The best approach was based on a random forest,  which is a collection of decision tree regressors. This method achieved: (i) the highest global correlation between the ground-truth and the corresponding deconvoluted data, (ii) the lowest inter-replicate variation, and (iii), compared to the benchmark, was able to detect more than a thousand additional extremely modulated genes, while reducing the false positives at the same time. Our analysis further showed that these gains are considerable when the gene populations were sampled in different proportions (here, genes in high and low bead proportions), with the k-means benchmark approach being systematically less accurate because it does little to mitigate the discrepancy in variability between the genes measured with high and low bead proportion. 

In addition, the random-forest approach achieved  these improvements with only 10 trees on 60 features. Thus, the algorithm is also relatively fast and easy to implement. By comparison, the fastest approach used a more traditional Gaussian mixture model (with plate-level adjustments), which turned out to be less accurate. Hence, and overall, our analysis provided evidence of the tremendous potential of using random-forest methods for deconvolution problems in biology.

<!-- 
Summary of the results presented in the methods section. 
Discussion generality of the solutions
- Novel? Have any of these solutions previously been applied to deconvolution problems?
- Specific to this problem or general to others?
Discuss implications of these methods for CMap production
- Preliminary results on past data conversion
- Directions for pipeline integration and generation of future data
- Cost savings
- Implementation strategy and outcomes
- Increase in data processing throughput
-->
