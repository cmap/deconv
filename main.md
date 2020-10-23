---
title: "Improving Deconvolution Methods in Biology through Open Innovation Competitions: \
        An Application to the Connectivity Map"

abstract: "Do machine learning methods improve standard deconvolution techniques for gene expression data? This paper uses a unique new dataset combined with an open innovation competition to evaluate a wide range of gene-expression deconvolution approaches developed by 294 competitors from 20 countries. The objective of the competition was to separate the expression of individual genes from composite measures of gene pairs. Outcomes were evaluated using direct measurements of single genes from the same samples. Results indicate that the winning algorithm based on random forest regression outperformed the other methods in terms of accuracy and reproducibility. More traditional gaussian-mixture methods performed well and tended to be faster. The best deep learning approach yielded outcomes slightly inferior to the above methods. We anticipate researchers in the field will find the dataset and algorithms developed in this study to be a powerful research tool for benchmarking their deconvolution methods and a useful resource for multiple applications."
    
date: 'Last updated: \today'

keywords:
  - connectivity map
  - crowdsourcing
  - deconvolution algorithm
  - gene expression
  - machine learning
  - open innovation competition

geometry: margin=1in
bibliography: dpeak.bib
biblio-title: References
biblio-style: unsrtnat
natbiboptions: numbers,super,compress
linestretch: 1.5
toc: false 

header-includes:
  - \usepackage{bbm}
  - \usepackage{longtable,booktabs}
  - \newcommand\Fig[1]{\nameref{#1}}
  - \usepackage{authblk}
  - \author[1,2,3]{Andrea Blasco\thanks{ablasco@hbs.edu}}
  - \author[3]{Ted Natoli\thanks{tnatoli@broadinstitute.org}}
  - \author[1,2]{Michael G. Endres}
  - \author[1,2]{Rinat A. Sergeev}
  - \author[1,2]{Steven Randazzo}
  - \author[1,2]{Jin H. Paik}
  - \author[3]{N. J. Maximilian Macaluso}
  - \author[3]{Rajiv Narayan}
  - \author[3]{Xiaodong Lu}
  - \author[3]{David Peck}
  - \author[1,2,4]{Karim R. Lakhani}
  - \author[3]{Aravind Subramanian}
  - \affil[1]{Laboratory for Innovation Science at Harvard, Cambridge, MA, USA}
  - \affil[2]{Harvard Business School, Harvard University, Boston, MA, USA}
  - \affil[3]{Broad Institute of Harvard and MIT, Cambridge, MA, USA}
  - \affil[4]{National Bureau of Economic Research, Cambridge, MA, USA}
  #- \usepackage{titlesec}
  #- \titleformat{\subsection}[runin]{\normalfont\bfseries}{\thesection.\quad}{.5em}{}
  
---
