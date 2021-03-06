# cluster_diffex

This is a pipeline for identifying clusters of similar cells in single cell RNA-Seq data and conducting pairwise statistical analysis of differential expression between the single cell profiles in those clusters. The program is designed to take similarity matrices computing using the similarity_pipeline as input. This similarity matrix is used to compute a k-nearest neighbors graph and identify clusters using the Phenograph algorithm (from Jacob Levine and Dana Pe'er). The pipeline then uses the SCDE algorithm (from Peter Kharchenko) for differential expression analysis. SCDE models the dispersion of RNA-Seq counting data using the negative binomial distribution and the dropout using the Poisson distribution. The algorithm uses Bayesian statistics to combine these two pieces of information and infer differential expression. We use SCDE to conduct differential expresion analysis between all possible pairs of clusters found using Phenograph. We also superimpose the expression of the differentially expressed genes on t-SNE visualizations of the cells (as output by the similarity_pipeline).

# Requirements

Python3.4 or higher                                                                                                                         
Numpy                                                                                                                                           
Matplotlib Python module                                                                                                                      
Phenograph version 1.5.2 Python module (from https://github.com/jacoblevine/PhenoGraph)                                                      
R version 3.2 or higher                                                                                                                                                                                                                                               
SCDE version 1.2.1 and its dependencies (from http://hms-dbmi.github.io/scde/package.html)                      
