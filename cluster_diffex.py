#! /usr/bin/python3.5
import sys
import os

########### Peter Sims, Columbia University, 2016 ###########

# This is the main program in a pipeline for identifying clusters of similar cells in single cell RNA-Seq data and 
# pairwise statistical analysis of differential expression between those clusters. The program is designed to take
# similarity matrices computed using the similarity_pipeline (DropoutDE_analysis5.py is main pipeline program) as
# input. This similarity matrix is used to compute a k-nearest neighbors graph and identify clusters using the
# Phenograph algorithm (from Jacob Levine and Dana Pe'er). The pipeline then uses the SCDE algorithm (from Peter
# Kharchenko) for differential expresison analysis. SCDE models the dispersion of RNA-Seq counting data using the
# negative binomial distribution and dropout using the Poisson distribution. The algorithm uses a Bayesian method
# to combine these two pieces of information and infer differential expression.  We use SCDE to conduct differential
# expression analysis between all possible pairs of clusters found using Phenograph.  We also superimpose the expression
# of the differentially expressed genes on t-SNE visualizations of the cells (as output by the similarity_pipeline).

run_NAME = sys.argv[1] # name of sequencing sample, e.g. PJ018
matrix_INFILE = sys.argv[2] # name of counts matrix file - typically, col 1 is gene id, col 2 is gene symbol, col 3 contains molecular counts of first cell
norm_INFILE = sys.argv[3] # three-column file - col 1 is cell-identifying barcode, col 2 is number of molecules detected, col 3 is number of genes detected 
corr_INFILE = sys.argv[4] # file containing similarity matrix
tsne_INFILE = sys.argv[5] # file containing t-SNE coordinates (or coordinates for any other two-dimensional representation of the cells)
proc_N = int(sys.argv[6]) # number of threads

# Compute k-nearest neighbors graph and identify clusters using Phenograph and output new counts matrices for each cluster.
# The program also outputs a t-SNE visualization of the clusters.
print('Running Phenograph on correlation matrix...')
pg_OUTFILE = run_NAME+'/'+run_NAME+'.pg.txt'
pg_PREFIX = run_NAME+'/'+run_NAME+'.pg'
tsne_PDF = run_NAME+'/'+run_NAME+'.tsne5pg.pdf'
cmd = 'python3.5 phenograph_tsne.py %(corr_INFILE)s %(tsne_INFILE)s %(matrix_INFILE)s %(pg_OUTFILE)s %(pg_PREFIX)s %(tsne_PDF)s' % vars()
print(cmd)
os.system(cmd)


print('Running SCDE on pairwise communities...')
# Get the community (cluster) number of each cell.
communities = []
with open(pg_OUTFILE) as f:
	for line in f:
		communities.append(int(line.split()[0]))

mn = min(communities)
mx = max(communities)
N = len(set(communities))
print('Found %(N)d communities...' % vars())
# For each pairwise combination of clusters, run SCDE to assess differential expression and color t-SNE plots with differentially expressed genes.
for i in range(mn,mx+1):
	if len(str(i)) == 1: # some stuff for managing file names and headers
		c1 = 'comm0'+str(i)
	else:
		c1 = 'comm'+str(i)
	for j in range(i+1,mx+1):
		print('Running SCDE on %(i)d vs. %(j)d...' % vars())
		if len(str(j)) == 1:
			c2 = 'comm0'+str(j)
		else:
			c2 = 'comm'+str(j)
		cmd = 'python3.5 scde_diffex.py %(pg_PREFIX)s_%(i)d.matrix.txt %(pg_PREFIX)s_%(j)d.matrix.txt %(run_NAME)s.%(c1)s_%(c2)s.diffex.txt %(run_NAME)s_%(c1)s_%(c2)s.R %(run_NAME)s/%(run_NAME)s_%(c1)s_%(c2)s %(c1)s %(c2)s gene_names.txt %(proc_N)d' % vars()
		print(cmd)
		os.system(cmd)
		print('t-SNE output for SCDE on %(i)d vs. %(j)d...' % vars())
		cmd = 'python3.5 diffex_tsne.py %(run_NAME)s %(run_NAME)s/%(run_NAME)s_%(c1)s_%(c2)s/%(run_NAME)s.%(c1)s_%(c2)s.diffex.txt.conv.txt %(matrix_INFILE)s %(norm_INFILE)s %(tsne_INFILE)s 3 0.0001' % vars()
		print(cmd)
		os.system(cmd)
		

