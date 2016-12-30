#! /usr/bin/python3.4
import sys
import os

# This program is part of the cluster_diffex pipeline for single cell RNA-Seq differential expession analysis.  The purpose of this program
# is to generate t-SNE visualization plots color-coded by the differentially expressed genes found in each pairwise comparison of single cell
# clusters found by the pipeline.  The t-SNE coordinates generally come from the similarity_pipeline.

run_NAME = sys.argv[1]         # name of sequencing sample, e.g. PJ018
diffex_INFILE = sys.argv[2]    # name of SCDE output file containing differential expression analysis
matrix_INFILE = sys.argv[3]    # file containing matrix of molecular counts
norm_INFILE = sys.argv[4]      # three-column file - col 1 is cell-identifying barcode, col 2 is number of molecules detected, col 3 is number of genes detected 
tsne_INFILE = sys.argv[5]      # file containing two-dimensional coordinates of t-SNE visualization
fc_THRESH = float(sys.argv[6]) # required absolute log2(fold-change) for display
pv_THRESH = float(sys.argv[7]) # required adjusted p-value for display
tsne_OUTFILE = diffex_INFILE+'.tsne.pdf'

# Get differentially expressed genes (based on fold-change and p-value thresholds).
genes = [] 
with open(diffex_INFILE) as f:
	for line in f:
		llist = line.split()
		if float(llist[5]) < pv_THRESH:
			if abs(float(llist[3])) > fc_THRESH:
				genes.append(llist[0])
Ngenes = len(genes)
print(Ngenes)
if Ngenes > 100: # If there are more than 100 differentially expressed genes, take the top 100 based on adjusted p-value
	genes = genes[0:100]

# Run t-SNE plotting program to color-code the t-SNE visualization by the differentially expessed genes.
cmd = 'python3.5 gene_plot.py %(run_NAME)s %(matrix_INFILE)s %(norm_INFILE)s 2 3 %(tsne_INFILE)s %(tsne_OUTFILE)s' % vars()
cmd=cmd+' '+' '.join(genes)
print(cmd)
os.system(cmd)
