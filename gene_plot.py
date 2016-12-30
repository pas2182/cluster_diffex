#! /usr/bin/python3.4
import numpy as np
from math import exp,log10,sqrt,ceil
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os


# This program is part of the cluster_diffex pipeline. The purpose of this program is to plot t-SNE visualizations (generally those)
# generated using the similarity_pipeline) colored based on the differentially expressed genes found between clusters by the 
# cluster_diffex pipleine.

run = sys.argv[1]                 # name of sequencing sample, e.g. PJ018
matrix_INFILE = sys.argv[2]       # name of counts matrix file - typically, col 1 is gene id, col 2 is gene symbol, col 3 contains molecular counts of first cell
norm_INFILE = sys.argv[3]         # three-column file - col 1 is cell-identifying barcode, col 2 is number of molecules detected, col 3 is number of genes detected 
gene_COLUMN = int(sys.argv[4])    # column number of matrix_INFILE to use for gene names
ct_COLUMN = int(sys.argv[5])      # column number of matrix_INFILE where the counts for first cell occur
tsne_OUTFILE = sys.argv[6]        # file containing two-dimensional coordinates of t-SNE visualization
tsne_PDF = sys.argv[7]            # PDF output for colored t-SNE plots
first_markers = []
for i in range(8,len(sys.argv)):  # names of differentially expressed genes
        first_markers.append(sys.argv[i])

markers = set(first_markers)

# Load the matrix of molecular counts into memory along with gene names in a separate list.
print('Loading count matrix...')
genes = []
matrix = []
with open(matrix_INFILE) as f:
	for line in f:
		llist = line.split()
		gene = llist[gene_COLUMN-1]
		if gene in markers:
			matrix.append([int(llist[i]) for i in range(ct_COLUMN-1,len(llist))])
			genes.append(gene)

Ngenes = len(genes)
Ncells = len(matrix[0])
fNcells = float(Ncells)
print('Found %(Ngenes)d genes in %(Ncells)d cells...' % vars())

# Load the total molecular counts for each cell into memory.
print('Normalizing matrix...')
norm = []
with open(norm_INFILE) as f:
	for line in f:
		llist = line.split()
		norm.append(float(llist[1]))	

# Normalize the matrix of molecular counts by total counts for each cell.
nmatrix = [[float(cts[i])/norm[i] for i in range(0,Ncells)] for cts in matrix]

# Load the t-SNE coordinates into memory.
tsne_out = []
with open(tsne_OUTFILE) as h:
	for line in h:
		llist = line.split()
		tsne_out.append([float(llist[0]),float(llist[1])])
tsne_out = np.array(tsne_out)
tsne_x = tsne_out[:,0]
tsne_y = tsne_out[:,1]
xmx = max(abs(np.array(tsne_x)))
ymx = max(abs(np.array(tsne_y)))
axis_mx = max((xmx,ymx))*1.1 # adjust for adequate, square axes

# Plot t-SNE visualizations color-coded by expression of each differentially expressed gene
with PdfPages(tsne_PDF) as pdf:
	plt.plot(tsne_out[:,0],tsne_out[:,1],'ko',markersize=4) # first plot a black-and-white t-SNE
	plt.xlabel('t-SNE Axis 1')
	plt.ylabel('t-SNE Axis 2')
	plt.xlim([-axis_mx,axis_mx])
	plt.ylim([-axis_mx,axis_mx])
	ax=plt.gca()
	ax.set_aspect('equal')
	pdf.savefig()
	plt.close()
	i=0
	for marker in first_markers: # for each differentially expressed gene
		i = genes.index(marker)
		expression = np.log2(np.array(nmatrix[i])+1.0e-5) # Use log-scale, normalized expression for coloring
		gene = marker
		plt.scatter(tsne_out[:,0],tsne_out[:,1],c=expression,s=10,lw=0,edgecolor='none',cmap='cool')
		plt.xlabel('t-SNE Axis 1')
		plt.ylabel('t-SNE Axis 2')
		plt.xlim([-axis_mx,axis_mx])
		plt.ylim([-axis_mx,axis_mx])
		ax=plt.gca()
		ax.set_aspect('equal')
		plt.title(gene)
		pdf.savefig()
		plt.close()
		i+=1













	
