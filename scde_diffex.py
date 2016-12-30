#! /usr/bin/python3.5
import sys
import os

# This program is part of the cluster_diffex pipeline. It is designed to run the SCDE algorithm from Kharchenko between a pair of
# single cell clusters (each of which is comprised of multiple single cell profiles). The SCDE algorithm is written in R, and this
# Python program writes a custom R script for each job and runs it.


matrix1_INFILE = sys.argv[1]   # matrix of molecular counts for first cluster of single cells
matrix2_INFILE = sys.argv[2]   # matrix of molecular counts for second cluster of single cells
diffex_OUTFILE = sys.argv[3]   # SCDE output file
run_OUTFILE = sys.argv[4]      # custom R script for running SCDE
output_DIRECTORY = sys.argv[5] # name of output directory
group1_NAME = sys.argv[6]      # name of first cluster
group2_NAME = sys.argv[7]      # name of second cluster
genenames_INFILE = sys.argv[8] # file containing gene ID-gene symbol conversion table
proc_N = int(sys.argv[9])      # number of threads

# Prepare and format input matrix for SCDE.
cmd = 'mkdir %(output_DIRECTORY)s' % vars() # create output directory
print(cmd)
os.system(cmd)
cmd = 'cut -f2- %(matrix2_INFILE)s > %(matrix2_INFILE)s_tmp' % vars() # isolate counts from matrix2_INFILE in temp file (remove columns with gene ID, gene symbol info)
print(cmd)
os.system(cmd)
cmd = 'paste %(matrix1_INFILE)s %(matrix2_INFILE)s_tmp > %(output_DIRECTORY)s/%(group1_NAME)s_%(group2_NAME)s.cts.txt' % vars() 
print(cmd)
os.system(cmd) # concatenate counts matrices for two clusters
cmd = 'rm %(matrix2_INFILE)s_tmp' % vars() # remove temporary file
print(cmd)
os.system(cmd)

# Generate custom R script for running SCDE. 
run_INFILE = output_DIRECTORY+'/'+run_OUTFILE
run_INPUT = open(run_INFILE,'w')
run_INPUT.write('library(scde)\n')
run_INPUT.write('counts<-read.delim(\'%(output_DIRECTORY)s/%(group1_NAME)s_%(group2_NAME)s.cts.txt\',skip=0,header=TRUE,sep=\"\\t\",stringsAsFactors=FALSE,row.names=1)\n' % vars())
run_INPUT.write('cell.labels<-substr(colnames(counts),0,6)\n')
run_INPUT.write('groups<-factor(cell.labels,levels=c(\"%(group1_NAME)s\",\"%(group2_NAME)s\"))\n' % vars())
run_INPUT.write('table(groups)\n')
run_INPUT.write('counts<-counts[rowSums(counts)>0,]\n') 
run_INPUT.write('ord<-order(groups)\n')
run_INPUT.write('n.cores<-%(proc_N)d\n' % vars())
run_INPUT.write('scde.fitted.model <- scde.error.models(counts=counts,groups=groups,n.cores=n.cores,save.model.plots=FALSE)\n')
run_INPUT.write('scde.prior <- scde.expression.prior(models=scde.fitted.model,counts=counts)\n')
run_INPUT.write('ediff <- scde.expression.difference(scde.fitted.model,counts,scde.prior,groups=groups,n.cores=n.cores)\n')
run_INPUT.write('p.values<-2*pnorm(abs(ediff$Z),lower.tail=FALSE)\n')
run_INPUT.write('p.values.adj<-2*pnorm(abs(ediff$cZ),lower.tail=FALSE)\n')
run_INPUT.write('significant.genes<-which(p.values.adj<0.05)\n')
run_INPUT.write('length(significant.genes)\n')
run_INPUT.write('ord<-order(p.values.adj[significant.genes])\n')
run_INPUT.write('de<-cbind(ediff[significant.genes,1:3],p.values.adj[significant.genes])[ord,]\n')
run_INPUT.write('colnames(de)<-c(\"Lower bound\",\"log2 fold change\",\"Upper bound\",\"p-value\")\n')
run_INPUT.write('write.table(de,file=\"%(output_DIRECTORY)s/%(diffex_OUTFILE)s\")\n' % vars())
run_INPUT.close()

# Run custom R script for SCDE.
cmd = 'Rscript %(run_INFILE)s' % vars()
print(cmd)
os.system(cmd)

# Generate dictionary from genenames_INFILE table with gene IDs as keys
# and gene symbols as values.
gene_dict = {}
with open(genenames_INFILE) as f:
	for line in f:
		llist = line.split()
		gene_dict[llist[0]] = llist[1]

# Open output of SCDE and add gene symbols in a second output file.
diffex_OUTFILE = output_DIRECTORY+'/'+diffex_OUTFILE
diffex2_OUTFILE = diffex_OUTFILE+'.conv.txt'
diffex2_OUTPUT = open(diffex2_OUTFILE,'w')
i=0
with open(diffex_OUTFILE) as f:
	for line in f:
		if i == 1:
			gid = line.split()[0][1:-1]
			gene = gene_dict[gid]
			newline = gene+'\t'+line
			diffex2_OUTPUT.write(newline)
		i = 1
diffex2_OUTPUT.close()
