import os
import sys
import fileinput
import argparse
import numpy as np
import scipy as sp
import scipy.linalg as splin 
import gc
import gzip
from sklearn import covariance

##############FUNCTIONS

def open_file(filename):
	'''
	Function to open a (potentially gzipped) file.
	'''
	file_connection = open(filename)
	file_header = file_connection.readline()
	file_connection.close()
	if file_header.startswith("\x1f\x8b\x08"):
		file_connection = gzip.open(filename, 'rb')
	else:
		file_connection = open(filename)
	return file_connection

def make_genpos_dict(POS_fh, CHROM):
	'''
	Function to read in SNPs and their positions and make a dict.
	Keys are SNP IDs; values are positions.
	Only stores the information of the SNPs from the given chromosome.
	'''
 
	POS = open_file(POS_fh)
	POS.readline()

	pos_dict = {}
	
	for line in POS:
		line = line.rstrip().split()
		if line[1] == CHROM:
			pos_dict[line[0]] = float(line[2])
	POS.close()
	return pos_dict

def make_phepos_dict(POS_fh, CHROM):
	'''
	Function to read in phenotypes (probes, genes, peaks) with their start and end positions and make a dict.
	Keys are phenotype IDs; values are start and end positions.
	Only stores the information of the phenotypes from the given chromosome.
	'''
	POS = open_file(POS_fh)
	POS.readline()

	pos_dict = {}
	
	for line in POS:
		line = line.rstrip().split()
		if line[1] == CHROM:
			pos_array = np.array(line[2:4])
			pos_dict[line[0]] = np.array(map(float, pos_array))
	POS.close()
	return pos_dict

def make_gen_dict(GEN_fh, pos_dict):
	'''
	Function to read in genotype matrix from MatrixEQTL and make dict.
	Keys are SNP positions; values are genotypes.
	'''
	GEN = open_file(GEN_fh)
	GEN.readline()

	gen_dict = {}

	for line in GEN: #Go through each line of the genotype matrix and add line to gen_dict
		line = line.rstrip().split()
		snp = pos_dict[line[0]]
		genos = np.array(line[1:])
		genos[genos == 'NA'] = -1
		genos = np.array(map(float, genos))
		genos[genos == -1] = np.mean(genos[genos != -1])
		gen_dict[snp] = genos
	GEN.close()
	return gen_dict

def make_test_dict(QTL_fh, gen_dict, genpos_dict, phepos_dict, cis_dist):
	'''
	Function to make dict of SNP-gene tests. Also returns the header of the file.
	Keys are gene IDs. Values are dict giving list of tested SNPs and the best SNP and it's p-value.
	SNPs are identified by their positions.
	Assumes the first column of the input file is the SNP ID and the second column is the GENE.
	The column with the p-value is determined by the function, which allows for a more flexible input format.
	'''
	QTL = open_file(QTL_fh)
	header = QTL.readline().rstrip().split()
	# find the column with the p-value
	if 'p-value' in header:
		pvalIndex = header.index('p-value') 
	elif 'p.value' in header:
		pvalIndex = header.index('p.value')
	elif 'pvalue' in header:
		pvalIndex = header.index('pvalue')
	else:
		sys.exit('Cannot find the p-value column in the tests file.')

	test_dict = {}

	for line in QTL:
		line = line.rstrip().split()
		if line[0] in genpos_dict:
			snp = genpos_dict[line[0]]
			gene = line[1]
			if snp in gen_dict and gene in phepos_dict:
				phepos = phepos_dict[gene]
				distance = min(abs(phepos - snp))
				if distance <= cis_dist:
					pval = line[pvalIndex]
					pval = float(pval)
					if gene not in test_dict:
						test_dict[gene] = {'snps' : [snp], 'best_snp' : snp, 'pval' : pval, 'line' : '\t'.join(line)}
					else:
						if pval < test_dict[gene]['pval']:
							test_dict[gene]['best_snp'] = snp
							test_dict[gene]['pval'] = pval
							test_dict[gene]['line'] = '\t'.join(line)
						test_dict[gene]['snps'].append(snp)

	QTL.close()
	return test_dict, "\t".join(header)

def make_test_dict_external(QTL_fh, gen_dict, genpos_dict, phepos_dict, cis_dist):
	'''
	Same arguments and output as the function above.
	Main difference from the previous function is that it assumes the genotype matrix and position file
	are separate from that used in the Matrix-eQTL run. This function is to be used with the external option 
	to allow calculation of the effective number of tests using a different, preferably larger, genotype sample.
	'''
	QTL = open_file(QTL_fh)
	header = QTL.readline().rstrip().split()
	# find the column with the p-value
	if 'p-value' in header:
		pvalIndex = header.index('p-value') 
	elif 'p.value' in header:
		pvalIndex = header.index('p.value')
	elif 'pvalue' in header:
		pvalIndex = header.index('pvalue')
	else:
		sys.exit('Cannot find the p-value column in the tests file.')

	test_dict = {}

	for line in QTL:
		line = line.rstrip().split()
		if line[0] in genpos_dict:
			snp = genpos_dict[line[0]]
			gene = line[1]
			if snp in gen_dict and gene in phepos_dict:
				phepos = phepos_dict[gene]
				distance = min(abs(phepos - snp))
				if distance <= cis_dist:
					pval = line[pvalIndex]
					pval = float(pval)
					if gene not in test_dict:
						test_dict[gene] = {'best_snp' : snp, 'pval' : pval, 'line' : '\t'.join(line)}
					else:
						if pval < test_dict[gene]['pval']:
							test_dict[gene]['best_snp'] = snp
							test_dict[gene]['pval'] = pval
							test_dict[gene]['line'] = '\t'.join(line)
					
	QTL.close()
	snps = np.array(genpos_dict.values())
	for gene in test_dict:
		phepos = phepos_dict[gene]
		#Calculate distances to phenotype start and end positions
		is_in_cis_start = abs(snps - phepos[0]) <= cis_dist
		is_in_cis_end = abs(snps - phepos[1]) <= cis_dist
		test_dict[gene]['snps'] = snps[is_in_cis_start | is_in_cis_end] 
	return test_dict, "\t".join(header)


def bf_eigen_windows(test_dict, gen_dict, phepos_dict, OUT_fh, input_header, var_thresh, window):
	'''
	Function to process dict of SNP-gene tests.
	Calculates the genotype correlation matrix for the SNPs tested for each gene using windows around the gene.
	Will calculate a regularized correlation matrix from the Ledoit-Wolf estimator of the covariance matrix.
	Finds the eigenvalues for this matrix.
	Calculates how many eigenvalues are needed to reach the variance threshold.
	This final value will be the effective Bonferroni correction number.
	Will output to file corrected pvalue for best SNP per gene as it goes along.
	'''
	OUT = open(OUT_fh, 'w')
	OUT.write(input_header + '\tBF\tTESTS\n')

	counter = 1.0
	genes = test_dict.keys()
	numgenes = len(genes)
	TSSs = []
	for gene in genes:
		TSSs.append(phepos_dict[gene][0])
	zipped = zip(TSSs, genes)
	zipped.sort()
	TSSs, genes = zip(*zipped)
	for gene in genes:
		perc = (100 * counter / numgenes)
		if (counter % 100) == 0:
			print str(counter) + ' out of ' + str(numgenes) + ' completed ' + '(' + str(round(perc, 3)) + '%)' 
		sys.stdout.flush()
		counter += 1
		snps = np.sort(test_dict[gene]['snps'])
		start = 0
		stop = window
		M = len(snps)
		m_eff = 0
		window_counter = 0
		while start < M:
			if stop - start == 1:
				m_eff += 1
				break ##can't compute eigenvalues for a scalar, so add 1 to m_eff and break from the while loop
			snps_window = snps[start:stop]
			genotypes = []
			for snp in snps_window:
				if snp in gen_dict:
					genotypes.append(gen_dict[snp])
			genotypes = np.matrix(genotypes)
			m, n = np.shape(genotypes)
			gen_corr, alpha = lw_shrink(genotypes)
			window_counter += 1
			eigenvalues = splin.eigvalsh(gen_corr)
			eigenvalues[eigenvalues < 0] = 0
			m_eff += find_num_eigs(eigenvalues, m, var_thresh)
			start += window
			stop += window
			if stop > M:
				stop = M
		OUT.write(test_dict[gene]['line'] + '\t' + str(min(test_dict[gene]['pval'] * m_eff, 1)) + '\t' + str(m_eff) + '\n')
		OUT.flush()
		gc.collect()
	OUT.close()
	
def lw_shrink(genotypes):
	'''
	Function to obtain smoothed estimate of the genotype correlation matrix.
	Uses the method proposed by Ledoit and Wolf to estimate the shrinkage parameter alpha.
	Input: genotype matrix.
	Output: smoother correlation matrix and the estimated shrinkage parameter.
	'''
	lw = covariance.LedoitWolf()
	m, n = np.shape(genotypes)
	try:	
		fitted = lw.fit(genotypes.T)
		alpha = fitted.shrinkage_
		shrunk_cov = fitted.covariance_
		shrunk_precision = np.mat(np.diag(np.diag(shrunk_cov)**(-.5)))
		shrunk_cor = shrunk_precision * shrunk_cov * shrunk_precision
	except: #Exception for handling case where SNPs in the window are all in perfect LD
		row = np.repeat(1, m)
		shrunk_cor = []
		for i in range(0,m):
			shrunk_cor.append(row)
		shrunk_cor = np.mat(shrunk_cor)
		alpha = 'NA'
	return shrunk_cor, alpha

def find_num_eigs(eigenvalues, variance, var_thresh):
	'''
	Function to find the number of eigenvalues required to reach a certain threshold of variance explained.
	'''
	eigenvalues = np.sort(eigenvalues).tolist()
	eigenvalues.reverse()
	running_sum = 0
	counter = 0
	while running_sum < variance * var_thresh:
		running_sum += eigenvalues[counter]
		counter += 1
	return counter


##############MAIN
USAGE = """
	Takes in SNP-gene tests from MatrixEQTL output and performs gene level Bonferroni correction using eigenvalue decomposition of
	the genotype correlation matrix. Picks best SNP per gene.
    """

parser = argparse.ArgumentParser(description = USAGE)
parser.add_argument('--QTL', dest = 'QTL', required = True, help = 'Matrix-EQTL output file for one chromosome')
parser.add_argument('--GEN', dest = 'GEN', required = True, help = 'genotype matrix file')
parser.add_argument('--var_thresh', dest = 'var_thresh', default = 0.99, help = 'variance threshold')
parser.add_argument('--OUT', dest = 'OUT', required = True, help = 'output filename')
parser.add_argument('--window', dest = 'window', default = 200, help = 'SNP window size')
parser.add_argument('--GENPOS', dest = 'GENPOS', required = True, help = 'map of genotype to chr and position (as required by Matrix-eQTL)')
parser.add_argument('--PHEPOS', dest = 'PHEPOS', required = True, help = 'map of measured phenotypes to chr and position (eg. gene expression to CHROM and TSS; as required by Matrix-eQTL)')
parser.add_argument('--CHROM', dest = 'CHROM', required = True, help = 'Chromosome that is being processed (must match format of chr in POS)')
parser.add_argument('--cis_dist', dest = 'cis_dist', default = 1e6, help = 'threshold for bp distance from the gene TSS to perform multiple testing correction (default = 1e6)')
parser.add_argument('--external', dest = 'external', action = 'store_true', help = 'indicates whether the provided genotype matrix is different from the one used to call cis-eQTLs initially (default = False)')

args = parser.parse_args()

QTL_fh = args.QTL 
GEN_fh = args.GEN
OUT_fh = args.OUT
GENPOS_fh = args.GENPOS
PHEPOS_fh = args.PHEPOS
CHROM = args.CHROM
var_thresh = float(args.var_thresh)
window = int(args.window)
cis_dist = float(args.cis_dist)
external = args.external

##Make SNP position dict
print 'Processing genotype position file.'
sys.stdout.flush()
genpos_dict = make_genpos_dict(GENPOS_fh, CHROM)

##Make phenotype position dict
print 'Processing phenotype position file.'
sys.stdout.flush()
phepos_dict = make_phepos_dict(PHEPOS_fh, CHROM)

##Make genotype dict
print 'Processing genotype matrix.'
sys.stdout.flush()
gen_dict = make_gen_dict(GEN_fh, genpos_dict)

##Make SNP-gene test dict
if not external:
	print 'Processing Matrix-eQTL tests file.'
	sys.stdout.flush()
	test_dict, input_header = make_test_dict(QTL_fh, gen_dict, genpos_dict, phepos_dict, cis_dist)
else:
	print 'Processing Matrix-eQTL tests file. External genotype matrix and position file assumed.'
	sys.stdout.flush()
	test_dict, input_header = make_test_dict_external(QTL_fh, gen_dict, genpos_dict, phepos_dict, cis_dist)

##Perform BF correction using eigenvalue decomposition of the correlation matrix
print 'Performing eigenMT correction.'
sys.stdout.flush()
bf_eigen_windows(test_dict, gen_dict, phepos_dict, OUT_fh, input_header, var_thresh, window)



