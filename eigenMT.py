from __future__ import print_function
import os
import sys
import fileinput
import argparse
import numpy as np
import pandas as pd
import scipy.linalg as splin
import gc
import gzip
from sklearn import covariance

##############FUNCTIONS

def open_file(filename):
    """Function to open a (potentially gzipped) file."""
    with open(filename, 'rb') as file_connection:
        file_header = file_connection.readline()
    if file_header.startswith(b"\x1f\x8b\x08"):
        opener = gzip.open(filename, 'rt')
    else:
        opener = open(filename)
    return opener

def load_tensorqtl_output(tensorqtl_parquet, group_size_s=None):
    """Read tensorQTL output"""
    df = pd.read_parquet(tensorqtl_parquet)
    if 'gene_id' not in df:
        if ':' in df['phenotype_id'].iloc[0]:
            df['gene_id'] = df['phenotype_id'].apply(lambda x: x.rsplit(':',1)[1] if ':' in x else x)
        else:
            df.rename(columns={'phenotype_id':'gene_id'}, inplace=True)
    # eigenMT requires a 'p-value' column (see make_test_dict); first column must be variant, second gene/phenotype
    df = df[['variant_id', 'gene_id']+[i for i in df.columns if i not in ['variant_id', 'gene_id']]]
    # select p-value column
    if 'pval_nominal' in df.columns:
        df['p-value'] = df['pval_nominal'].copy()
    elif 'pval_gi' in df.columns:  # interaction model
        df['p-value'] = df['pval_gi'].copy()
    if group_size_s is not None:
        print('  * adjusting p-values by phenotype group size')
        df['p-value'] = np.minimum(df['p-value']*df['gene_id'].map(group_size_s), 1.0)
    return df

def make_genpos_dict(POS_fh, CHROM):
    '''
    Function to read in SNPs and their positions and make a dict.
    Keys are SNP IDs; values are positions.
    Only stores the information of the SNPs from the given chromosome.
    '''
    pos_dict = {}
    with open_file(POS_fh) as POS:
        POS.readline()  # skip header
        for line in POS:
            line = line.rstrip().split()
            if line[1] == CHROM:
                pos_dict[line[0]] = float(line[2])
    return pos_dict

def make_phepos_dict(POS_fh, CHROM):
    '''
    Function to read in phenotypes (probes, genes, peaks) with their start and end positions and make a dict.
    Keys are phenotype IDs; values are start and end positions.
    Only stores the information of the phenotypes from the given chromosome.
    '''
    pos_dict = {}
    with open_file(POS_fh) as POS:
        POS.readline()  # skip header
        for line in POS:
            line = line.rstrip().split()
            if line[1] == CHROM:
                pos_array = np.array(line[2:4])
                pos_dict[line[0]] = np.float64(pos_array)
    return pos_dict

def make_gen_dict(GEN_fh, pos_dict, sample_ids=None):
    '''
    Function to read in genotype matrix from MatrixEQTL and make dict.
    Keys are SNP positions; values are genotypes.
    '''
    gen_dict = {}
    with open_file(GEN_fh) as GEN:
        header = GEN.readline().rstrip().split()
        if sample_ids is not None:
            ix = [header[1:].index(i) for i in sample_ids]
        for line in GEN: #Go through each line of the genotype matrix and add line to gen_dict
            line = line.rstrip().split()
            snp = pos_dict[line[0]]
            genos = np.array(line[1:])
            if sample_ids is not None:
                genos = genos[ix]
            genos[genos == 'NA'] = -1  # no effect if already -1
            genos = np.float64(genos)
            genos[genos == -1] = np.mean(genos[genos != -1])
            gen_dict[snp] = genos
    return gen_dict  # pos->genotypes

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

def make_test_dict_tensorqtl(QTL_fh, gen_dict, genpos_dict, cis_dist, group_size_s=None):
    """
    Same arguments and output as make_test_dict, for output from tensorQTL.
    QTL_fh: parquet file with variant-gene pair associations
    """
    qtl_df = load_tensorqtl_output(QTL_fh, group_size_s=group_size_s)
    qtl_df = qtl_df[qtl_df['tss_distance'].abs()<=cis_dist]

    gdf = qtl_df.groupby('gene_id')
    test_dict = {}
    for gene_id,g in gdf:
        g0 = g.loc[g['p-value'].idxmin()]
        test_dict[gene_id] = {
            'snps': [genpos_dict[i] for i in g['variant_id']],  # variant positions
            'best_snp':genpos_dict[g0['variant_id']],
            'pval':g0['p-value'],
            'line':'\t'.join([i if isinstance(i, str) else '{:.6g}'.format(i)  for i in g0.values])
        }
    return test_dict, '\t'.join(qtl_df.columns)

def make_test_dict_external(QTL_fh, gen_dict, genpos_dict, phepos_dict, cis_dist):
    '''
    Same arguments and output as make_test_dict.
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
    TSSs, genes = [list(x) for x in zip(*sorted(zip(TSSs, genes), key=lambda p: p[0]))]
    for gene in genes:
        perc = (100 * counter / numgenes)
        if (counter % 100) == 0:
            print(str(counter) + ' out of ' + str(numgenes) + ' completed ' + '(' + str(round(perc, 3)) + '%)', flush=True)
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
    eigenvalues = np.sort(eigenvalues)[::-1]
    running_sum = 0
    counter = 0
    while running_sum < variance * var_thresh:
        running_sum += eigenvalues[counter]
        counter += 1
    return counter


##############MAIN
if __name__=='__main__':
    USAGE = """
    Takes in SNP-gene tests from MatrixEQTL output and performs gene level Bonferroni correction using eigenvalue decomposition of
    the genotype correlation matrix. Picks best SNP per gene.
    """

    parser = argparse.ArgumentParser(description = USAGE)
    parser.add_argument('--QTL', required = True, help = 'Matrix-EQTL output file for one chromosome')
    parser.add_argument('--GEN', required = True, help = 'genotype matrix file')
    parser.add_argument('--var_thresh', type=float, default = 0.99, help = 'variance threshold')
    parser.add_argument('--OUT', required = True, help = 'output filename')
    parser.add_argument('--window', type=int, default = 200, help = 'SNP window size')
    parser.add_argument('--GENPOS', required = True, help = 'map of genotype to chr and position (as required by Matrix-eQTL)')
    parser.add_argument('--PHEPOS', required = True, help = 'map of measured phenotypes to chr and position (eg. gene expression to CHROM and TSS; as required by Matrix-eQTL)')
    parser.add_argument('--CHROM', required = True, help = 'Chromosome that is being processed (must match format of chr in POS)')
    parser.add_argument('--cis_dist', type=float, default = 1e6, help = 'threshold for bp distance from the gene TSS to perform multiple testing correction (default = 1e6)')
    parser.add_argument('--external', action = 'store_true', help = 'indicates whether the provided genotype matrix is different from the one used to call cis-eQTLs initially (default = False)')
    parser.add_argument('--sample_list', default=None, help='File with sample IDs (one per line) to select from genotypes')
    parser.add_argument('--phenotype_groups', default=None, help='File with phenotype_id->group_id mapping')
    args = parser.parse_args()

    ##Make SNP position dict
    print('Processing genotype position file.', flush=True)
    genpos_dict = make_genpos_dict(args.GENPOS, args.CHROM)

    ##Make phenotype position dict
    print('Processing phenotype position file.', flush=True)
    phepos_dict = make_phepos_dict(args.PHEPOS, args.CHROM)

    ##Make genotype dict
    print('Processing genotype matrix.', flush=True)
    if args.sample_list is not None:
        with open(args.sample_list) as f:
            sample_ids = f.read().strip().split('\n')
        print('  * using subset of '+str(len(sample_ids))+' samples.')
    else:
        sample_ids = None
    gen_dict = make_gen_dict(args.GEN, genpos_dict, sample_ids)

    ##Make SNP-gene test dict
    if not args.external:
        if args.QTL.endswith('.parquet'):
            print('Processing tensorQTL tests file.', flush=True)
            if args.phenotype_groups is not None:
                group_s = pd.read_csv(args.phenotype_groups, sep='\t', index_col=0, header=None, squeeze=True)
                group_size_s = group_s.value_counts()
            else:
                group_size_s = None
            test_dict, input_header = make_test_dict_tensorqtl(args.QTL, gen_dict, genpos_dict, args.cis_dist, group_size_s=group_size_s)
        else:
            print('Processing Matrix-eQTL tests file.', flush=True)
            test_dict, input_header = make_test_dict(args.QTL, gen_dict, genpos_dict, phepos_dict, args.cis_dist)
    else:
        print('Processing Matrix-eQTL tests file. External genotype matrix and position file assumed.', flush=True)
        test_dict, input_header = make_test_dict_external(args.QTL, gen_dict, genpos_dict, phepos_dict, args.cis_dist)

    ##Perform BF correction using eigenvalue decomposition of the correlation matrix
    print('Performing eigenMT correction.', flush=True)
    bf_eigen_windows(test_dict, gen_dict, phepos_dict, args.OUT, input_header, args.var_thresh, args.window)
