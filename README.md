# eigenMT
An Efficient Multiple-Testing Adjustment for eQTL Studies that Accounts for Linkage Disequilibrium between Variants

Introduction
------------
eigenMT is a computationally efficient multiple testing correction method for cis-eQTL studies. Typically, modern cis-eQTL studies correct for multiple testing across thousands of variants for a single gene using permutations. This correction method, though accurate, requires a high computational cost on the order of days to weeks for large numbers of permutations and/or large sample sizes. Our method reduces this computational burden by orders of magnitude while maintaining high accuracy when compared to permutations. To accomplish this, our method attempts to estimate the number of effective tests performed for a given gene by directly accounting for the LD relationship among the tested variants.  

Download
------------
Begin by downloading or cloning this respository. eigenMT runs as a stand-alone python script. It requires an active installaton of python (version 2.7 or higher) to be installed. To download and install a python distribution, there are a few convenient options:
- [Anaconda](https://store.continuum.io/cshop/anaconda/)
- [Enthought Canopy](https://www.enthought.com/products/canopy/)

These bundled installations already include a number of modules required for running eigenMT including:
- [numpy](http://www.numpy.org/)
- [scipy](http://www.scipy.org/)
- [Scikit-learn](http://scikit-learn.org/stable/)

Again, these modules should come pre-packaged in one of the bundled python installations above. If you decide to install these packages yourself, please see the listed websites for detailed instructions.

Input
------------
A typical run will have the following command line format:
```
python eigenMT.py --CHROM <chromosome ID>
	--QTL <Matrix-eQTL SNP-gene tests file> \
	--GEN <Matrix-eQTL genotype matrix> \
	--POS <Matrix-eQTL genotype position file> \
	--var_thresh [variance explained threshold, default 0.99] \
	--OUT <output filename> \
	--window [window size, default 200]
```

The input fields are defined as follows.

Argument        | Description
---------------------------  |-------------
CHROM                        | Chromosome ID. Indicates which chromosome the analysis will be performed on. Must match the ID used in the Matrix-eQTL SNP-gene tests file.
QTL                          | Filename for Matrix-eQTL SNP-gene tests file. See `example/qtls.txt` for an example.
GEN                          | Genotype matrix in Matrix-eQTL format. Can accept either hardcoded or dosage based genotypes. Missing genotypes must be encoded as NA and will be imputed to the mean genotype during correction. See `example/genotypes.txt` for an example.
GENPOS                       | Matrix-eQTL genotype position file. See `example/positions.txt` for an example.
var_thresh                   | Threshold for amount of variance explained in the genotype correlation matrix. Default is 99% variance explained. Increasing this threshold will increase estimates of effective number of tests (M_eff) and decrease accuracy of the approximation.
OUT                          | Output filename. Output format is described below.
window                       | Window size parameter. Determines what size of disjoint windows to split genotype matrices for each gene into. Default is 200 SNPs. We recommend using a window size of at least 50 SNPs up to 200 SNPs to balance accuracy and speed.
cis_dist                     | Threshold for bp distance from the gene TSS to perform multiple testing correction. Default is 1e6.
external                     | Logical indicating whether the provided genotype matrix is different from the one originally used to test for cis-eQTLs. 

Output
------------
The output file is in tab-separated format with the following columns.

Column          |  Description
--------------- |   ------------
1               |  variant ID
2               |  Gene ID
3               |  T-statistic from Matrix-eQTL
4               |  p-value from Matrix-eQTL 
5               |  FDR estimate from Matrix-eQTL
6               |  estimate of effect size BETA from Matrix-eQTL
7               |  eigenMT corrected p-value
8               |  estimated number of independent tests for the gene

Note: The first 5 columns of the output will correspond to the typical Matrix-eQTL output. Also, each tested gene will appear once in the output file with its most significant SNP and the eigenMT corrected p-value.


Example
------------
We offer a small example dataset for testing. We provide a genotype matrix and corresponding position file in Matrix-eQTL format. This genotype matrix is for chromosome 19 for the EUR373 samples as part of the [GEUVADIS](http://www.nature.com/nature/journal/v501/n7468/full/nature12531.html?WT.ec_id=NATURE-20130926) study. We also provide a sample of the cis-eQTL tests performed for these samples using Matrix-eQTL. This sample SNP-gene tests file represents 100 genes. To extract the example dataset, run
```
tar -zxvf example.tar.gz
```

To run eigenMT on the example data, use the following command:
```
python eigenMT.py --CHROM 19 \
	--QTL example/qtls.txt \
	--GEN example/genotypes.txt \
	--GENPOS example/gen.positions.txt \
	--PHEPOS example/phe.positions.txt \
	--OUT example/exampleOut.txt
```
Note: this example uses the default settings for window size and variance explained threshold.

External genotype data
------------
Our method offers the ability to perform multiple testing correction using a genotype matrix from a separate sample population than the one initially used for cis-eQTL testing. It is important to note that this external genotype matrix should come from the same background population as the one under study. For example, if cis-eQTL testing is performed in samples of European ancestry, then any external genotype matrix used should come from a similar European population. We have shown that using genotype data from studies with larger sample sizes can improve the accuracy of our method compared to using the genotype data for the study. Genotype data from larger studies will provide better estimates of the LD structure for variants around each gene, improving our estimates of the effective number of tests.

Population stratification and other covariates
------------
We recommend first removing the effects of population stratification (genotype PCs, population or ancestry assignments) and other covariates (age, gender, PEER factors, etc) from the expression matrix. We have shown that by first removing these effects and then performing cis-eQTL calling on the inverse rank normalized residuals, we maintain conservativeness and accuracy. We will soon provide an R script with example code for how to perform this normalization step prior to a run of eigenMT. 


Citation
------------
Davis JR, Fresard L, Knowles DA, Pala M, Bustamante CD, Battle A, Montgomery SB (2015). An Efficient Multiple-Testing Adjustment for eQTL Studies that Accounts for Linkage Disequilibrium between Variants. The American Journal of Human Genetics, 98(1), 216–224. doi: (http://doi.org/10.1016/j.ajhg.2015.11.021)

Contact
------------
- Joe Davis: joed3@stanford.edu
- Laure Fresard: lfresard@stanford.edu
- Stephen Montgomery: smontgom@stanford.edu
