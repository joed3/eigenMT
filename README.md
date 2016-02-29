# eigenMT
An Efficient Multiple-Testing Adjustment for eQTL Studies that Accounts for Linkage Disequilibrium between Variants

Introduction
------------


Download
------------
Begin by downloading or cloning this respository. eigenMT runs as a stand-alone python script. It requires an active installaton of python (version 2.7 or higher) to be installed. To download and install a python distribution, there are a few convenient options:
- Anaconda: https://store.continuum.io/cshop/anaconda/
- Enthought Canopy: https://www.enthought.com/products/canopy/

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

--CHROM
Chromosome ID. Indicates which chromosome the analysis will be performed on.
Must match the ID used in the Matrix-eQTL SNP-gene tests file.

--QTL
Filename for Matrix-eQTL SNP-gene tests file.
See qtls.txt for an example.

--GEN
Genotype matrix in Matrix-eQTL format. Can accept either hardcoded or dosage based genotypes. 
Missing genotypes must be encoded as NA and will be imputed to the mean genotype during correction.
See genotypes.txt for an example.

--GENPOS
Matrix-eQTL genotype position file.
See positions.txt for an example.

--var_thresh
Threshold for amount of variance explained in the genotype correlation matrix. 
Default is 99% variance explained. Increasing this threshold will increase estimates of effective number of tests (M_eff) and decrease accuracy of the approximation. 

--OUT
Output filename. Output format is described below.

--window
Window size parameter. Determines what size of disjoint windows to split genotype matrices for each gene into.
Default is 200 SNPs. We recommend using a window size of at least 50 SNPs up to 200 SNPs to balance accuracy and speed.

Output
------------
The output file is in tab-separated format with the following columns:
- Col 1: SNP ID
- Col 2: GENE ID
- Col 3: estimate of effect size BETA from Matrix-eQTL
- Col 4: T-statistic from Matrix-eQTL
- Col 5: p-value from Matrix-eQTL
- Col 6: eigenMT corrected p-value
- Col 7: estimated number of independent tests for the gene

Note: each tested gene will appear once in the output file with it's most significant SNP and the eigenMT corrected p-value.


Example
------------
We offer a small example dataset for testing. We provide a genotype matrix and corresponding position file in Matrix-eQTL format. This genotype matrix is for chromosome 19 for the EUR373 samples as part of the GEUVADIS study. We also provide a sample of the cis-eQTL tests performed for these samples using Matrix-eQTL. This sample SNP-gene tests file represents 100 genes. To run eigenMT on the example data, use the following command:
```
python eigenMT.py --CHROM 19 \
	--QTL qtls.txt \
	--GEN genotypes.txt \
	--POS positions.txt \
	--OUT exampleOut.txt
```
Note: this example uses the default settings for window size and variance explained threshold.


Population stratification and other covariates
------------
We recommend first removing the effects of population stratification (genotype PCs, population or ancestry assignments) and other covariates (age, gender, PEER factors, etc) from the expression matrix. We have shown that by first removing these effects and then performing cis-eQTL calling on the inverse rank normalized residuals, we maintain conservativeness and accuracy. We will soon provide an R script with example code for how to perform this normalization step prior to a run of eigenMT. 


Citation
------------
Davis JR, Fresard L, Knowles DA, Pala M, Bustamante CD, Battle A, Montgomery SB (2015). An Efficient Multiple-Testing Adjustment for eQTL Studies that Accounts for Linkage Disequilibrium between Variants. The American Journal of Human Genetics, 98(1), 216–224. doi: (http://doi.org/10.1016/j.ajhg.2015.11.021)

Contact
------------
-Joe Davis: joed3@stanford.edu
-Laure Fresard: lfresard@stanford.edu
-Stephen Montgomery: smontgom@stanford.edu
