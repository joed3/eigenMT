#### eigenMT
#### Authors: Joe Davis and Laure Fresard
#### Date: 7/16/15

#### 1. INSTALLATION

## eigenMT can be downloaded from: http://montgomerylab.stanford.edu/resources/eigenMT/eigenMT.html
##There are two files available for download:
	## 1. eigenMT.tgz
	## 2. eigenMTwithTestData.tgz
## To open eigenMT, run the following command from the directory holding the downloaded file:

tar -xzvf eigenMT.tgz

## OR

tar -xzvf eigenMTwithTestData.tgz

## This will create a directory named eigenMT
## Inside this directory will be the eigenMT python script (eigenMT.py) and this README.
## In the version with test data, there will also be example input files for a test run of the algorithm.

## eigenMT runs as a stand-alone python script. 
## It requires an active installaton of python (version 2.7 or higher) to be installed.
## To download and install a python distribution, there are a few convenient options
	## Anaconda: https://store.continuum.io/cshop/anaconda/
	## Enthought Canopy: https://www.enthought.com/products/canopy/
## These bundled installations already include a number of modules for running eigenMT
## The following are modules required for running eigenMT
	## numpy: http://www.numpy.org/
	## scipy: http://www.scipy.org/
	## Scikit-learn: http://scikit-learn.org/stable/
## Again, these modules should come pre-packaged in one of the bundled python installations above.
	## If you decide to install these packages yourself, please see the listed websites for detailed instructions.

#### 2. INPUT

## An eigenMT run will have the following command format:
## python eigenMT.py --CHROM <chromosome ID>
	## --QTL <Matrix-eQTL SNP-gene tests file> \
	## --GEN <Matrix-eQTL genotype matrix> \
	## --POS <Matrix-eQTL genotype position file> \
	## --var_thresh [variance explained threshold, default 0.99] \
	## --OUT <output filename> \
	## --window [window size, default 200]

## The input fields are defined as follows.

## --CHROM
## Chromosome ID. Indicates which chromosome the analysis will be performed on.
## Must match the ID used in the Matrix-eQTL SNP-gene tests file.

## --QTL
## Filename for Matrix-eQTL SNP-gene tests file.
## See qtls.txt for an example.

## --GEN
## Genotype matrix in Matrix-eQTL format. Can accept either hardcoded or dosage based genotypes. 
## Missing genotypes must be encoded as NA and will be imputed to the mean genotype during correction.
## See genotypes.txt for an example.

## --GENPOS
## Matrix-eQTL genotype position file.
## See positions.txt for an example.

## --var_thresh
## Threshold for amount of variance explained in the genotype correlation matrix. 
## Default is 99% variance explained. Increasing this threshold will increase estimates of effective number of tests (M_eff)
	## and decrease accuracy of the approximation. 

## --OUT
## Output filename. Output format is described below.

## --window
## Window size parameter. Determines what size of disjoint windows to split genotype matrices for each gene into.
## Default is 200 SNPs. We recommend using a window size of at least 50 SNPs up to 200 SNPs to balance accuracy and speed.

#### 3. OUTPUT

## The output file is in tab-separated format with the following columns:
	## Col 1: SNP ID
	## Col 2: GENE ID
	## Col 3: estimate of effect size BETA from Matrix-eQTL
	## Col 4: T-statistic from Matrix-eQTL
	## Col 5: p-value from Matrix-eQTL
	## Col 6: eigenMT corrected p-value
	## Col 7: estimated number of independent tests for the gene
## Note: each tested gene will appear once in the output file with it's most significant SNP and the eigenMT corrected p-value.

#### 4. EXAMPLE

## We offer a small example for use of eigenMT.
## We provide a genotype matrix and corresponding position file in Matrix-eQTL format.
## This genotype matrix is for chromosome 19 for the EUR373 samples as part of the GEUVADIS study.
## We also provide a sample of the cis-eQTL tests performed for these samples using Matrix-eQTL.
## This sample SNP-gene tests file represents 100 genes.
## The full SNP-gene tests file for chromosomne 19 is also provided for replication of results from our paper. 
## To run eigenMT on the example data, use the following command from the eigenMT directory:

python eigenMT.py --CHROM 19 \
	--QTL qtls.txt \
	--GEN genotypes.txt \
	--POS positions.txt \
	--OUT exampleOut.txt

## Note: this example uses the default settings for window size and variance explained threshold.

####5. Population stratification and other covariates

## We recommend first removing the effects of population stratification 
## (genotype PCs, population or ancestry assignments) and other covariates 
## (age, gender, PEER factors, etc) from the expression matrix. 
## We have shown that by first removing these effects and then performing cis-eQTL 
## calling on the inverse rank normalized residuals, we ensure the conservativeness 
## and accuracy of eigenMT. We will soon provide an R script with example code for 
## how to perform this normalization step prior to a run of eigenMT. 


#### 6. Citation
## When using eigenMT, please cite:
	## Davis, JR, Fresard, L, et al. eigenMT: efficient multiple hypothesis testing correction for eQTLs.




