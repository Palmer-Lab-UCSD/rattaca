# RATTACA :rat:

:construction: Under Construction :construction:

The RATTACA method aims at predicting rat
traits from low coverage sequencing measurements.
Using these trait predictions, animals with predicted
extreme phenotypes will be selected for study.  Here
we provide helper tools for fitting a Linear
Mixed Model (LMM) with the 'rrBLUP' package to
data, assess the LMM performance under the
stated goals of the RATTACA project, and scripts to
integrate in bioinformatic pipelines.


## Installation

First you'll need to make sure that:

* R >= 4.1 is installed

as well as the following `R` packages

* genio (>=1.1.2)
* rrBLUP (>=4.6.2)


### Cloning the repository

First, simply clone the repository.  Assuming that
the cloned repository is in the current working
directory

```bash
R CMD INSTALL [-l library_path] rattaca
```

Or if you prefer to install from R, you'll need to
build the package.  This can be done from the
command line by:

```bash
R CMD build rattaca
```

which will generate a tar ball
`rattaca-<version_number>.tar.gz` in the
current working directory.  Then start `R`
and in the `R` command line

```R
install.packages('rattaca-<version_number>.tar.gz', repos=NULL)
```

where we have set `repos` to `NULL` as we are
installing from a local package.


### Download a tagged version

On this GitHub page, on the right-hand side, under
`Releases` select `tags`.  Choose a tagged version,
and click on the compressed source code link for
`.tar.gz`.  Given that this downloaded in the
`Downloads` folder

```bash
cd ~/Downloads
R CMD INSTALL [-l library_path] rattaca-<version_number>.tar.gz
```

## Package Features

Here we introduce `rattaca` package features with a series
of examples.

### Example: fitting model to simulation data


Using the `R` command line we can generate arbitrary simulated
biallelic genotypes

```R
library(rattaca)

q_markers <- 100
n_samples <- 210

genotypes <- sample(c(0,1), q_markers * n_samples,
                    replace=TRUE)

genotypes <- matrix(genotype, ncol=q_markers, nrow=n_samples)
``` 

from these data we can generate a simulation closure to 
sample trait values given a trait heritability and 
genotypes

```R
heritability <- 0.4

simulator <- rattaca::gen_sim_closure_r_sq(heritability, genotypes)

trait_data <- simulator(genotypes)
```

Using the genotype and trait pairs, we can fit the data using
the `rattaca` wrapper for rrBLUP's `rrBLUP::mixed.solve`.

```R
out_pars <- rattaca::fit(trait_data,
                         genotypes)
```

Then by using the `ls` function we see that the `out_pars`
list contains 

```R
ls(out_pars)
[1] "beta"  "beta.SE"   "LL"  "pearson_corr"
[5] "r_sq"  "spearman_corr    "u"
[9] "u.SE"  "Ve  "Vu"
```

where the parameter definitions are as follows:

| parameter | description |
| --------- | ------------|
| beta        |   Fixed effect size, i.e. the intercept term |
| beta.SE     |   Standard error of intercept |
| LL          |   log-likelihood of fit |
| pearson_corr    | Pearson correlation between predictions and training data |
| r_sq        |   Coefficient of determination, R^2, between predictions and training data |
| spearman_corr   | Spearman correlation between predictions and training data |
| u   |   BLUP for each marker random effect |
| u.SE | Standard error BLUP for each marker random effect |
| Ve  | Variance of error term |
| Vu  | Variance of random effects |


Note, that fitting is computationally expensive for large amounts
of fitting data.  Test on your system prior to running to estimate
time.


### Loading genotypes and trait data from R
```R
genotypes <- rattaca::load_and_prepare_plink_data("genotypes")
trait_data <- rattaca::load_and_prepare_trait_data("traits.csv",
    "rfid", "mass")
```


## Example of using RATTACA Scripts

`rattaca` scripts can be found in the `inst` directory of this repository.


### Fitting marker BLUPs

To fit the model to data, you'll need to have on hand:

* path to plink bed, bim, and fam files, e.g. `data/genotypes.{bed,bim,fam}`
* path to trait file
    - csv table, rows being animals and columns being trait values
    - headers are required
    - first column should be rfid
* name of trait 


```bash
$ Rscript rat_fit.R \
--plink_genotypes_prefix PATH_TO_GENE\
--
```

### Power analysis





## Contribute

Please clone or fork the repository.  Submit changes as a new
branch and open a `Pull Request`.  Also, please add new unit tests
and run them prior to submitting a recommended change.  Recall that
an R package can be tested by, assuming that the cloned `rattaca`
package directory is in your current working directory,

```bash
R CMD check rattaca
```

This command will:

* print several helpful messages to the standard out,
* create directory `rattaca.Rcheck` with helpful log files,
* build documentation, and 
* run unit-tests.


## Citation

* Johnson, Benjamin B., et al. "RATTACA: Genetic predictions in
Heterogeneous Stock rats offer a new tool for genetic correlation
and experimental design." bioRxiv (2023): 2023-09.

## References

* Endelman, J.B. 2011. Ridge regression and other kernels
for genomic selection with R package rrBLUP. Plant Genome 4:250-255.
beta
beta.SE
LL
pearson_corr
r_sq
spearman_corr
u
u.SE
Ve
Vu
