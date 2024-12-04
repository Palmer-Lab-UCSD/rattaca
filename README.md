# RATTACA :rat:

:construction: Under Construction :construction:

The RATTACA method uses SNP genotypes produced from low-coverage DNA sequencing 
to predict phenotypes in HS rats, then uses these trait predictions to assign 
for study animals with predicted extreme phenotypes. This package provides 
helper tools for fitting linear mixed models (LMM) for prediction, and assessing 
LMM performance under the goals of the [RATTACA](https://ratgenes.org/rattaca/) 
project.  

Predictions are calculated by estimating marker effects using the 'RR-BLUP' 
method of the [rrBLUP package](https://cran.r-project.org/web/packages/rrBLUP/index.html). 
The primary output of this package - a set of phenotype predictions - can then 
be used to assign rats using the [RATTACA assignment package](https://github.com/Palmer-Lab-UCSD/rattaca_assignment).


## Installation

RATTACA predictions are conducted in `R`. You will need to install R >= 4.1 
with the following `R` packages:

* genio (>=1.1.2) to read/write PLINK genotype data
* rrBLUP (>=4.6.2) to fit prediction models
* roxygen2 (>=2.7.0) (optional: to build package documentaion)
* viridis (>=0.6.2) (for optional plotting)
* scales (>=1.2.0) (for optional plotting)

First, simply clone the repository:
```bash
git clone git@github.com:Palmer-Lab-UCSD/rattaca.git
```

Optionally, to build the package documentation (which allows dynamic generation
of RATTACA documentation in `R`). Assuming the cloned repository is in the 
current working directory, start `R` and from the `R` command line:

```R
library(roxygen2)
roxygenize()
```

Then install from the command line. Assuming the cloned repository is in the 
current working directory:

```bash
R CMD INSTALL [-l library_path] rattaca
```

Alternatively, to install from within R, you'll first need to
build the package. First, from the command line:

```bash
R CMD build rattaca
```

... which will generate a tar ball `rattaca-<version_number>.tar.gz` in the
current working directory.  Then start `R`
and in the `R` command line:

```R
install.packages('rattaca-<version_number>.tar.gz', repos=NULL)
```

... ensuring `repos` is set to `NULL` to install from a local package.


### Download a tagged version

On this GitHub page, on the right-hand side, under
`Releases` select `tags`.  Choose a tagged version,
and click on the compressed source code link for
`.tar.gz`.  Given that this downloaded in the
`Downloads` folder, you can install or build as above:

```bash
cd ~/Downloads
R CMD INSTALL [-l library_path] rattaca-<tagged_version_number>.tar.gz
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

Suppose that our current working direcotry contained a directory
`data` for which the plink and trait data files are located,

```bash
-data/
    |- genotype.bed
    |- genotype.bim
    |- genotype.fam
    |- traits.csv
```

Then to load the genotypes and trait data for `mass` where
`rfid` is the sample id column,

```R
genotypes <- rattaca::load_and_prepare_plink_data("data/genotypes")
trait <- rattaca::load_and_prepare_trait_data("data/traits.csv",
                "rfid", "mass")
```


## RATTACA Scripts

`rattaca` scripts can be found in the `inst` directory of this repository.


### Example: fitting marker BLUPs

To fit the model to data, you'll need to have on hand:

* genotype files
* trait measurements as csv
* `rat_fit.R` script

Note that the animal identifies in the genotype
files must match that of the trait measurements.  Let us
assume that the current working directory


```bash
- rat_fit.R
- data/
    |- genotype.bed
    |- genotype.bim
    |- genotype.fam
    |- traits.csv
```

Then to run the script

```bash
$ Rscript rat_fit.R \
--plink_genotypes_prefix 'data/genotype' \
--trait_file 'data/traits.csv' \
--trait_rat_id_colname 'rfid' \
--trait 'mass' \
--out_dir 'results'
```

which will print the results to file `results/mass.bpar`.

### Example: performing power analysis


## Output specification

The `.bpar` file contains records separated by line.  Records
come in three flavors and are distinguishable by line prefex:

| Prefix | Description |
| -------| ----------- |
| `##`   | Meta data containing key value pairs delimited by `=`. |
| `#`    | Header, defines the values for each variant record, comma delimited. |
| No prefix | Variant record, comma delimited. |

The header includes

| Column Name | Description |
| ----------- | ----------- |
| var_id | As specified in |
| random_effect | BLUP of the random effect for variant produced by `rrBLUP:mixed.solve` |
| random_effect_se | Standard error of BLUP produced by `rrBLUP::mixed.solve` |


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
