# RATTACA :rat:

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

First, simply clone the repository.  Assuming that the cloned
repository is in the current working directory

```bash
$ R CMD INSTALL [-l library_path] rattaca
```

Or if you prefer to install from R, you'll need to build the package

```bash
$ R CMD build rattaca
```

which will generate a tar ball `rattaca-<version_number>.tar.gz` in the
current working directory

```bash
$ R
> install.packages('rattaca-<version_number>.tar.gz', repos=NULL)
```


### Download a tagged version

On this GitHub page, on the right-hand side, under `Releases` select `tags`.
Choose a tagged version, and click on the compressed source code link for
`.tar.gz`.  Given that this downloaded in the `Downloads` folder

```bash
$ cd ~/Downloads
$ R CMD INSTALL [-l library_path] rattaca-[version number].tar.gz
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



## Package Feature Example


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
