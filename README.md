# RATTACA :rat:

:construction: Under Construction :construction:

test

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
R CMD INSTALL [-l library_path] rattaca
```

Or if you prefer to install from R, you'll need to build the package.  This
can be done from the command line by:

```bash
R CMD build rattaca
```

which will generate a tar ball `rattaca-<version_number>.tar.gz` in the
current working directory.  Then start `R` and in the `R` command line

```R
install.packages('rattaca-<version_number>.tar.gz', repos=NULL)
```

where we have set `repos` to `NULL` as we are installing from a
local package.


### Download a tagged version

On this GitHub page, on the right-hand side, under `Releases` select `tags`.
Choose a tagged version, and click on the compressed source code link for
`.tar.gz`.  Given that this downloaded in the `Downloads` folder

```bash
cd ~/Downloads
R CMD INSTALL [-l library_path] rattaca-<version_number>.tar.gz
```

## Package Features

### Fitting model to data

We'll need:

* plink genotype files bed, bim, and fam. e.g. genotypes.{bed,bim,fam}
* table of trait values in csv format, *Note* a header is required. e.g.

rfid | project_id | mass | bmi | color |
-----------------------------------
00383402 | new_proj | 423.21 | 24. | brown |
03324225 | diff_proj | 232.35 | 25. | brown|


Using the `R` command line we can read in the data by

```R
library(rattaca)

genotypes <- rattaca::load_and_prepare_plink_data("genotypes")
trait_data <- rattaca::load_and_prepare_trait_data("traits.csv",
    "rfid", "mass")
``` 

where when loading the trait data we referenced the column with
animal id's and that of the measurements we are interested in.

Next we need to make sure that there exists a genotype and
trait measurement for each sample.

```R
rat_ids <- intersect(rownames(trait_data), rownames(genotypes))
genotypes <- genotypes[rat_ids,]
trait_data <- trait_data[rat_ids,,drop=FALSE]

out_pars <- rattaca::fit(trait_data[,"mass"],
                         genotypes)
```

Note, that fitting is very expensive for large amount of SNPs and
animals.  If you have more than 1,000 SNPs, you'll want to submit
a job to HPC.



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
