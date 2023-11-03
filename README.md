# RATTACA

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



## Example of using RATTACA Scripts


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


## Citation

* Johnson, Benjamin B., et al. "RATTACA: Genetic predictions in
Heterogeneous Stock rats offer a new tool for genetic correlation
and experimental design." bioRxiv (2023): 2023-09.

## References

* Endelman, J.B. 2011. Ridge regression and other kernels
for genomic selection with R package rrBLUP. Plant Genome 4:250-255.
