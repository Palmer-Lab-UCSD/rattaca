# Understanding rrBLUP, and marker vs. kinship LMM's
#
# By: Robert Vogel
# Date: 2023-08-08

.libPaths("~/.local/R-4.2.1/lib")

library(readr)
library(genio)
library(rrBLUP)

N_RATS <- 100
Q_SNPS <- 250

GENO_PREFIX <- file.path("data", "ldprune")
TRAIT_FILENAME <- file.path("data", "test_phenotypes.csv")

SUBSET_DIR <- "subset"
SUBSET_TRAIT <- file.path(SUBSET_DIR, "trait.csv")
SUBSET_GENO <- file.path(SUBSET_DIR, "geno")
SUBSET_SNP_FILENAME <- file.path(SUBSET_DIR, "snp_ids.txt")
SUBSET_RAT_FILENAME <- file.path(SUBSET_DIR, "rat_ids.txt")


# make data set with 100 rats and 250 snps

if (!dir.exists(SUBSET_DIR))
{
    dir.create(SUBSET_DIR, mode="0750")
}


# To simplify analysis, I randomly subset data and write
# to file.  If all requisite subset files exists then
# continue, else make a new set of subset files.

if (file.exists(SUBSET_TRAIT) &&
    file.exists(paste(SUBSET_GENO, "bim", sep=".")) &&
    file.exists(paste(SUBSET_GENO, "bed", sep=".")) &&
    file.exists(paste(SUBSET_GENO, "fam", sep=".")))
{} else {

    # make snp subset list
    snps <- genio::read_bim(paste(GENO_PREFIX, "bim", sep="."))
    
    write.table(snps$id[sample.int(nrow(snps), size=Q_SNPS, replace=FALSE)],
                file=SUBSET_SNP_FILENAME,
                quote=FALSE,
                row.names=FALSE,
                col.names=FALSE)


    # make rat_id subset list
    trait <- read.csv(TRAIT_FILENAME)
    rat_genos <- genio::read_fam(paste(GENO_PREFIX, "fam", sep="."))

    trait <- merge(trait, rat_genos, by.x="rfid", by.y="id")
    trait <- trait[!apply(is.na(trait), 1, any),]

    if (nrow(trait) > N_RATS){
        trait <- trait[sample.int(nrow(trait), 
                                  size=N_RATS,
                                  replace=FALSE), ]; 
    }

    # to be used for rrBLUP, trait values for rats
    write.table(trait,
                file=SUBSET_TRAIT,
                quote=FALSE,
                row.names=FALSE,
                sep=",")

    # to be used for subsetting plink files
    write.table(trait[, c("fam", "rfid")],
                file=SUBSET_RAT_FILENAME,
                quote=FALSE,
                row.names=FALSE,
                col.names=FALSE)


    # make subset plink files
    system2("plink1.9",
            args=c(paste("--bfile", GENO_PREFIX),
                   paste("--keep", SUBSET_RAT_FILENAME),
                   paste("--extract", SUBSET_SNP_FILENAME),
                   "--make-bed",
                   paste("--out", SUBSET_GENO))
                   )
}


