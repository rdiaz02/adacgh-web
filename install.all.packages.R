
## For LAM first install Rmpi doing

## whateverpath R CMD INSTALL


## source this in R, and things should work

source("http://bioconductor.org/biocLite.R")
my.repos <- c("http://stat.ethz.ch/CRAN", biocinstallRepos())


## or http://cran.ch.r-project.org/src/contrib/

## we will change the call to ADaCGH; for now, that helps to
## get the whole set of needed packages

install.packages("Rmpi") ## if not already installed. Or compile against appropriate lib (LAM, whatever)
install.packages(c("ff", "multicore", "fork", "ADaCGH", "snowfall", "papply", "rsprng", "rlecuyer", "snapCGH"), repos = my.repos)


install.packages(c("Rcmdr", "JGR", "Deducer", "getopt", "lme4", "multcomp")) ## for classes and myself

## when using a previously installed dir, such as binary, might need to use:
# update.packages(repos = my.repos, ask = FALSE)
