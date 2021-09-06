# cll_scrnaseq

This is an analysis project repository for the BTKi resistance scRNA-seq project.

## Installation

You should download the source data from the link provided or have access to it on the local network.

The source data is a compiled R binary package with the prefix cll.scrnaseq.datapkg followed by version number.

The code in R/dependencies.R will set up an "renv" environment from the lockfile in the project repository.  Basically this copies the required R packages from your library into the local project directory.  Run 

```r
renv::restore()
```

to do this.  If there are missing packages it should prompt you to install.  If running on the blaserlab server, uncommenting and running the block of renv::install() commands should install anything that isn't automatically installed.  

It is highly recommended that you install the package "blaseRtools" and all required dependencies.  

If not running on the blaserlab server this can be installed with 

```r
renv::install("blaserlab/blaseRtools")
```

Next, run

``` r
blaseRtools::bb_renv_datapkg("<path to directory holding the datapkg>")
```

This will install or update the latest version of the data package in your renv environment for your analysis project.  

Load the libraries by running the library() commands.

Finally run

```r
lazyData::requireData("cll.scrnaseq.datapkg")
```

as shown.  This is a workaround to put the data in a hidden environment that can be called like lazy data would normally be called.  ScRNA-seq data packages are too large to load via normal "lazy loading" in R.

If you modify the data during analysis, for example by adding a new metadata column, the modified object will appear in the global environment.  If we like the change we can add it to the next version of the data package.

## Example

For example to create a umap plot:

``` r
requireData("cll.scrnaseq.datapkg")
bb_var_umap(cds_main, var = "partition")
```

