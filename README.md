# Single-cell RNA-seq analysis reveals distinct tumor and immunosuppressive T-cell phenotypes in CLL patients treated with ibrutinib.

This repository is an R project holding the code necessary to reproduce all figures from this manuscript.

Processed data is available as an R package from the link included in the manuscript. This analysis project requires the data package to function. The data package contains R objects, documentation, and code used generate the R objects.

## System requirements

1.  R version 4.5.0 or greater
2.  8 GB memory

## Installation

### Using a standard workflow

The data package can be installed from the .tar.gz file like any other R package.

Because of the size of the objects they are no lazy-loaded so you will have to load them manually using

```         
data(cds_main)
```

Because the single cell objects are stored with the matrix on disk, you will have to use monocle3 functions to load them from the package installation directory with something like

```         
load_monocle_objects("</path/to/package/directory/>cll.scrnaseq.datapkg/extdata/cds_main")
```

### Using blaseRtemplates (recommended)

Alternatively, you can load the data package using functions from the blaseRtemplates package. This will handle loading the on-disk single cell objects and will load other objects as promises which will come into memory when and if they are called.

blaseRtemplates is a full computing environment but the project_data function can be used as a standalone function.

To install blaseRtemplates, run:

```         
install.packages('blaseRtemplates', repos = c('https://blaserlab.r-universe.dev', 'https://cloud.r-project.org'))
```

Additional documentation is available at <https://blaserlab.github.io/blaseRtemplates/>.

### Install dependencies

The tsv file library_catalogs/blas02_cll_scrnaseq.tsv lists all R package dependencies for the analysis.  Filtering the status column for status == "active" will list the direct dependencies.  All others aer available at the time of publication and include indirect dependencies and unrelated/unused packages.  You should install the "active" packages and their dependencies using your method of choice.

## Reproducing figures

1.  Edit R/dependencies.R to point the project_data function to the data package.
2.  Source R/dependencies.R
3.  Edit R/configs.R to point to the appropriate output directories.
4.  Source R/configs.R
5.  Open any R file within R/fig_composition.  Make sure source_all_figs is TRUE.  Source that file.
6.  To reproduce individual panels, refer to the files in R/fig_staging.  Each will produce 1 or occasionally a few panels.  The figure composition files list the required figure staging files.

Note:  The main single cell object is cds_main.  This is a monocle3 cell_data_set.  It can be visualized and explored using functions from that package.

## Reproducing tables

1.  See 1-4 in the previous section.
2.  Source R/tables.R

