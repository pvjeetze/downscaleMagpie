# Downscaling of MAgPIE Outputs to 0.5 Degree

R package **downscaleMagpie**, version **0.1.0.9001**

[![CRAN status](https://www.r-pkg.org/badges/version/downscaleMagpie)](https://cran.r-project.org/package=downscaleMagpie)     

## Purpose and Functionality

Prepares the outputs of a MAgPIE run for the model-agnostic downscaling kernels of
    luscale and writes the gridded results: readers of the gdx, the configuration and the 0.5 degree
    files of the output folder, the MAgPIE-specific products (land with sub-types and transitions by
    stage, peatland, land conservation, age classes, cropland split, biodiversity intactness) and a
    writer of the output files.


## Installation

For installation of the most recent package version an additional repository has to be added in R:

```r
options(repos = c(CRAN = "@CRAN@", pik = "https://rse.pik-potsdam.de/r/packages"))
```
The additional repository can be made available permanently by adding the line above to a file called `.Rprofile` stored in the home folder of your system (`Sys.glob("~")` in R returns the home directory).

After that the most recent version of the package can be installed using `install.packages`:

```r 
install.packages("downscaleMagpie")
```

Package updates can be installed using `update.packages` (make sure that the additional repository has been added before running that command):

```r 
update.packages()
```

## Questions / Problems

In case of questions / problems please contact Patrick von Jeetze <vjeetze@pik-potsdam.de>.

## Citation

To cite package **downscaleMagpie** in publications use:

von Jeetze P, Sauer P (2026). "downscaleMagpie: Downscaling of MAgPIE Outputs to 0.5 Degree." Version: 0.1.0.9001.

A BibTeX entry for LaTeX users is

 ```latex
@Misc{,
  title = {downscaleMagpie: Downscaling of MAgPIE Outputs to 0.5 Degree},
  author = {Patrick {von Jeetze} and Pascal Sauer},
  date = {2026-10-16},
  year = {2026},
  note = {Version: 0.1.0.9001},
}
```
