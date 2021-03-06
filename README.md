# tess3r

This has been hacked by Eric C. Anderson to provide a function that lets you pull
out the raster stack underlying the map plots in tess3r.  That function is called
`tess3Q_map_rasters()`.  Enjoy!

tess3r is an R package for estimating and visualizing spatial population
structure based on geographically constrained non-negative matrix
factorization and population genetics. tess3r has fast and efficient
algorithms for estimating ancestry coefficients and for running genome
scans for selection (see [Overview](https://bioshock38.github.io/TESS3_encho_sen/articles/main-vignette.html)).


## Installation

To install this fork do this:

Install the latest version from github (requires [devtools](https://github.com/hadley/devtools)):
```R
# install.packages("devtools")
devtools::install_github("eriqande/TESS3_encho_sen")
```

## References

- Kevin Caye, Timo Deist, Helena Martins, Olivier Michel, Olivier Francois. TESS3: fast inference of spatial population structure and genome scans for selection. Molecular Ecology Resources, Blackwell, 2015, [<10.1111/1755-0998.12471>](http://dx.doi.org/10.1111/1755-0998.12471). [<hal-01222555>](https://hal.archives-ouvertes.fr/hal-01222555)

- Kevin Caye, Flora Jay, Olivier Michel, Olivier Francois. Fast Inference of Individual Admixture Coefficients Using Geographic Data. bioRxiv, 2016, doi: [http://dx.doi.org/10.1101/080291](http://dx.doi.org/10.1101/080291)

- Helena Martins, Kevin Caye, Keurcien Luu, Michael GB Blum, Olivier Francois. Identifying outlier loci in admixed and in continuous populations using ancestral population differentiation statistics. bioRxiv, 2016, doi: [http://dx.doi.org/10.1101/054585](http://dx.doi.org/10.1101/054585)
