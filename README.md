<p align="center">
<img src="static/logo.png" width="200">
</p>

Learns the topology of a class of graphs such as K-component, tree, bipartite,
and regular.

See the vignette for a detailed documentation with several illustrative examples.

The package is based on the paper:

## Installation

```r
# Installation from GitHub
devtools::install_github("dppalomar/spectralGraphTopology")
```

Alternatively, one can install the development version as follows:
```
$ git clone https://github.com/dppalomar/spectralGraphTopology.git
$ cd spectralGraphTopology
$ make build && make install
```


```r
# Getting help
library(spectralGraphTopology)
help(package = "spectralGraphTopology")
package?spectralGraphTopology
?spectralGraphTopology

# Citing this work
citation("spectralGraphTopology")
```

## Tests
To run unit tests, type on your favourite terminal:
```
$ cd spectralGraphTopology
$ Rscript -e "devtools::test()"
```

If any of the tests fail, please open a ticket [here](https://github.com/dppalomar/spectralGraphTopology/issues).

## Vignette
For more detailed information, please check the vignette in
[html](https://rawgit.com/dppalomar/spectralGraphTopology/master/vignettes/SpectralGraphTopology-vignette.html) or
[pdf](https://rawgit.com/dppalomar/spectralGraphTopology/master/vignettes/SpectralGraphTopology-vignette.pdf).


## Links
Package: [GitHub](https://github.com/dppalomar/spectralGraphTopology)

README file: [GitHub-readme](https://rawgit.com/dppalomar/spectralGraphTopology/master/README.html)
