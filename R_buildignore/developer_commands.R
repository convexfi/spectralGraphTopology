##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type="source")
# Installation from GitHub
devtools::install_github("dppalomar/spectralGraphTopology")
# Installation from CRAN
install.packages("spectralGraphTopology")
# Getting help
library(spectralGraphTopology)
help(package="spectralGraphTopology")
package?spectralGraphTopology
?learnGraphTopology


##
## Developer commands (http://r-pkgs.had.co.nz/)
##
library(devtools)
devtools::load_all()  #or Ctrl-Shift-L
devtools::install()
#devtools::build()  # to generate the installation file

# Documentation
devtools::document()  #to generate all documentation via roxygen
?learnGraphTopology

# Code tests
devtools::test()
#covr::package_coverage()  #coverage of tests


# CRAN check and submission (http://r-pkgs.had.co.nz/release.html)
#  checklist: https://kalimu.github.io/post/checklist-for-r-package-submission-to-cran/
devtools::check()
rcmdcheck::rcmdcheck()
devtools::build()
#devtools::build_win()  #to check under windows
#R CMD build .  # this is to generate tarball
#R CMD check spectralGraphTopology_0.1.0.tar.gz --as-cran  # this is before submission to CRAN
#R CMD install spectralGraphTopology_0.1.0.tar.gz
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html

