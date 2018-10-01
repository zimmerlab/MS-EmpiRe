# MS-EmpiRe

**M**ass **S**pectrometry analysis using **Empi**rical and **Re**plicate based statistics.

MS-EmpiRe is a R package for quantitative analyses of Mass Spectrometry proteomics data.
It allows highly sensitive and specific identitification of differentially
abundant proteins between different experimental conditions.


## Installation

### Dependencies

MS-EmpiRe requires the R package `Biobase` from [Bioconductor](https://www.bioconductor.org).
`Biobase` can be installed from the R command line using the following
commands:

	source("https://bioconductor.org/biocLite.R")
	biocLite("Biobase")


### Installing MS-EmpiRe

You can install MS-EmpiRe directly from github using the R package `devtools`.

Installing `devtools`:

	install.packages("devtools")

Loading `devtools`:

	library(devtools)

Installing MS-EmpiRe:

	install_github("zimmerlab/MS-EmpiRe")

Loading MS-EmpiRe:

	library(msEmpiRe)

## License
MS-EmpiRe is released under the GNU Affero General Public License. See LICENSE for further details.
