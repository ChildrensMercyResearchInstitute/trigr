library(magrittr)
library(pander)
library(yaml)
library(tibble)
#library(glue)

panderOptions("table.split.table", Inf)

figure_path <- function(filename="") {
  file.path(getOption("knitr.figure_dir"), filename)
}

# Format number with commas
pn <- function(i, ...) {
  prettyNum(i, big.mark=",", ...)
}

yesno <- function(booleans) {
  ifelse(booleans == TRUE, "Yes", "No")
}

# Force knitr to stop evaluation when an error is encountered
knitr::opts_chunk$set(error=FALSE)

# Don't show code blocks by default
knitr::opts_chunk$set(echo=FALSE)

# Don't reformat R code
knitr::opts_chunk$set(tidy=FALSE)

# Set up figure defaults
knitr::opts_chunk$set(fig.cap="", fig.width=7, fig.height=5, fig.path=figure_path())

# Create output directory if it doesn't exist
if(!dir.exists(getOption("knitr.figure_dir"))) {
  dir.create(getOption("knitr.figure_dir"))
}

#source("/gpfs/data/home/wacheung/analysis/single_cell/Laurel_scPig/SRC/caching.r")
#source("/gpfs/data/home/wacheung/analysis/single_cell/Laurel_scPig/SRC/parallel.r")
source("SRC/caching.r")
source("SRC/parallel.r")
