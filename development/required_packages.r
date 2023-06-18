# Check that required packages are installed ----

is.installed <- function(mypkg) mypkg %in% installed.packages()[,1]

install.missing <- function(mypkg, 
                            repos = "https://ftp.acc.umu.se/mirror/CRAN/", 
                            dependencies = TRUE) {
  for (i in 1:length(mypkg)) {
    if (is.installed(mypkg[i]) == FALSE) {
      install.packages(
        pkgs = mypkg[i], 
        lib =.Library, 
        repos = repos, 
        dependencies = dependencies)
    }
  }
}

pkgs <- c(
  "bbmle", 
  "pastecs", 
  "MASS", 
  "mgcv", 
  "splines", 
  "nlme",
  "emdbook", 
  "TeachingDemos", 
  "zoo")

# install.missing(pkgs)

# Load required packages ----
# pkgs <- as.list(c(pkgs))
# lapply(pkgs, require, character.only = TRUE)