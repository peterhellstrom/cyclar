# https://r-pkgs.org

library(devtools)

p <- "W:/projects/R/cyclar"
# usethis::create_package(p, check_name = FALSE)

usethis::use_mit_license()

use_git_config(user.name = "peterhellstrom", user.email = "peter.hellstrom@nrm.se")
usethis::use_git()
usethis::use_github()
# GitHub API error (401): Bad credentials

create_github_token()
load_all()

# Must run document() to add export functions to NAMESPACE
document()
install()

test()

# Ignore ----
use_build_ignore(c("data-raw", "development", "examples"))

# Document data:
# https://r-pkgs.org/data.html

install_github("peterhellstrom/cyclar")

## Load package ----
library(cyclar)

## Data sets ----
usethis::use_data_raw()

lynx_elton
timetk::tk_tbl(lynx_elton, rename_index = "year")

rodents <- rodents |>
  janitor::clean_names() |>
  dplyr::rename(year = ar, season = arstid) |>
  dplyr::mutate(
    season_text = dplyr::case_when(
      season == 2 ~ "fall",
      season == 1 ~ "spring",
      TRUE ~ NA_character_
    ),
    .after = season
  )

rodents |>
  dplyr::filter(site == "sf")

rodents_long |>
  dplyr::mutate(species = tolower(species)) |>
  dplyr::rename(year = ar, season = arstid)

voles_kilpis
