# https://r-pkgs.org

library(devtools)
library(usethis)

p <- "W:/PROJEKT/R/cyclar"
#usethis::create_package(p, check_name = FALSE)
usethis::create_github_token()

usethis::use_git()
usethis::use_github_links()
usethis::use_mit_license()

use_git_config(user.name = "peterhellstrom", user.email = "peter.hellstrom@nrm.se")


devtools::load_all()

# Must run document() to add export functions to NAMESPACE
devtools::document()
devtools::install()

devtools::test()

# Document data:
# https://r-pkgs.org/data.html

## Load package ----
library(cyclar)

## Data sets ----
lynx_elton

rodents <- rodents |>
  janitor::clean_names() |>
  dplyr::rename(year = ar, season = arstid) |>
  dplyr::mutate(
    season_text = dplyr::case_when(
      season == 2 ~ "fall",
      season == 1 ~ "spring",
      TRUE ~ NA_character_)) |>
  dplyr::relocate(season_text, .after = season)

rodents |>
  dplyr::filter(site == "sf")

rodents_long |>
  dplyr::mutate(species = tolower(species)) |>
  dplyr::rename(year = ar, season = arstid)

voles_kilpis

## Functions ----

