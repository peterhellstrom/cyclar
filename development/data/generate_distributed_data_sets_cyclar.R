# Data sets ----

## (Parts of the) Kilpisjärvi vole data set 1949-1997/2007 ----

# Source:
# Hansen, T. F., Stenseth, N. C. & Henttonen, H. (1999) Multiannual vole cycles and population regulation during long winters: an analysis of seasonal density dependence. American Naturalist, 154, 129-139.
# Figure 1, p. 132

# Hansen et al covered 1949-1997 (spring). For 1997 (fall)-2007 (fall) data was provided by Heikki Henttonen (pers. comm).

# 3 variables:
# av = all voles, cr = grey-sided voles, ov = other voles
# 2 points per year, starting at 1=spring, 2=fall

voles_kilpis <-
  structure(c(NA, 20, NA, 15.7, NA, 0.4, NA, 1, 0.2, 5.7, 1.9,
              10.8, 12.9, 15.1, 2.4, 5.6, 0.2, 0.6, 4, 3.1, 5.9, 9.9, 5.9,
              6.3, 0.1, 0.5, 0.2, 3, 2.6, 12.6, 8.8, 13.1, 0.7, 0.4, 0, 1.5,
              0.7, 1.5, 0.4, 5.3, 5.3, 13.9, 6.8, 3.4, 0, 0.1, 0, 1, 3.2, 15.1,
              14, 14.7, 8.2, 3.8, 0, 0.5, 2.8, 11.1, 9.9, 18.2, 0, 0, 0, 1.9,
              2.1, 17.1, 18.1, 24.1, 15.3, 26.3, 0.7, 0.3, 0.2, 1, 1.5, 15.5,
              10.7, 16, 3, 6.4, 0.6, 0.8, 0.1, 2.5, 0.1, 2.4, 1.2, 7.1, 1.8,
              3.6, 1.1, 1.6, 0.2, 0.3, 0.2, 0.5, 0.5, 3.6, 1, 1.6, 0, 2, 0.4,
              3.7, 3.5, 15.3, 10.3, 10.2, 0.2, 0.4, 0.1, 2.2, 1, 3, 1, 4.8,
              2, 7.3, NA, 28.3, NA, 22.1, NA, 0.5, NA, 1.3, 0.4, 6.9, 2.1,
              13.8, 13.4, 23.3, 2.4, 7.4, 0.4, 0.9, 4.2, 3.8, 7.1, 12.1, 15.4,
              18.4, 0.5, 1, 0.2, 3.2, 5.2, 14.1, 16.6, 15.6, 1, 0.5, 0, 1.8,
              1.2, 2.1, 1.2, 7.8, 5.8, 38.1, 8.2, 9.8, 0.1, 0.8, 0.5, 4.6,
              6.1, 19.1, 23, 17.2, 9.7, 4.6, 0.1, 0.6, 3.6, 15.2, 12.5, 32.2,
              0, 0.3, 0.4, 2.1, 3.1, 19.9, 21.2, 27, 16.2, 27.1, 1, 0.4, 0.3,
              1.7, 1.8, 19.2, 12.2, 18.9, 3.9, 7.5, 1.2, 1.6, 0.2, 10.3, 2,
              15.3, 5.3, 16.7, 4.4, 12.4, 1.8, 2.7, 0.3, 6.5, 1.5, 8.2, 6.2,
              23, 4.8, 12, 0.4, 12, 2.9, 15.4, 7.9, 30.2, 23.2, 23.7, 0.5,
              6.5, 0.8, 7.5, 4, 10.1, 2.5, 12.2, 5.5, 20.7, NA, 8.3, NA, 6.4,
              NA, 0.1, NA, 0.3, 0.2, 1.2, 0.2, 3, 0.5, 8.2, 0, 1.8, 0.2, 0.3,
              0.2, 0.7, 1.2, 2.2, 9.5, 12.1, 0.4, 0.5, 0, 0.2, 2.6, 1.5, 7.8,
              2.5, 0.3, 0.1, 0, 0.3, 0.5, 0.6, 0.8, 2.5, 0.5, 24.2, 1.4, 6.4,
              0.1, 0.7, 0.5, 3.6, 2.9, 4, 9, 2.5, 1.5, 0.8, 0.1, 0.1, 0.8,
              4.1, 2.6, 14, 0, 0.3, 0.4, 0.2, 1, 2.8, 3.1, 2.9, 0.9,
              0.8, 0.3, 0.1, 0.1, 0.7, 0.3, 3.7, 1.5, 2.9, 0.9,
              1.1, 0.6, 0.8, 0.1, 7.8, 1.9, 12.9, 4.1, 9.6, 2.6, 8.8, 0.7,
              1.1, 0.1, 6.2, 1.3, 7.7, 5.7, 19.4, 3.8, 10.4, 0.4, 10, 2.5,
              11.7, 4.4, 14.9, 12.9, 13.5, 0.3, 6.1, 0.7, 5.3, 3, 7.1, 1.5,
              7.4, 3.5, 13.4), .Dim = c(118L, 3L), .Dimnames = list(NULL, c("cr",
                                                                            "av", "ov")), .Tsp = c(1949, 2007.5, 2), class = c("mts", "ts"
                                                                            ))

# usethis::use_data(voles_kilpis)

voles_kilpis <- as_tibble(voles_kilpis) %>%
  mutate(
    time_var = as.numeric(time(kilpis_s)),
    year = as.integer(time_var),
    season = case_when(
      (time_var - year) == 0 ~ "spring",
      (time_var - year) == 0.5 ~ "fall",
      TRUE ~ NA_character_)) %>%
  select(-time_var)

voles_kilpis

voles_kilpis %>%
  pivot_wider(names_from = season, values_from = cr:ov)


## Soay sheep ----
# Berryman, A., and Lima, M. 2006. Am. Nat. 168: 784-795. (Appendix, table A1)
# See also Berryman, A., and Lima, M. 2007. Ecology 88(8): 2121-2123.

soay_sheep <- stats::ts(cbind(
	N = c(1296,710,1038,1447,694,889,1449,957,1284,1520,1176,1826,1751,1968,933,1409,1889,907,1568,1996,1365),
	R = c(-0.602,0.38,0.332,-0.735,0.248,0.488,-0.415,0.294,0.169,-0.257,0.44,-0.042,0.117,-0.746,0.412,0.293,-0.734,0.547,0.241,-0.38,NA),
	NAO = c(-0.02,-0.64,-0.17,2.61,2.08,0.63,1.95,1.92,0.86,1.36,-1.43,0.8,0.04,0.68,1.72,-1.03,0.71,0.12,-0.08,0.48,NA),
	NDVI = c(6,6.07,5.78,5.92,5.14,5.23,5.29,4.84,4.84,4.78,5.68,6.25,6.34,5.85,6.09,5.71,3.72,3.85,3.73,NA,NA)),
	start = 1985, deltat = 1)

usethis::use_data(soay_sheep)

as_tibble(soay_sheep) %>%
  mutate(Year = as.numeric(time(soay_sheep))) %>%
  relocate(Year)

## Arctic fox, Sweden 1974-2010 ----
# Angerbjörn, A., Tannerfeldt, M., Bjärvall, A., Ericson, M., From, J. & Norén, E. (1995) Dynamics of the arctic fox population in Sweden. Annales Zoologici Fennici, 32, 55-68.
# Angerbjörn (personal communication)
arctic_fox <- stats::ts(c(179, 59, 0, 23, 175, 0, 30, 155, 183, 4, 25, 73, 50, 49, 28, 66, 90,
	60, 42, 23, 31, 12, 25, 0, 12, 3, 4, 18, 8, 2, 28, 52, 6, 48, 66, 4, 58), start = 1974, deltat = 1)

usethis::use_data(arctic_fox)
