library(zoo)
library(tidyverse)

# Oikos 59(1):115-120
# From Table 1, p. 116

# 0 = low year, below average
# 1 = peak year, above average
# Regions: sn: South Norway, cn: Central Norway, cb: Combined, nn: North Norway

sn <- "01100110011011100110110011001100110011111010001011101100100110011001100100110011"
cn <- "01100110011001100100110011001100110010001100000001101100100110001001100100010011"
cb <- "01100110011011100110110011001100110011111110001011101100100110011001100100110011"
nn <- "00100110011011000100100011001000110011001100111011100100100011011001100110011011"

b <- map(
  list(sn, cn, cb, nn),
  ~ str_split(.x, "") |>
    unlist() |>
    ts(start = 1871)
) |>
  setNames(c("sn", "cn", "cb", "nn"))

ts.union(b)

x <- ts.union(
  "southeast_norway" = b$sn,
  "central_norway" = b$cn,
  "combined" = b$cb,
  "north_norway" = b$nn
)

plot(x)
plot(as.zoo(x), main = "Wildhagen's data as analyzed by Steen et al 1990")

x_df <- as_tibble(x) |>
  mutate(
    across(1:4, as.integer),
    year = 1871:1950
  ) |>
  relocate(year)

writexl::write_xlsx(
  list(wildhagen = x_df),
  "Steen_et_al_1990_Wildhagen_data.xlsx"
)
