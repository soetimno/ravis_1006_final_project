library(tidyverse)

x1 <- read_csv(
  "TitlingDeforestation/BHM_sample20.csv")

str(x1)

x1 <- as_tibble(x1)

