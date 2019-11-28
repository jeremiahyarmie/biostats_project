library(tidyverse)

hp_meta_raw <- read.csv("hp_gwas_isolates_used.csv")
hp_meta <- as_tibble(select(hp_meta_raw, isolate, gwas_group))