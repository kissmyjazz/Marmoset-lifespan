library(here)
library(tidyverse)
library(readr)
library(ggstatsplot)

set.seed(3)


df <- readr::read_csv(here("processed_data", "cell_counts.csv"))

df_census <- readr::read_csv(here("processed_data", "marmoset_census.csv")) |> 
  dplyr::select(marmoset = name, sex, age = group) |> 
  dplyr::mutate(sex = factor(sex),
                age = factor(age))

df <- df |> 
  dplyr::mutate(marmoset = factor(marmoset),
                measurement_region = factor(measurement_region),
                celltype = factor(celltype))

df <- df |> left_join(df_census)


# NeuN --------------------------------------------------------------------

df_neun_sum <- df |> 
  dplyr::filter(celltype == "NeuN") |> 
  dplyr::group_by(marmoset, measurement_region, celltype) |> 
  summarise(cell_density = mean(cell_density),
            NeuN_density = mean(NeuN_density),
            NeuN_counts = mean(NeuN_counts),
            NeuN_area = mean(NeuN_area),
            densities_ratio = mean(densities_ratio),
            sex = first(sex),
            age = first(age), .groups = "drop")

# area 46 infragranualar
df_neun_46i_sum <- df_neun_sum |> 
  dplyr::filter(measurement_region == "46-infragranular")


ggbetweenstats(
  data  = df_neun_46i_sum,
  x     = age,
  y     = NeuN_density,
  title = "Distribution of NeuN-ir densities"
)


# Parvalbumin -------------------------------------------------------------

df_pv_sum <- df |> 
  dplyr::filter(celltype == "PV") |> 
  dplyr::group_by(marmoset, measurement_region, celltype) |> 
  summarise(cell_density = mean(cell_density),
            NeuN_density = mean(NeuN_density),
            NeuN_counts = mean(NeuN_counts),
            NeuN_area = mean(NeuN_area),
            densities_ratio = mean(densities_ratio),
            sex = first(sex),
            age = first(age), .groups = "drop")

# area 46 infragranualar
df_pv_46i_sum <- df_pv_sum |> 
  dplyr::filter(measurement_region == "46-infragranular")

# unadjusted densities
ggbetweenstats(
  data  = df_pv_46i_sum,
  x     = age,
  y     = cell_density,
  title = "Distribution of PV-ir relative densities versus NeuN-ir densities"
)

# NeuN-adjusted density ratio
ggbetweenstats(
  data  = df_pv_46i_sum,
  x     = age,
  y     = densities_ratio,
  title = "Distribution of PV-ir relative densities versus NeuN-ir densities"
)
