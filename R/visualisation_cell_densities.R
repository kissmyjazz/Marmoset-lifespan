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
df_neun_46s_sum <- df_neun_sum |> 
  dplyr::filter(measurement_region == "46-supragranular")


ggbetweenstats(
  data  = df_neun_46s_sum,
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

# area 46 infragranualar
df_pv_46s_sum <- df_pv_sum |> 
  dplyr::filter(measurement_region == "46-supragranular")

# unadjusted densities
ggbetweenstats(
  data  = df_pv_46s_sum,
  x     = age,
  y     = cell_density,
  title = "Distribution of PV-ir densities"
)

# NeuN-adjusted density ratio
ggbetweenstats(
  data  = df_pv_46s_sum,
  x     = age,
  y     = densities_ratio,
  title = "Distribution of PV-ir relative densities adjusted for NeuN-ir densities"
)


# Somatostatin ------------------------------------------------------------

df_sst_sum <- df |> 
  dplyr::filter(celltype == "SST") |> 
  dplyr::group_by(marmoset, measurement_region, celltype) |> 
  summarise(cell_density = mean(cell_density),
            NeuN_density = mean(NeuN_density),
            NeuN_counts = mean(NeuN_counts),
            NeuN_area = mean(NeuN_area),
            densities_ratio = mean(densities_ratio),
            sex = first(sex),
            age = first(age), .groups = "drop")

# area 46 supragranualar
df_sst_46s_sum <- df_sst_sum |> 
  dplyr::filter(measurement_region == "46-supragranular")

# unadjusted densities
ggbetweenstats(
  data  = df_sst_46s_sum,
  x     = age,
  y     = cell_density,
  title = "Distribution of SST-ir densities"
)


# Neuropeptide Y ----------------------------------------------------------

df_npy_sum <- df |> 
  dplyr::filter(celltype == "NPY") |> 
  dplyr::group_by(marmoset, measurement_region, celltype) |> 
  summarise(cell_density = mean(cell_density),
            NeuN_density = mean(NeuN_density),
            NeuN_counts = mean(NeuN_counts),
            NeuN_area = mean(NeuN_area),
            densities_ratio = mean(densities_ratio),
            sex = first(sex),
            age = first(age), .groups = "drop")

# 8Av white matter
df_npy_8av_wm_sum <- df_npy_sum |> 
  dplyr::filter(measurement_region == "8Av-wm")

# unadjusted densities
ggbetweenstats(
  data  = df_npy_8av_wm_sum,
  x     = age,
  y     = cell_density,
  title = "Distribution of NPY-ir densities"
)


# Somatostatin/neuropeptide Y  --------------------------------------------

df_sst_npy_sum <- df |> 
  dplyr::filter(celltype == "NPY_SST") |> 
  dplyr::group_by(marmoset, measurement_region, celltype) |> 
  summarise(cell_density = mean(cell_density),
            NeuN_density = mean(NeuN_density),
            NeuN_counts = mean(NeuN_counts),
            NeuN_area = mean(NeuN_area),
            densities_ratio = mean(densities_ratio),
            sex = first(sex),
            age = first(age), .groups = "drop")

# 8Av white matter
df_sst_npy_v1i_sum <- df_sst_npy_sum |> 
  dplyr::filter(measurement_region == "V1-infragranular")

# unadjusted densities
ggbetweenstats(
  data  = df_sst_npy_v1i_sum,
  x     = age,
  y     = cell_density,
  title = "Distribution of SST+/NPY+-ir densities"
)


# VIP ---------------------------------------------------------------------

df_vip_sum <- df |> 
  dplyr::filter(celltype == "VIP") |> 
  dplyr::group_by(marmoset, measurement_region, celltype) |> 
  summarise(cell_density = mean(cell_density),
            NeuN_density = mean(NeuN_density),
            NeuN_counts = mean(NeuN_counts),
            NeuN_area = mean(NeuN_area),
            densities_ratio = mean(densities_ratio),
            sex = first(sex),
            age = first(age), .groups = "drop")

# 8Av white matter
df_vip_v1s_sum <- df_vip_sum |> 
  dplyr::filter(measurement_region == "V1-supragranular")

# unadjusted densities
ggbetweenstats(
  data  = df_vip_v1s_sum,
  x     = age,
  y     = cell_density,
  title = "Distribution of VIP-ir densities"
)

