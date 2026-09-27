library(here)
library(tidyverse)
library(readr)
library(glmmTMB)
library(performance)

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


# Parvalbumin -------------------------------------------------------------

df_pv <- df |> 
  dplyr::filter(celltype == "PV")

m_pv <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset),
  family = betabinomial(link = "logit"),
  data = df_pv
)

m_pv2 <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset) +
    (1 | marmoset:measurement_region),
  family = betabinomial(link = "logit"),
  data = df_pv
)

summary(m_pv)
summary(m_pv2)
compare_performance(m_pv, m_pv2)
