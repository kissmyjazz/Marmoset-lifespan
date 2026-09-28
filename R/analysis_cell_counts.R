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


# NeuN --------------------------------------------------------------------

df_neun <- df |> 
  dplyr::filter(celltype == "NeuN")

m_neun <- glmmTMB(
  Num.Detections ~
    age * measurement_region +
    offset(log(Area.µm.2)) +
    (1 | marmoset),
  family = nbinom2(link = "log"),
  data = df_neun
)

m_neun2 <- glmmTMB(
  Num.Detections ~
    age * measurement_region +
    offset(log(Area.µm.2)) +
    (1 | marmoset) +
    (1 | marmoset:measurement_region),
  family = nbinom2(link = "log"),
  data = df_neun
)

m_neun3 <- glmmTMB(
  Num.Detections ~
    age * measurement_region +
    offset(log(Area.µm.2)) +
    (1 | marmoset),
  family = nbinom1(link = "log"),
  data = df_neun
)

m_neun4 <- glmmTMB(
  Num.Detections ~
    age * measurement_region +
    offset(log(Area.µm.2)) +
    (1 | marmoset) +
    (1 | marmoset:measurement_region),
  family = nbinom1(link = "log"),
  data = df_neun
)

m_neun5 <- glmmTMB(
  cell_density ~
    age * measurement_region +
    (1 | marmoset),
  family = Gamma(link = "log"),
  data = df_neun
)

m_neun6 <- glmmTMB(
  cell_density ~
    age * measurement_region +
    (1 | marmoset) +
    (1 | marmoset:measurement_region),
  family = Gamma(link = "log"),
  data = df_neun
)

compare_performance(m_neun, m_neun2, m_neun3, m_neun4, m_neun5, m_neun6)

summary(m_neun4)
summary(m_neun6)

# change the reference group to compare adult with old marmosets

df_neun_old <- df_neun

df_neun_old$age <- relevel(
  factor(df_neun$age),
  ref = "old"
)

m_neun4_old <- update(
  m_neun4,
  data = df_neun_old
)

summary(m_neun4_old)

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

# Somatostatin ------------------------------------------------------------

df_sst <- df |> 
  dplyr::filter(celltype == "SST")

m_sst <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset),
  family = betabinomial(link = "logit"),
  data = df_sst
)

m_sst2 <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset) +
    (1 | marmoset:measurement_region),
  family = betabinomial(link = "logit"),
  data = df_sst
)

summary(m_sst)
summary(m_sst2)
compare_performance(m_sst, m_sst2)

df_sst_old <- df_sst

df_sst_old$age <- relevel(
  factor(df_sst$age),
  ref = "old"
)

m_sst2_old <- update(
  m_sst2,
  data = df_sst_old
)

summary(m_sst2_old)

# Neuropeptide Y ----------------------------------------------------------

df_npy <- df |> 
  dplyr::filter(celltype == "NPY")

m_npy <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset),
  family = betabinomial(link = "logit"),
  data = df_npy
)

m_npy2 <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset) +
    (1 | marmoset:measurement_region),
  family = betabinomial(link = "logit"),
  data = df_npy
)

summary(m_npy)
summary(m_npy2)
compare_performance(m_npy, m_npy2)

# Somatostatin/neuropeptide Y ---------------------------------------------

df_sst_npy <- df |> 
  dplyr::filter(celltype == "NPY_SST")

m_sst_npy <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset),
  family = betabinomial(link = "logit"),
  data = df_sst_npy
)

m_sst_npy2 <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset) +
    (1 | marmoset:measurement_region),
  family = betabinomial(link = "logit"),
  data = df_sst_npy
)

summary(m_sst_npy)
summary(m_sst_npy2)
compare_performance(m_sst_npy, m_sst_npy2)

# Vasoactive intestinal peptide -------------------------------------------

df_vip <- df |> 
  dplyr::filter(celltype == "VIP")

m_vip <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset),
  family = betabinomial(link = "logit"),
  data = df_vip
)

m_vip2 <- glmmTMB(
  cbind(Num.Detections, NeuN_counts) ~
    age * measurement_region +
    offset(log(Area.µm.2 / NeuN_area)) +
    (1 | marmoset) +
    (1 | marmoset:measurement_region),
  family = betabinomial(link = "logit"),
  data = df_vip
)

summary(m_vip)
summary(m_vip2)
compare_performance(m_vip, m_vip2)

df_vip_old <- df_vip

df_vip_old$age <- relevel(
  factor(df_vip$age),
  ref = "old"
)

m_vip2_old <- update(
  m_vip2,
  data = df_vip_old
)

summary(m_vip2_old)
