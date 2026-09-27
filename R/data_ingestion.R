library(here)
library(tidyverse)
library(readr)
library(arrow)

path <- here("raw_data") |> normalizePath()

# regex patterns
marmoset_rgx <- "^[^_]+"

# Get all csv files in the folder
data_files_ann <- list.files(path = path, pattern = "*_ann.tsv$",
                             full.names = TRUE, recursive = TRUE)

named_data_files_ann <- purrr::set_names(data_files_ann, nm = str_extract(data_files_ann, 
                                                                          pattern = "/(.+)\\.tsv$"))


df <- named_data_files_ann |>
  map_dfr(read_tsv, skip_empty_rows = TRUE, col_types = "ccddddiiiiiii", name_repair = "universal_quiet") |>
  dplyr::mutate(marmoset = str_extract(Image, marmoset_rgx) |> factor(), .after = Image) |>
  dplyr::mutate(celltype = dplyr::case_when(!is.na(`Num.NPY_SST`) ~ "NPY_SST",
                                          !is.na(`Num.NPY`) ~ "NPY",
                                          !is.na(`Num.SST`) ~ "SST",
                                          !is.na(`Num.PV`) ~ "PV",
                                          !is.na(`Num.NeuN`) ~ "NeuN",
                                          !is.na(`Num.VIP`) ~ "VIP", .ptype = "factor"),
                layer = dplyr::case_when(str_detect(Name, "-supragranular") ~ "supragranular",
                                        str_detect(Name, "-infragranular") ~ "infragranular",
                                        str_detect(Name, "-wm") ~ "wm",
                                        .default = Name),
                brain_region = str_extract(Name, "[^-]+")) |>
  dplyr::rename(measurement_region = Name) |>
  dplyr::select(-all_of(c("Num.NPY_SST", "Num.NPY", "Num.SST", "Num.PV", "Num.NeuN", "Num.VIP"))) |>
  dplyr::distinct() 

# For some marmosets more than 3 sections were processed. In these few cases I
# remove extra sections from the final data set.

exclude <- c("Atticus_W4_S1", "Atticus_W4_S2", "Atticus_W4_S3", "Triton_W5_S2",
             "Usopp_W17_S1", "Usopp_W17_S2", "Usopp_W17_S3", "Walnut_W10_S5_V1missing",
             "Atticus_W4_S1_NeuN", "Atticus_W4_S2_NeuN", "Atticus_W4_S3_NeuN", 
             "Usopp_W17_S1_NeuN", "Usopp_W17_S2_NeuN", "Usopp_W17_S4_NeuN",
             "Cake_W12_S2_NeuN", "Salacia_W9_S2_VIP", "Salacia_W9_S3_VIP",
             "Walnut_W10_S2_VIP", "Walnut_W10_S3_VIP")

# order the file so that NeuN cell counts are first for each section within marmoset 
df2 <- df |> dplyr::arrange(Image) |> 
  dplyr::filter(!Image %in% exclude) |> 
  dplyr::mutate(section = str_extract(Image, "S\\d+.*")) |> 
  dplyr::arrange(marmoset, celltype, section, measurement_region) |> 
  # order sections for each marmoset and cell type and assign them sequential numbers 
  # to match sections across series
  dplyr::group_by(marmoset, celltype) |> 
  dplyr::mutate(section_id = dense_rank(parse_number(section))) |> 
  dplyr::ungroup()


# to calculate cell type to NeuN ration I will make another column where NeuN Values are carried forward

celltype_order <- c("NeuN", "PV", "NPY", "SST", "NPY_SST", "VIP")

df3 <- df2 |> 
  dplyr::arrange(marmoset, section_id, measurement_region, factor(celltype, levels = celltype_order)) |> 
  # area is converted to square mm
  dplyr::mutate(cell_density = (Num.Detections/(Area.µm.2/1e6))) |> 
  dplyr::group_by(marmoset, section_id, measurement_region) |> 
  dplyr::mutate(NeuN_density = first(cell_density), 
                NeuN_counts = first(Num.Detections),
                NeuN_area = first(Area.µm.2)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(densities_ratio = (cell_density/NeuN_density)) |> 
  dplyr::select(-c("Centroid.X.µm", "Centroid.Y.µm", "Perimeter.µm"))


readr::write_csv(df3, here("processed_data", "cell_counts.csv"))
