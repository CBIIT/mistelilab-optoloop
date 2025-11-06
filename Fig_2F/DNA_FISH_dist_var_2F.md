# DNA FISH spot distance & cell-to-cell variability calculations
Martin Stortz/Gianluca Pegoraro
November 6, 2025

### Load packages

### User variables input and settings specification

``` r
GLOB_C <- "*nuclei_information_well*" # Pattern for Single Cell data files
GLOB_S <- "*spots_locations_well*" # Pattern for Single Spot data files

XY_RES <- 0.108 # pixel size in um
```

### Metadata: Experimental conditions and plate layout

Provide experimental conditions and plate layout

``` r
#|
# Provide experimental conditions in desired order:
light_levs <- c("clone1-Dark", "clone1-Light", "clone3-Dark", "clone3-Light", "clone7-Dark", "clone7-Light")

plate_layout <- tibble(
  row = c(8L, 8L, 9L, 9L, 10L, 10L),
  column = c(6L, 20L, 6L, 20L, 6L, 20L),
  well_index = c("H6", "H20", "I6", "I20", "J6", "J20"),
  light = c("clone1-Dark", "clone1-Light", "clone3-Dark", "clone3-Light", "clone7-Dark", "clone7-Light")
) |>
  mutate(light = factor(light, levels = light_levs))
```

### Read cells and spots data

HITIPS output text files must be in the `data` directory and the file
names must match the GLOB pattern.

Download data if needed

### Read cells data

``` r
cell_tbl <- dir_info("data",
  recurse = T,
  glob = GLOB_C
) |>
  filter(size > "1K") |>
  rowwise() |>
  reframe(fread(path), file_name = as.character(path)) |>
  mutate(well_index = paste0(LETTERS[row], column)) |>
  select(file_name, column, row, well_index,
    field_index, time_point, cell_index,
    area, solidity,
    x = `centroid-0`,
    y = `centroid-1`
  ) |>
  arrange(
    column, row, field_index,
    time_point, cell_index
  ) |>
  left_join(plate_layout, by = c("well_index", "row", "column"))

cell_tbl <- cell_tbl[!is.na(cell_tbl$light), ]
```

Read and process the spot level data. Convert the x, y and z coordinates
to microns (They were originally calculated in pixels) and filter spots
that do not belong to any nucleus. These have a `cell_index` value of
`0`.

### Read spots data

``` r
spots_tbl <- dir_info("data",
  recurse = T,
  glob = GLOB_S
) |>
  filter(size > "1K") |>
  rowwise() |>
  reframe(fread(path), file_name = as.character(path)) |>
  mutate(well_index = paste0(LETTERS[row], column)) |>
  select(file_name, column, row,
    well_index, field_index, time_point,
    cell_index, channel,
    spot_index = V1,
    x = x_location, y = y_location
  ) |>
  # Convert xy distances from pixels to microns
  mutate(across(x:y, list(mic = ~ .x * XY_RES))) |>
  # Filter out spots not in any nucleus
  filter(cell_index != 0) |>
  arrange(
    column, row, field_index,
    time_point, cell_index, channel
  ) |>
  left_join(plate_layout, by = c("well_index", "row", "column"))

spots_tbl <- spots_tbl[!is.na(spots_tbl$light), ]
```

### Filter cells by size and shape

``` r
cell_filter_tbl <- cell_tbl |>
  filter(area >= 100, solidity >= 0.875)
```

### Spot number per cell calculations

``` r
cell_n_spots_tbl <- spots_tbl |>
  group_by(
    column, row, well_index, light, field_index,
    time_point, cell_index, channel
  ) |>
  count() |>
  pivot_wider(
    names_from = channel,
    names_prefix = "n_spots_channel_",
    values_from = n,
    values_fill = 0
  ) |>
  ungroup()
```

### Filter cells by spot number

``` r
#|
spot_filt_criterion <- quote(n_spots_channel_2 > 1 & n_spots_channel_3 >= n_spots_channel_2) # set criterion for filtering

cell_n_spots_filt_tbl <- cell_n_spots_tbl |>
  semi_join(cell_filter_tbl, by = c(
    "row", "column", "well_index", "light", "field_index",
    "time_point", "cell_index"
  )) |>
  filter(!!spot_filt_criterion)
```

### Filter out spots according to previous cell filtering by size, shape and spot number

``` r
spot_filt_tbl <- spots_tbl |>
  semi_join(cell_n_spots_filt_tbl,
    by = c(
      "column", "row", "field_index", "light",
      "time_point", "cell_index"
    )
  ) |>
  as.data.table()
```

### Calculate distances between all possible spot pairs for each cell

``` r
gf_dist_2D <- spot_filt_tbl[channel %in% 2:3, # Only  Red and Green spots
  reshape2::melt(
    SpatialTools::dist2(
      as.matrix(.SD[
        channel == 2,
        list(x_mic, y_mic)
      ]),
      as.matrix(.SD[
        channel == 3,
        list(x_mic, y_mic)
      ])
    ),
    value.name = "gf_dist"
  ),
  by = list(
    row,
    column,
    well_index,
    light,
    field_index,
    time_point,
    cell_index
  )
]

setnames(gf_dist_2D, c("Var1", "Var2"), c("g_index", "f_index"))
```

### Find minimum spot distance Minimum distance is calculated based on channel 3 spots as reference

``` r
setkey(
  gf_dist_2D, row, column, well_index,
  field_index, time_point, cell_index, f_index
)

gf_dist_min_2D <- gf_dist_2D[, .SD[which.min(gf_dist), ],
  by = key(gf_dist_2D)
] |>
  left_join(select(cell_tbl, column:well_index, light, field_index:area),
    by = c(
      "column", "row",
      "well_index", "light", "field_index",
      "time_point", "cell_index"
    )
  ) |>
  # Filter out distances > 2 um
  filter(gf_dist <= 2)
```

### Analysis of cell-to-cell contact variability

Calculate observed distributions

``` r
# Keep only cells with 2 alleles per cell
gf_dist_pairwise <- gf_dist_min_2D |>
  group_by(light, field_index, time_point, cell_index) |>
  filter(n() == 2) |>
  ungroup()

# Pivot distances so each cell has gf_dist_1 and gf_dist_2
gf_dist_wide <- gf_dist_pairwise |>
  group_by(light, field_index, time_point, cell_index) |>
  mutate(spot_id = paste0("gf_dist_", row_number())) |>
  select(light, field_index, time_point, cell_index, spot_id, gf_dist) |>
  pivot_wider(names_from = spot_id, values_from = gf_dist) |>
  drop_na(gf_dist_1, gf_dist_2) # Ensure no NAs

# Classify cells according to distances
gf_dist_classified <- gf_dist_wide |>
  mutate(
    category = case_when(
      gf_dist_1 < 0.27 & gf_dist_2 < 0.27 ~ "Both close",
      gf_dist_1 < 0.27 | gf_dist_2 < 0.27 ~ "One close",
      TRUE ~ "Both distant"
    )
  )
observed_distr <- gf_dist_classified |>
  group_by(light, category) |>
  summarise(n_observed = n(), .groups = "drop") |>
  group_by(light) |>
  mutate(fraction_observed = n_observed / sum(n_observed))
```

Calculate expected distributions assuming complete allele independence

``` r
# Count all distances per condition
allele_freqs <- gf_dist_wide |>
  pivot_longer(cols = c(gf_dist_1, gf_dist_2), names_to = "allele", values_to = "distance") |>
  group_by(light) |>
  summarise(
    n_alleles = n(),
    n_close = sum(distance < 0.27),
    p_close = n_close / n_alleles,
    .groups = "drop"
  )

# Get expected counts
expected_distr <- allele_freqs |>
  mutate(
    both_close = p_close^2,
    one_close = 2 * p_close * (1 - p_close),
    both_distant = (1 - p_close)^2
  ) |>
  pivot_longer(
    cols = c(both_close, one_close, both_distant),
    names_to = "category_raw",
    values_to = "expected_fraction"
  ) |>
  mutate(
    category = case_when(
      category_raw == "both_close" ~ "Both close",
      category_raw == "one_close" ~ "One close",
      category_raw == "both_distant" ~ "Both distant"
    ),
    n_expected = expected_fraction * n_alleles / 2 # divide by 2 because each cell has 2 alleles
  ) |>
  select(light, category, n_expected, expected_fraction)
```

### Chi-square test comparing observed vs expected

``` r
# Merge observed and expected
chi_data <- observed_distr |>
  inner_join(expected_distr, by = c("light", "category"))

# Run chi-square test separately for each condition
chi_results <- chi_data |>
  group_by(light) |>
  summarise(
    chisq_test = list(
      chisq.test(
        x = n_observed,
        p = n_expected / sum(n_expected),
        rescale.p = TRUE
      )
    ),
    .groups = "drop"
  )

# Extract test results
chi_results_summary <- chi_results |>
  mutate(
    statistic = map_dbl(chisq_test, "statistic"),
    df = map_dbl(chisq_test, "parameter"),
    p_value = map_dbl(chisq_test, "p.value")
  ) |>
  select(light, statistic, df, p_value)
```

Write output files

``` r
write.csv(chi_results_summary, "output/chi_results_summary.csv", row.names = FALSE)
write.csv(observed_distr, "output/observed_distr_summary.csv", row.names = FALSE)
write.csv(expected_distr, "output/expected_distr_summary.csv", row.names = FALSE)
```

Session info

``` r
sessionInfo()
```

    R version 4.5.1 (2025-06-13)
    Platform: aarch64-apple-darwin20
    Running under: macOS Sequoia 15.7.1

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: America/New_York
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] SpatialTools_1.0.5 data.table_1.17.8  fs_1.6.6           lubridate_1.9.4   
     [5] forcats_1.0.1      stringr_1.6.0      dplyr_1.1.4        purrr_1.2.0       
     [9] readr_2.1.5        tidyr_1.3.1        tibble_3.3.0       ggplot2_4.0.0     
    [13] tidyverse_2.0.0   

    loaded via a namespace (and not attached):
     [1] magic_1.6-1        generics_0.1.4     stringi_1.8.7      lattice_0.22-7    
     [5] hms_1.1.4          digest_0.6.37      magrittr_2.0.4     evaluate_1.0.5    
     [9] grid_4.5.1         timechange_0.3.0   RColorBrewer_1.1-3 fastmap_1.2.0     
    [13] plyr_1.8.9         jsonlite_2.0.0     Matrix_1.7-3       Formula_1.2-5     
    [17] scales_1.4.0       spBayes_0.4-8      abind_1.4-8        cli_3.6.5         
    [21] rlang_1.1.6        withr_3.0.2        yaml_2.3.10        tools_4.5.1       
    [25] reshape2_1.4.4     tzdb_0.5.0         coda_0.19-4.1      vctrs_0.6.5       
    [29] R6_2.6.1           lifecycle_1.0.4    pkgconfig_2.0.3    pillar_1.11.1     
    [33] gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0         xfun_0.54         
    [37] tidyselect_1.2.1   rstudioapi_0.17.1  knitr_1.50         farver_2.1.2      
    [41] htmltools_0.5.8.1  rmarkdown_2.30     compiler_4.5.1     S7_0.2.0          
    [45] sp_2.2-0          
