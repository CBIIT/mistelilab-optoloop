# DNA FISH spot distance & RNA FISH calculations integrated with GFP
intensity analysis
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
# Provide experimental conditions in desired order:
light_levs <- c("Ctrl", "gRNA")

plate_layout <- tibble(
  row = c(8L, 9L),
  column = c(13L, 13L),
  well_index = c("H13", "I13"),
  light = c("Ctrl", "gRNA")
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
    area, perimeter, solidity,
    GFP_int = `ch5_mean_intensity`,
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
    x = x_location, y = y_location, integrated_intensity
  ) |>
  # Convert xy distances from pixels to microns
  mutate(across(x:y, list(mic = ~ .x * XY_RES))) |>
  # filter spots not in the nucleus
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
  filter(area >= 50, solidity >= 0.92)
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
spot_filt_criterion <- quote(n_spots_channel_4 >= 2 & n_spots_channel_2 >= n_spots_channel_4) # set criterion for filtering

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

### Calculate distances between all possible DNA spot pairs for each cell

``` r
gf_dist_2D <- spot_filt_tbl[channel %in% 2:4, # Only far-Red and Green spots
  reshape2::melt(
    SpatialTools::dist2(
      as.matrix(.SD[
        channel == 2,
        list(x_mic, y_mic)
      ]),
      as.matrix(.SD[
        channel == 4,
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

### Find minimum DNA spot distance Minimum distance is calculated based on channel 4 spots as reference

``` r
setkey(
  gf_dist_2D, row, column, well_index,
  field_index, time_point, cell_index, f_index
)

gf_dist_min_2D <- gf_dist_2D[, .SD[which.min(gf_dist), ],
  by = key(gf_dist_2D)
] |>
  left_join(select(cell_tbl, column:well_index, light, field_index:GFP_int),
    by = c(
      "column", "row",
      "well_index", "light", "field_index",
      "time_point", "cell_index"
    )
  ) |>
  mutate(gf_dist_norm = gf_dist^2 / area) |>
  filter(gf_dist <= 1.5)
```

### Per well summary calculations of DNA distances

Calculate fraction of close alleles

``` r
overlap_dist <- 0.27 # Threshold distance (in um) for contacts

overlaps <- gf_dist_min_2D |>
  group_by(
    well_index, light,
    row, column
  ) |>
  summarize(close_alleles = mean(gf_dist < overlap_dist))
```

    `summarise()` has grouped output by 'well_index', 'light', 'row'. You can
    override using the `.groups` argument.

``` r
overlaps$light <- factor(overlaps$light, levels = c(light_levs))
```

Summarize the number of DNA spots per cell data into well-level data.

``` r
well_tbl_n_spots <- cell_n_spots_filt_tbl |>
  group_by(
    well_index, light,
    row, column
  ) |>
  summarise(N_cells = sum(!!spot_filt_criterion)) |>
  ungroup() |>
  select(
    row, column,
    well_index, light, N_cells
  ) |>
  arrange(
    well_index,
    row,
    column
  )
```

    `summarise()` has grouped output by 'well_index', 'light', 'row'. You can
    override using the `.groups` argument.

Summarize the DNA spot distance level data into well-level data.

``` r
well_tbl_gf_dist <- gf_dist_min_2D |>
  group_by(
    well_index, light,
    row, column
  ) |>
  summarise(across(
    c(gf_dist),
    list(
      mean = ~ mean(.x, rm.na = TRUE),
      median = ~ median(.x, rm.na = T)
    )
  )) |>
  ungroup() |>
  select(
    row, column, light,
    well_index, gf_dist_mean:gf_dist_median
  ) |>
  arrange(
    well_index,
    row,
    column
  )
```

    `summarise()` has grouped output by 'well_index', 'light', 'row'. You can
    override using the `.groups` argument.

Bring all the DNA FISH well-level data information together.

``` r
well_tbl <- well_tbl_n_spots |>
  left_join(well_tbl_gf_dist, by = c("row", "column", "well_index", "light")) |>
  left_join(overlaps, by = c("row", "column", "well_index", "light")) |>
  rename(
    Condition = light,
    mean_dist = gf_dist_mean,
    median_dist = gf_dist_median
  )
```

Bin cells according to GFP intensity

``` r
gf_dist_binned <- gf_dist_min_2D |>
  mutate(GFP_bin = case_when(
    GFP_int < 100 ~ "<100",
    GFP_int >= 100 & GFP_int < 300 ~ "100–299",
    GFP_int >= 300 & GFP_int < 700 ~ "300–699",
    GFP_int >= 700 & GFP_int < 2000 ~ "700–1999",
    GFP_int >= 2000 ~ "2000+"
  )) |>
  mutate(GFP_bin = factor(GFP_bin, levels = c("<100", "100–299", "300–699", "700–1999", "2000+")))
```

Summarize fraction of overlapping alleles per GFP bin

``` r
overlaps_by_bin <- gf_dist_binned |>
  group_by(GFP_bin, light) |>
  summarise(
    n_alleles = n(),
    close_alleles = mean(gf_dist < overlap_dist, na.rm = TRUE),
    .groups = "drop"
  ) |>
  rename(Condition = light)
```

### Compute per-cell DNA FISH metrics

``` r
# Calculate per-cell mean DNA distance and fraction of close alleles
gf_per_cell_metrics <- gf_dist_min_2D |>
  group_by(well_index, light, column, row, field_index, time_point, cell_index, GFP_int) |>
  summarise(
    mean_gf_dist = mean(gf_dist, na.rm = TRUE),
    frac_close = mean(gf_dist < overlap_dist, na.rm = TRUE),
    .groups = "drop"
  )
```

### RNA FISH data analysis

Define RNA FISH channel

``` r
# Define RNA FISH channel
channel_rna <- 3
```

Count RNA FISH spots per cell

``` r
rna_spots_tbl <- spots_tbl |>
  filter(channel == channel_rna) |>
  semi_join(cell_n_spots_filt_tbl,
    by = c(
      "column", "row", "field_index", "light",
      "time_point", "cell_index"
    )
  )

rna_spots_per_cell <- rna_spots_tbl |>
  group_by(well_index, light, column, row, field_index, time_point, cell_index) |>
  summarise(n_rna_spots = n(), .groups = "drop")
```

Merge per-cell DNA and RNA FISH data

``` r
rna_dna_tbl <- cell_n_spots_filt_tbl |>
  select(well_index, light, column, row, field_index, time_point, cell_index) |>
  left_join(rna_spots_per_cell,
    by = c("well_index", "light", "column", "row", "field_index", "time_point", "cell_index")
  ) |>
  left_join(gf_per_cell_metrics,
    by = c("well_index", "light", "column", "row", "field_index", "time_point", "cell_index")
  ) |>
  filter(!is.na(GFP_int)) |>
  mutate(n_rna_spots = replace_na(n_rna_spots, 0))
```

Bin cells with RNA data according to GFP intensity

``` r
rna_dna_tbl_GFPbinned <- rna_dna_tbl |>
  mutate(GFP_bin = case_when(
    GFP_int < 100 ~ "<100",
    GFP_int >= 100 & GFP_int < 300 ~ "100–299",
    GFP_int >= 300 & GFP_int < 700 ~ "300–699",
    GFP_int >= 700 & GFP_int < 2000 ~ "700–1999",
    GFP_int >= 2000 ~ "2000+"
  )) |>
  mutate(GFP_bin = factor(GFP_bin, levels = c("<100", "100–299", "300–699", "700–1999", "2000+")))
```

Summarize RNA spot number per GFP intensity bin

``` r
rna_dna_tbl_GFP_summary <- rna_dna_tbl_GFPbinned |>
  group_by(GFP_bin, light) |>
  summarise(
    n_cells = n(),
    Avg_n_RNA_spots = mean(n_rna_spots, na.rm = TRUE),
    .groups = "drop"
  ) |>
  rename(Condition = light)
```

Save CSV output files

``` r
write.csv(rna_dna_tbl_GFP_summary, "output/RNA_GFP-bin_summary.csv", row.names = FALSE)
write.csv(overlaps_by_bin, "output/Contacts_GFP-bin_summary.csv", row.names = FALSE)
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
     [1] knitr_1.50         SpatialTools_1.0.5 ggthemes_5.1.0     data.table_1.17.8 
     [5] fs_1.6.6           lubridate_1.9.4    forcats_1.0.1      stringr_1.6.0     
     [9] purrr_1.2.0        readr_2.1.5        tidyr_1.3.1        tibble_3.3.0      
    [13] ggplot2_4.0.0      tidyverse_2.0.0    dplyr_1.1.4       

    loaded via a namespace (and not attached):
     [1] magic_1.6-1        generics_0.1.4     stringi_1.8.7      lattice_0.22-7    
     [5] hms_1.1.4          digest_0.6.37      magrittr_2.0.4     evaluate_1.0.5    
     [9] grid_4.5.1         timechange_0.3.0   RColorBrewer_1.1-3 fastmap_1.2.0     
    [13] plyr_1.8.9         jsonlite_2.0.0     Matrix_1.7-3       Formula_1.2-5     
    [17] scales_1.4.0       spBayes_0.4-8      abind_1.4-8        cli_3.6.5         
    [21] rlang_1.1.6        withr_3.0.2        yaml_2.3.10        tools_4.5.1       
    [25] reshape2_1.4.4     tzdb_0.5.0         coda_0.19-4.1      vctrs_0.6.5       
    [29] R6_2.6.1           lifecycle_1.0.4    pkgconfig_2.0.3    pillar_1.11.1     
    [33] gtable_0.3.6       Rcpp_1.1.0         glue_1.8.0         xfun_0.54         
    [37] tidyselect_1.2.1   rstudioapi_0.17.1  farver_2.1.2       htmltools_0.5.8.1 
    [41] rmarkdown_2.30     compiler_4.5.1     S7_0.2.0           sp_2.2-0          
