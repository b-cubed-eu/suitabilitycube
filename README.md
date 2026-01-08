# Suitability Cube 
The **Suitability Cube (SC)** is a reproducible, multidimensional structure designed to critically evaluate the outcomes of species distribution models (SDMs) in relation to the environmental data on which they are built. By organizing information across space, time, and taxa, the SC links species occurrence records with environmental predictors within a unified analytical framework, allowing a more transparent assessment of habitat suitability patterns and model reliability.

SC incorporates measures of environmental distance, which quantify how different the environmental conditions at prediction sites are from those represented in the training data. This allows identifying where model projections are environmentally supported and where they involve extrapolation into unfamiliar conditions. It also integrates metrics of niche breadth (derived from hypervolumes) to describe the ecological space occupied by each species and the Area of Applicability (AOA) to delineate the environmental domain within which model predictions can be considered valid.

By integrating these components, the SC provides a coherent system for comparing species and time periods, highlighting environmentally uncertain regions, and identifying where models can be meaningfully interpreted. Ultimately, the Suitability Cube enhances the critical evaluation of SDMs by linking what the model predicts to where those predictions should be trusted, fostering more transparent, reproducible, and ecologically grounded modelling practices. 

## Key concepts 
* **Hypervolume (HV)**: it represents a species’ ecological niche as an n-dimensional region defined by environmental variables, following the concept introduced by Hutchinson (1957). Each axis corresponds to an independent and ecologically relevant factor, such as temperature or precipitation, and the resulting hypervolume encompasses all combinations of conditions under which the species can persist.
* **Environmental Distance (Dissimilarity Index, DI)**: It quantifies how different the environmental conditions at a prediction site are from those represented in the training (occurrence) data. Following the method of Meyer & Pebesma (2021), DI is computed as a standardized distance in a multidimensional predictor space that has been (i) standardized, (ii) weighted by variable importance, and (iii) expressed relative to typical distances between training samples. This provides a unitless, comparable measure of environmental novelty
* **Area of Applicability (AOA)**: the Area of Applicability (AOA) defines the spatial domain where model predictions can be considered environmentally supported by the training data. It is computed directly from the Dissimilarity Index (DI) and provides a binary, spatially explicit assessment of whether a prediction site lies inside or outside the environmental space represented during model training (Meyer & Pebesma 2021).

## Conceptual workflow 
The conceptual workflow underlying the Suitability Cube (SC) consists of 4 main phases: **data acquisition**, **pre-processing**, **indicator computation**, and **cube building**. The workflow is organized into sequential, reproducible steps that ensure a coherent progression from data acquisition to cube construction while remaining flexible enough to accommodate different user needs and modelling contexts.
1. **Data download**. In the first stage, the necessary input data are collected and harmonized: the **national boundary**, **bioclimatic predictors** and **species occurrence data**.
To reduce collinearity among predictors, a correlation-based variable selection is performed on the present dataset, and the same subset of variables is applied to future scenarios to maintain temporal consistency.
2. **Pre-processing**. Following data acquisition, environmental and occurrence datasets are cleaned and harmonized. 
3. **Indicators**. This phase focuses on deriving three diagnostic indicators (HV, DI, AOA) that form the informational core of the Suitability Cube. All of them are computed from the same environmental predictors and occurrence data used to train the SDMs, ensuring consistency between model inputs and their evaluation.
4. **Cube building**. This phase integrates all indicators into a unified, three-dimensional structure that organizes information across space, species, and time.

<p align="center">
  <img src="images/GBIF%20occurrences%20from%20rgbif.png" alt="GBIF occurrences from rgbif" width="600"><br>
  <em> Workflow</em>
</p>


The final output is a **three-dimensional data cube** (cell × species × time) that integrates all indicators derived from SDM inputs. This cube provides a reproducible and transparent framework for exploring how species occupy environmental space, identifying areas of high uncertainty or extrapolation, and evaluating the robustness of model predictions through space and time.

<p align="center">
  <img src="images/vsc_page-0014.jpg" alt="cube" width="800"><br>
  <em> Final output with the 3 dimensions </em>
</p>


## Installation

You can install the development version of suitabilitycube from GitHub with:

```r
install.packages("remotes")
remotes::install_github("b-cubed-eu/suitabilitycube")
```
Then load the package: 
```r
library(suitabilitycube)
```


## Example

### Setup
This configuration defines the spatial, temporal, and taxonomic parameters used throughout the workflow.
```r
## 1. User inputs (EDIT here)
params <- list(
  species      = c("Bufo bufo", "Bufotes viridis", "Bombina variegata"),
  country_name = "Italy",
  country_iso  = "IT",
  res_arcmin   = 2.5,
  ssp_code     = "245",
  gcm_model    = "BCC-CSM2-MR",
  period       = "2041-2060",
  outdir       = tempdir(),      # change to a persistent path if desired
  gbif_years   = c(2010, 2020),
  gbif_limit   = 20000,
  cor_thr      = 0.7,
  cor_frac     = 0.10,
  grid_cellsize_deg = 0.25,      # ~25 km
  grid_square      = FALSE       # FALSE = hex, TRUE = square
)
```
### Data download
Bioclimatic variables for the present period were downloaded from WorldClim, while those for the future were retrieved from CMIP6 for the selected Global Circulation Model (BCC-CSM2-MR), scenario (SSP245), and time window (2041–2060).  Both datasets were downloaded at a spatial resolution of 2.5′ (~5 km). 

After download, each dataset was cropped and masked to the national boundary of Italy, retrieved with ```geodata::gadm()```. To guarantee spatial comparability, the present raster stack was aligned to the future using the custom helper function ```align_to()```, which resamples the data to a common extent, resolution, and coordinate reference system. 

Following alignment, a correlation-based variable selection was applied to the present bioclimatic dataset to reduce multicollinearity among predictors. Pairwise correlations were computed on a 10% random sample of grid cells, and variables exceeding a threshold of ```r > 0.7``` were iteratively removed using the custom function ```drop_high_corr()```. 

The resulting subset of predictors, ```bio3```, ```bio4```, ```bio8```, ```bio11```, ```bio14```, and ```bio16```, was retained as the final set of environmental variables. The same subset was then applied to the future dataset to ensure temporal consistency.

Species occurrence records were retrieved from GBIF using the R package ```rgbif```:
* Three amphibians: Bufo bufo, Bufotes viridis, Bombina variegata) 
* Country = IT, 
* Years = 2010–2020, 
* Record cap = 20,000 per species.



```r
# country boundary
country_vec <- geodata::gadm(params$country_name, level = 0, path = params$outdir)

# present bioclimatic predictors
bio_present <- geodata::worldclim_country(
  country = params$country_name, var = "bio", res = params$res_arcmin, path = params$outdir
) |> terra::crop(country_vec) |> terra::mask(country_vec)

# future bioclimatic predictors
bio_future <- geodata::cmip6_world(
  model = params$gcm_model, ssp = params$ssp_code, time = params$period,
  var = "bio", res = params$res_arcmin, path = params$outdir
) |> terra::crop(country_vec) |> terra::mask(country_vec)

# align_to aligns a raster stack to a target grid (bilinear)
bio_present_aligned <- align_to(bio_present, bio_future)

## 2. Variable selection (on PRESENT)
cor_res   <- drop_high_corr(bio_present_aligned, thr = params$cor_thr, frac = params$cor_frac)
cmat      <- cor_res$cor
vars_keep <- cor_res$selected
vars_drop <- cor_res$dropped

# align names (do this once right after creating bio_future)
names(bio_future) <- names(bio_present_aligned)

# keep same variables in the future
bio_present_sel <- bio_present_aligned[[vars_keep]]
bio_future_sel  <- bio_future[[vars_keep]]

plot(bio_present_sel)

## 3. GBIF occurrences --------------------------------------------------------
occ_list <- gbif_occ_list(params$species, params$country_iso, params$gbif_years, params$gbif_limit)

```

### Hypervolume
In this workflow, the hypervolume is empirically estimated from observed species occurrences and their associated environmental predictors. The computation is performed using the R package ```hypervolume``` (Blonder et al., 2018), which applies a Gaussian kernel density estimation (KDE) method (```hypervolume_gaussian```) to model the probability density of species occurrences in environmental space. The hypervolume is computed only for the present period, as it depends on empirical occurrences that are not available for future conditions and could be affected by dispersal limitations or niche shifts.

* ```extract_predictors_at_points```: extracts raster values at point locations, keeps only selected predictor variables, removes rows with missing values, and drops numeric predictors with zero variance.
* ```z_transform```: numeric columns are standardized (mean 0, sd 1). Missing values are replaced with the column mean prior to standardization. Columns with zero variance become all zeros
* ```compute_global_bw```: uses ```hypervolume::estimate_bandwidth``` on the predictor matrix. Returns 1 if estimation fails or produces non-finite values.
* ```hyp_calc```: computes a Gaussian hypervolume using ```hypervolume::hypervolume_gaussian``` and returns the volume. Returns NA if there are insufficient unique points or if the computation fails.
* Wrapper loop: iterates over all species to compute individual hypervolume values for the present period

``` r
pred_vars_present <- names(bio_present_sel)
hv_by_species     <- setNames(vector("list", length(params$species)), params$species)

# hypervolume calculation
for (sp in params$species) {
  occ_sf <- occ_list[[sp]]
  if (is.null(occ_sf) || nrow(occ_sf) == 0) { hv_by_species[[sp]] <- NA_real_; next }
  train_df  <- extract_predictors_at_points(bio_present_sel, occ_sf, pred_vars_present)
  if (nrow(train_df) < (ncol(train_df)+1L)) { hv_by_species[[sp]] <- NA_real_; next }
  z_df <- z_transform(train_df)
  hv_by_species[[sp]] <- hyp_calc(z_df, compute_global_bw(z_df))
}

# output 
hv_by_species

# $`Bufo bufo`
# [1] 1728.868

# $`Bufotes viridis`
# [1] 1482.359

# $`Bombina variegata`
# [1] 1881.454
```

### AOA and DI
A set of functions was implemented to automate the computation of DI and the AOA for each species. Some functions overlap with those used for hypervolume estimation, ensuring consistency across indicators.
* ```extract_predictors_at_points```: extracts environmental variable values from raster layers at species occurrence locations (also used in Hypervolume computation)
* ```z_transform```: standardizes environmental predictors to make them comparable across variables (shared with Hypervolume workflow)
* ```compute_aoa_pair```: computes the 	DI for both present and future environmental rasters using ```CAST::aoa()```, which simultaneously derives the DI and the corresponding AOA.
* Wrapper loop: iterates over all species to compute per-species DI and AOA values for the present and future periods, skipping species with insufficient training data.

```r
aoa_di_by_species <- setNames(vector("list", length(params$species)), params$species)

for (sp in params$species) {
  occ_sf <- occ_list[[sp]]
  if (is.null(occ_sf) || nrow(occ_sf) == 0) { aoa_di_by_species[[sp]] <- NULL; next }
  train_df <- extract_predictors_at_points(bio_present_sel, occ_sf, pred_vars_present)
  if (nrow(train_df) < 5 || ncol(train_df) < 1) { aoa_di_by_species[[sp]] <- NULL; next }
  aoa_di_by_species[[sp]] <- compute_aoa_pair(train_df, new_present = bio_present_sel, new_future = bio_future_sel)
}

# output of one species and one scenario
aoa_di_by_species$`Bufo bufo`$present

# DI:
# class       : SpatRaster 
# dimensions  : 278, 285, 1  (nrow, ncol, nlyr)
# resolution  : 0.04166667, 0.04166667  (x, y)
# extent      : 6.625, 18.5, 35.5, 47.08333  (xmin, xmax, ymin, ymax)
# coord. ref. : lon/lat WGS 84 (EPSG:4326) 
# source(s)   : memory
# varname     : wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2041-2060 
# name        :       DI 
# min value   : 0.000000 
# max value   : 2.289293 

# AOA:
# class       : SpatRaster 
# dimensions  : 278, 285, 1  (nrow, ncol, nlyr)
# resolution  : 0.04166667, 0.04166667  (x, y)
# extent      : 6.625, 18.5, 35.5, 47.08333  (xmin, xmax, ymin, ymax)
# coord. ref. : lon/lat WGS 84 (EPSG:4326) 
# source(s)   : memory
# varname     : wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2041-2060 
# name        : AOA 
# min value   :   0 
# max value   :   1 


# Predictor Weights:
#  wc2.1_30s_bio_3 wc2.1_30s_bio_4 wc2.1_30s_bio_8 wc2.1_30s_bio_11 wc2.1_30s_bio_14 wc2.1_30s_bio_16
# 1               1               1               1                1                1                1

# AOA Threshold: 0.1816886
```

### Cube building
This phase organizes all previously computed indicators, HV, DI and AOA, into a unified, multidimensional structure that enables consistent spatial, temporal, and taxonomic analysis. The resulting cube is implemented in R using the ```stars``` package, which supports multidimensional data aligned by space, time, and species. Custom functions developed for this workflow are: 
* ```as_stars_on_grid```: aggregates raster values to polygon grid cells and converts outputs to stars format
* ```build_metric_cube```: builds multi-dimensional stars cubes (AOA or DI) by aggregating indicators across species and time
* ```build_hv_cube```: generates the hypervolume cube by inserting scalar niche-size values per species
* ```merge_cubes```: combines all indicator cubes into a single multi-attribute data cube aligned by dimensions





