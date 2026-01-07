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

# keep same variables in the future
bio_present_sel <- bio_present_aligned[[vars_keep]]
bio_future_sel  <- bio_future[[vars_keep]]

plot(bio_present_sel)

## 3. GBIF occurrences --------------------------------------------------------
occ_list <- gbif_occ_list(params$species, params$country_iso, params$gbif_years, params$gbif_limit)

```

``` r
library(suitabilitycube)
## basic example code
```
