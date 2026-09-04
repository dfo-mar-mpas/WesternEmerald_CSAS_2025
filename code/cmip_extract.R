#load libraries -------
library(raster)
library(sf)
library(sp)
library(dplyr)
library(tidyr)
library(R.matlab)
library(stars)

sf_use_s2(FALSE)

latlong <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm_mar <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

## ---- load site polygon ----
webmr <- read_sf("data/shapefiles/webmr.shp") %>%
  st_make_valid() %>%
  st_transform(utm_mar) %>%
  mutate(area = as.numeric(st_area(.)) / 1000 / 1000) %>%   # km², accurate planar area
  st_transform(latlong) %>%                                   # back to latlong for downstream use
  dplyr::select(area)

## ---- climate projection files ----
fls <- c(list.files("c:/Users/stanleyr/Documents/GitHub/MAR_thermal_emerg/data/climate_projections/2.6/", full.names = TRUE),
         list.files("c:/Users/stanleyr/Documents/GitHub/MAR_thermal_emerg/data/climate_projections/8.5/", full.names = TRUE))
fls <- fls[!grepl("CNRM", fls)]
fls <- fls[!grepl("GFDL", fls)]

if (!dir.exists("output/climate_extracts/")) dir.create("output/climate_extracts/", recursive = TRUE)

## ---- build the webmr raster mask once ----
data <- readMat(fls[1])
names(data) <- "datout"
bdata <- brick(data$datout, xmn = -83, xmx = -41, ymn = 38, ymx = 85, crs = latlong)

# no need to pull crs back out of bdata — we already know it
webmr_sp     <- webmr %>% st_transform(latlong) %>% as_Spatial()
webmr_extent <- extent(webmr_sp)

webmr_raster_mask <- rasterize(webmr_sp, crop(bdata[[1]], webmr_extent, snap = "out"), getCover = TRUE)
webmr_raster_mask[webmr_raster_mask == 0] <- NA

save(webmr_raster_mask, file = "data/webmr_raster_mask.RData")
rm(data, bdata)

## ---- reusable extraction function for webmr ----
extract_webmr <- function(bdata, mod, climate_proj, mask, site_sf) {
  
  # mask is already in latlong — do not reproject it, and crop bdata
  # to the mask's own extent so raster::mask() always lines up
  bdata_processed <- bdata %>%
    crop(extent(mask), snap = "out") %>%
    raster::mask(mask)
  
  # guard against 1-cell rounding mismatches after crop
  if (!compareRaster(bdata_processed, mask, stopiffalse = FALSE)) {
    bdata_processed <- raster::resample(bdata_processed, mask, method = "ngb")
  }
  
  data_extract <- bdata_processed %>%
    st_as_stars() %>%
    st_as_sf() %>%
    st_transform(utm_mar) %>%
    st_intersection(site_sf %>% st_transform(utm_mar) %>% dplyr::rename(site_area = area)) %>%
    st_transform(latlong) %>%
    mutate(cell_area = as.vector(st_area(.) / 1000 / 1000)) %>%
    gather("layer", "temp", starts_with("layer.")) %>%
    mutate(month = rep(rep(1:12, each = length(layer) / 86 / 12), 86),
           year  = rep(2015:2100, each = length(layer) / 86),
           mod = mod,
           climate_proj = climate_proj) %>%
    dplyr::select(mod, climate_proj, year, month, 
                  site_area, cell_area, temp, geometry) %>%
    suppressWarnings()
  
  list(
    data  = data_extract %>% data.frame() %>% dplyr::select(-geometry),
    shape = data_extract %>% filter(month == 1, year == 2015)
  )
}

## ---- STEP 1: ensemble mean extraction ----
emmission_scenarios <- c("2.6", "8.5")

# extent to crop to, taken from the mask, with a small buffer so edge
# cells aren't clipped before the raster::mask() step inside extract_webmr()
webmr_crop_ext <- extend(extent(webmr_raster_mask), 0.1)

for (i in emmission_scenarios) {
  
  cmip_files <- fls[grep(i, fls)]
  message(paste0("Loading data for RCP ", i))
  
  # read in each model, crop immediately to webmr's extent before any
  # further processing so averaging only happens on the small subset
  for (j in 1:3) {
    temp <- readMat(cmip_files[j])
    names(temp) <- "datout"
    temp <- brick(temp$datout, xmn = -83, xmx = -41, ymn = 38, ymx = 85, crs = latlong)
    temp <- crop(temp, webmr_crop_ext, snap = "out")
    assign(paste0("emmisiondata_", j), temp)
    rm(temp)
  }
  
  message("Conducting model averaging to generate an ensemble.")
  ensemble.list <- list()
  for (s in 1:dim(emmisiondata_1)[3]) {
    temp <- stack(emmisiondata_1[[s]], emmisiondata_2[[s]], emmisiondata_3[[s]])
    ensemble.list[[s]] <- calc(temp, fun = mean, na.rm = TRUE)
    rm(temp)
  }
  
  # let brick() inherit extent/crs from the already-cropped stack —
  # do NOT reassign the original full-domain extent here
  bdata <- ensemble.list %>% stack() %>% brick()
  
  message(paste0("Working on the climate extraction for RCP ", i))
  result <- extract_webmr(bdata, mod = "Ensemble", climate_proj = i,
                          mask = webmr_raster_mask, site_sf = webmr)
  
  cmip_extract       <- result$data
  cmip_extract_shape <- result$shape
  
  message(paste0("Saving output from Ensemble ", i))
  save(cmip_extract,       file = paste0("output/climate_extracts/climate_extracts_ensemble_", gsub("\\.", "-", i), ".RData"))
  save(cmip_extract_shape, file = paste0("output/climate_extracts/climate_extracts_shape_ensemble_", gsub("\\.", "-", i), ".RData"))
  
  rm(emmisiondata_1, emmisiondata_2, emmisiondata_3, ensemble.list, bdata,
     result, cmip_extract, cmip_extract_shape)
  gc()
}



## ---- STEP 2: per-model extraction ----
message("Running the climate-model based extractions.")

for (i in seq_along(fls)) {
  
  climate_proj <- basename(dirname(fls[i]))                 # "2.6" or "8.5", robust to path depth
  mod <- strsplit(basename(fls[i]), "_")[[1]][1]             # model name from filename
  
  message(paste0("Working on ", mod, " ", climate_proj))
  
  data <- readMat(fls[i])
  names(data) <- "datout"
  bdata <- brick(data$datout, xmn = -83, xmx = -41, ymn = 38, ymx = 85, crs = latlong)
  
  result <- extract_webmr(bdata, mod = mod, climate_proj = climate_proj,
                          mask = webmr_raster_mask, site_sf = webmr)
  
  cmip_extract       <- result$data
  cmip_extract_shape <- result$shape
  
  message(paste0("Saving output from ", mod, " ", climate_proj))
  save(cmip_extract,       file = paste0("output/climate_extracts/climate_extracts_", mod, "_", climate_proj, ".RData"))
  save(cmip_extract_shape, file = paste0("output/climate_extracts/climate_extracts_shape_", mod, "_", climate_proj, ".RData"))
  
  rm(data, bdata, result, cmip_extract, cmip_extract_shape)
  gc()
}