
# Change backgorund color grey to white, blue a bit darker
# Increase point size and surround with black supplementary material

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Ellis-Soto, Chapman, Koeltz 
# Biodiversity data shaped by socioeconomics risk biasing upcoming U.S. biodiversity assessments
# Figure 1
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Test cases for the Social, Political Consequences Manuscript
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
require(philentropy)
require(coin)
library(tidycensus)
library(tidyverse)
library(dismo)
library(dplyr)
library(sf)
require(lubridate)
require(rgbif)
require(stringi)
library(sf)

library(spData)
require(leaflet)
library(wesanderson)
require(gridExtra)
require(patchwork)
library(dplyr)
require(climateR)
library(tidyr)
library(alphahull)
require(LaplacesDemon)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

options(tigris_use_cache = TRUE)

cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)}
st_par <- function(sf_df, sf_func, n_cores, ...){
  # Create a vector to split the data set up by.
  split_vector <- rep(1:n_cores, each = nrow(sf_df) / n_cores, length.out = nrow(sf_df))
  # Perform GIS analysis
  split_results <- split(sf_df, split_vector) %>%
    parallel::mclapply(function(x) sf_func(x, ...), mc.cores = n_cores)
  # Combine results back together. Method of combining depends on the output from the function.
  if (class(split_results[[1]]) == 'list' ){
    result <- do.call("c", split_results)
    names(result) <- NULL
  } else {
    result <- do.call("rbind", split_results)
  }
  # Return result
  return(result)
}

lonlat_to_state <- function(pointsDF,
                            states = spData::us_states,
                            name_col = "NAME") {
  ## Convert points data.frame to an sf POINTS object
  pts <- st_as_sf(pointsDF, coords = 1:2, crs = 4326)
  
  ## Transform spatial data to some planar coordinate system
  ## (e.g. Web Mercator) as required for geometric operations
  states <- st_transform(states, crs = 3857)
  pts <- st_transform(pts, crs = 3857)
  
  ## Find names of state (if any) intersected by each point
  state_names <- states[[name_col]]
  ii <- as.integer(st_intersects(pts, states))
  state_names[ii]
}

# col_pal = wes_palette("Darjeeling2", 5)[1:2]
col_pal = c('#046C9A', 'white')


# Add if statements: 
# Lanternfly - Lycorma delicatula: New Jersey
# Mosquito - Aegis aegypti: Texas
# Monarch - Danaus plexippus: California

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Download the biodiversity data:
# lanternfly_usage_key = name_backbone("Lycorma delicatula")$usageKey
# monarch_usage_key = name_backbone("Danaus plexippus")$usageKey
# mosquito_usage_key = name_backbone("Aegis aegypti")$usageKey
# 
# # Occ downlaod begin:
# oc_d = occ_download(pred("taxonKey", lanternfly_usage_key),format = "SIMPLE_CSV", user = 'diego_ellis_soto', pwd = 'Atelopus1!', email = 'diego.ellissoto@yale.edu')
# 
# species_records <- occ_download_get(oc_d, overwrite = TRUE) %>%
#   occ_download_import() %>% filter(countryCode == 'US')
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Load stored biodiversity records:
# Lanternfly
load('/Users/diegoellis/projects/Proposals_funding/Yale_internal_grants/Redlining/Lanternfly_0021674-230224095556074.Rdata')
lanternfly = lantern_fly %>% dplyr::filter(stateProvince %in% c('Nj', 'New Jersey')) %>%
  mutate(states = lonlat_to_state(
    data.frame(
      decimalLongitude,
      decimalLatitude
    )
  ),
  states_abbrev = state.abb[match(states,state.name)]
  ) %>% filter(states_abbrev == 'NJ')
rm(lantern_fly)

# Monarch:
load('/Users/diegoellis/projects/Proposals_funding/Yale_internal_grants/Redlining/Danaus plexippus (Linnaeus, 1758).Rdata')
monarch = gbif_data_d %>% dplyr::filter(gbif_data_d$stateProvince %in% c('Ca', 'California')) %>%
  mutate(states = lonlat_to_state(
    data.frame(
      decimalLongitude,
      decimalLatitude
    )
  ),
  states_abbrev = state.abb[match(states,state.name)]
  )  %>% filter(states_abbrev == 'CA')
rm(gbif_data_d)

# Mosquito
load('/Users/diegoellis/projects/Proposals_funding/Yale_internal_grants/Redlining/Aedes aegypti (Linnaeus, 1762).Rdata')
mosquito = gbif_data_d  %>% dplyr::filter(gbif_data_d$stateProvince %in% c('Tx', 'Texas')) %>% 
  mutate(states = lonlat_to_state(
    data.frame(
      decimalLongitude,
      decimalLatitude
    )
  ),
  states_abbrev = state.abb[match(states,state.name)]
  )   %>% filter(states_abbrev == 'TX')
rm(gbif_data_d)


# Now make the function which inputs the 
# species = mosquito
# species = monarch
# species = lanternfly

plot_social_economics_terraclimate_species = function(species, outdir = '/Users/diegoellis/projects/Manuscripts_collab/EnvironmentalJustice/NationalNatureAssesment/PNAS_Fig/'){
  
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  # Get Income information
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  state_income <- get_acs(
    state = unique(species$states_abbrev),
    # county = "Orange",
    geography = "block group",
    variables = "B19013_001",
    geometry = TRUE,
    year = 2020
  )
  require(tigris)
  require(sf)
  
  # state_income <- get_acs(
  #   state = unique(species$states_abbrev),
  #   # county = "Orange",
  #   geography = "block group",
  #   variables = "B19013_001",
  #   geometry = TRUE,
  #   year = 2020,
  #   cb = FALSE
  # ) |>   erase_water(year = 2020)
  
  # 
  species_sf = species %>% st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = st_crs(4326))  %>% st_transform(st_crs(state_income))
  
  # St_join new jersey income with the species:
  gbif_state_USCENSUS = st_join(species_sf, state_income)
  
  # Permutation test:
  group1 = data.frame(income = round(gbif_state_USCENSUS$estimate),
                      group = as.factor('GBIF'))
  group2 = data.frame(income = round(state_income$estimate), 
                      group =as.factor('State'))
  
  income_perm = rbind(group1, group2) |> drop_na()
  
  message('Permutation income')
  
  
  print(
    oneway_test(income ~ group,
                data = income_perm)
  )
  
  
  # Median household income
  if(unique(species_sf$states_abbrev) == "NJ"){
    median_hh_income =   82545
    poverty = 35801
  }
  
  if(unique(species_sf$states_abbrev) == "TX"){
    median_hh_income =   64034
    poverty = 29950
  }
  
  if(unique(species_sf$states_abbrev) == "CA"){
    median_hh_income = 80440   
    poverty = 36900
  }
  
  poverty_line_lower_48 = 27750
  
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  # Plotting Income Map and Income Density Plot:
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  pal_estimate_state <- colorNumeric(cm.cols1(100), domain=state_income$estimate)
  
  species_state_income_map <- state_income %>%
    ggplot(aes(fill = estimate)) + 
    geom_sf(color = NA) + 
    scale_fill_viridis_c(option = "viridis", direction = -1) + # Reversed 'viridis' palette (dark = high, light = low)
    geom_sf(data = gbif_state_USCENSUS, aes(fill = estimate), size = 0.01, color = 'black', alpha = 0.7) + # Adding transparency to points
    theme_bw() + 
    ggtitle(paste0(unique(species_sf$states_abbrev), ' species records \n across income of census tracts')) +
    theme(axis.text.x = element_text(face = "bold", size = 16),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16)) +
    theme(legend.key.size = unit(2, "lines")) # Adjusting legend size
  
  ggsave(species_state_income_map, file = paste0(
  '/Users/diegoellis/projects/Manuscripts_collab/EnvironmentalJustice/NationalNatureAssesment/TREE_ScienceSociety/Resubmission/Figures/'
  , unique(species_sf$species),'.pdf') )
  # species_state_income_map <- state_income %>%
  #   ggplot(aes(fill = estimate)) + 
  #   geom_sf(color = NA) + 
  #   scale_fill_viridis_c(option = "viridis", direction = -1) + # Reversed 'viridis' palette (dark = high, light = low)
  #   geom_sf(data = gbif_state_USCENSUS, aes(fill = estimate), size = 0.01, color = 'black', alpha = 0.7) + # Adding transparency to points
  #   theme_bw() + 
  #   ggtitle(paste0(unique(species_sf$states_abbrev), ' species records \n across income of census tracts')) +
  #   theme(axis.text.x = element_text(face = "bold", size = 16),
  #         axis.title.x = element_text(face = "bold", size = 16),
  #         axis.text.y = element_text(face = "bold", size = 16),
  #         axis.title.y = element_text(face = "bold", size = 16)) +
  #   theme(legend.key.size = unit(2, "lines")) # Adjusting legend size
  # 
  # 
  
  # species_state_income_map = state_income %>%
  #   ggplot(aes(fill = estimate)) + 
  #   geom_sf(color = NA) + 
  #   scale_fill_viridis_c(option = "magma") +
  #   geom_sf(data = gbif_state_USCENSUS, aes(fill = estimate), size= 0.01, color = 'white') +
  #   theme_bw() + ggtitle(paste0(unique(species_sf$states_abbrev), ' species records \n across income of census tracts')) +
  #   theme(axis.text.x = element_text(face = "bold", size = 16),
  #         axis.title.x = element_text(face = "bold", size = 16),
  #         axis.text.y = element_text(face = "bold", size = 16),
  #         axis.title.y = element_text(face = "bold", size = 16)) +
  #   theme(legend.key.size = unit(2, "lines")) # You can adjust the value to make it larger or smaller
  # 
  # species_state_income_map <- state_income %>%
  #   ggplot(aes(fill = estimate)) + 
  #   geom_sf(color = NA) + 
  #   scale_fill_viridis_c(option = "magma", direction = -1) + # Reversed 'viridis' palette (dark = high, light = low)
  #   geom_sf(data = gbif_state_USCENSUS, aes(fill = estimate), size = 0.01, color = 'black', alpha = 0.7) + # Adding transparency to points
  #   theme_bw() + 
  #   ggtitle(paste0(unique(species_sf$states_abbrev), ' species records \n across income of census tracts')) +
  #   theme(axis.text.x = element_text(face = "bold", size = 16),
  #         axis.title.x = element_text(face = "bold", size = 16),
  #         axis.text.y = element_text(face = "bold", size = 16),
  #         axis.title.y = element_text(face = "bold", size = 16)) +
  #   theme(legend.key.size = unit(2, "lines")) # Adjusting legend size
  # 
  # species_state_income_map = state_income %>%
  #   ggplot(aes(fill = estimate)) + 
  #   geom_sf(color = NA) + 
  #   scale_fill_viridis_c(option = "magma") +
  #   geom_sf(data = gbif_state_USCENSUS, aes(fill = estimate), size= 0.01, color = 'white') +
  #   theme_bw() + ggtitle(paste0(unique(species_sf$states_abbrev), ' species records \n across income of census tracts')) +
  #   theme(axis.text.x = element_text(face = "bold", size = 16),
  #         axis.title.x = element_text(face = "bold", size = 16),
  #         axis.text.y = element_text(face = "bold", size = 16),
  #         axis.title.y = element_text(face = "bold", size = 16)) +
  #   theme(legend.key.size = unit(2, "lines")) # You can adjust the value to make it larger or smaller
  
  
  
  # species_state_income_map = state_income %>%
  #   ggplot(aes(fill = estimate)) + 
  #   geom_sf(color = NA) + 
  #   scale_fill_viridis_c(option = "viridis") +
  #   geom_sf(data = gbif_state_USCENSUS, aes(fill = estimate), size= 0.01, 
  #           color = 'black', alpha = 0.7) +
  #   theme_bw() + ggtitle(paste0(unique(species_sf$states_abbrev), ' species records \n across income of census tracts')) +
  #   theme(axis.text.x = element_text(face = "bold", size = 16),
  #         axis.title.x = element_text(face = "bold", size = 16),
  #         axis.text.y = element_text(face = "bold", size = 16),
  #         axis.title.y = element_text(face = "bold", size = 16)) +
  #   theme(legend.key.size = unit(2, "lines")) # You can adjust the value to make it larger or smaller
  
  gbif_state_USCENSUS$estimate_10k = gbif_state_USCENSUS$estimate / 10000
  
  
  # col_pal = c('#013D5D', 'white')
  col_pal = c('#046C9A', 'white')
  col_pal = c('#046C9A', 'white')
  col_pal = c('#046C9A', 'lightgoldenrod4')
  col_pal = c('#046C9A', 'beige')
  col_pal = c('#046C9A', 'beige')
  col_pal = c('#046C9A', 'bisque3')
  # Histogram of values 
  density_hist =   ggplot() +
    geom_density(aes(estimate/1000,
                     fill = "GBIF records"),
                 alpha = .2,
                 data = gbif_state_USCENSUS, linewidth = 0.8)  +
    geom_density(aes(estimate/1000, fill = "State background income"), alpha = .2, data = state_income, linewidth = 0.8) +
    ggtitle(paste0(unique(species_sf$states_abbrev), ' species records \n across income of census tracts ', unique(gbif_state_USCENSUS$species))) +
    geom_vline(xintercept = (median_hh_income/1000), linetype = 'dashed', size = .75) + # Median household income
    geom_vline(xintercept = (poverty/1000), linetype = 'dashed', size = .75) + # Poverty
    scale_fill_manual(values = col_pal) + theme_classic() + ylab('Sampling density') + xlab('Median household income in $') +
    theme(axis.text.x = element_text(face = "bold", size = 16 ,color='black'),
          axis.title.x = element_text(face = "bold", size = 16 ,color='black'),
          axis.text.y = element_text(face = "bold", size = 16 ,color='black'),
          axis.title.y = element_text(face = "bold", size = 16 ,color='black'))+xlim(0, 300) +
    theme(legend.position="none") # Remove legend
  
  ggsave(density_hist, file = paste0(
    '/Users/diegoellis/projects/Manuscripts_collab/EnvironmentalJustice/NationalNatureAssesment/TREE_ScienceSociety/Resubmission/Figures/'
    , unique(species_sf$species),'_density_hist.pdf') )
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  # Get Race information
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  racevars <- c(White = "P2_005N", 
                Black = "P2_006N", 
                Asian = "P2_008N", 
                Hispanic = "P2_002N")
  
  state_race <- get_decennial(
    geography = "block group",
    variables = racevars,
    state = unique(species$states_abbrev),
    geometry = TRUE,
    summary_var = "P2_001N",
    year = 2020
  ) 
  
  gbif_state_race = st_join(species_sf, state_race)
  
  gbif_racial_comp = gbif_state_race  %>%
    mutate(percent = 100 * (value / summary_value)) %>% 
    drop_na(variable) 
  
  state_race = state_race  %>%
    mutate(percent = 100 * (value / summary_value))
  
  gbif_racial_comp_comp_tmp = gbif_racial_comp %>% as_tibble() %>% dplyr::select(variable, percent) %>% mutate(type='gbif')
  # save(gbif_racial_comp_comp_tmp,file ='/Users/diegoellis/Desktop/mosquito_texas_race.Rdata')
  state_race_tmp = state_race %>% as_tibble() %>% dplyr::select(variable, percent)  %>% mutate(type='state')
  
  tmp = rbind(gbif_racial_comp_comp_tmp, state_race_tmp)
  
  tmp_race = tmp
  
  tmp_race_white = tmp_race[tmp_race$variable =='White',] |> drop_na()
  tmp_race_white$type = as.factor(tmp_race_white$type)
  
  message('Permutation race')
  print(
    oneway_test(percent ~ type,
                data = tmp_race_white)
  )
  
  # col_pal = c('#046C9A', 'azure3')
  col_pal = c('#046C9A', 'white')
  col_pal = c('#046C9A', 'lightgoldenrod4')
  col_pal = c('#046C9A', 'beige')
  col_pal = c('#046C9A', 'bisque3')
  # azure3", "gbif" = "#046C9A
  
  plot_racial_makeup_white= ggplot(tmp[tmp$variable %in% 'White',], aes(percent)) +
    geom_density(data = tmp[tmp$variable %in% 'White',], aes(fill = type), alpha = 0.4, linewidth = 0.8) +
    ggtitle(paste0(unique(species_sf$states_abbrev), ' species records \n across income of census tracts ', unique(species_sf$species))) +
    scale_fill_manual(values = col_pal) + theme_classic() + ylab('Density') + xlab('Percentage (in %) of census block group \n identifying as white') +
    theme(axis.text.x = element_text(face = "bold", size = 16 ,color='black'),
          axis.title.x = element_text(face = "bold", size = 16 ,color='black'),
          axis.text.y = element_text(face = "bold", size = 16 ,color='black'),
          axis.title.y = element_text(face = "bold", size = 16 ,color='black'))+
    theme(legend.position="none") # Remove legend
  
  ggsave(plot_racial_makeup_white, file = paste0(
    '/Users/diegoellis/projects/Manuscripts_collab/EnvironmentalJustice/NationalNatureAssesment/TREE_ScienceSociety/Resubmission/Figures/'
    , unique(species_sf$species),'_race_hist.pdf') )
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  count_centroid <- state_income %>% 
    # convert to projected coord system for better centroid
    st_centroid() %>%
    st_transform(4226)
  
  # mapview(sf_data)
  
  tmp = getTerraClim(
    count_centroid,
    varname = c('tmax','tmin', 'ppt'),
    startDate = '2021-01-01',
    endDate = '2021-12-31',
    verbose = FALSE,
    dryrun = FALSE
  ) 
  
  tx_tmax <- raster(terra::app(tmp$tmax, mean))
  tx_tmin <- raster(terra::app(tmp$tmin, mean))
  tx_ppt <- raster(terra::app(tmp$ppt, mean))
  
  
  count_centroid_sf = count_centroid %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                                                       lat = sf::st_coordinates(.)[,2]) %>% as.data.frame() %>% drop_na(lon,lat)
  
  
  count_centroid_sp = SpatialPointsDataFrame(count_centroid_sf,
                                             coords = count_centroid_sf[,c('lon', 'lat')],
                                             proj4string =CRS("+proj=longlat +datum=WGS84")
  )
  
  texas_tmax_2021 = raster::extract(tx_tmax, count_centroid_sp)
  texas_ppt_2021 = raster::extract(tx_ppt, count_centroid_sp)
  
  # Now for mosquito:
  
  mosquito_centroid_sf = state_income[state_income$GEOID %in% gbif_state_USCENSUS$GEOID,] %>%   st_centroid() %>%
    st_transform(4226) %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                                         lat = sf::st_coordinates(.)[,2]) %>% as.data.frame() %>% drop_na(lon,lat)
  
  mosquito_count_centroid_sp = SpatialPointsDataFrame(mosquito_centroid_sf,
                                                      coords = mosquito_centroid_sf[,c('lon', 'lat')],
                                                      proj4string =CRS("+proj=longlat +datum=WGS84")
  )
  
  mosquito_tmax = raster::extract(tx_tmax, mosquito_count_centroid_sp)
  mosquito_ppt = raster::extract(tx_ppt, mosquito_count_centroid_sp)
  
  tx = data.frame(tmax = texas_tmax_2021, precip = texas_ppt_2021, id = 'state')
  gbif = data.frame(tmax = mosquito_tmax, precip = mosquito_ppt, id = 'gbif')
  
  A_B <- rbind(tx, gbif)
  
  A_B_env = A_B
  A_B_env$precip = round(A_B_env$precip)
  A_B_env$id = factor(A_B_env$id)
  
  A_B_env = A_B_env |> drop_na()
  message('Permutation precipitation')
  
  print(
    oneway_test(precip ~ id,
                data = A_B_env)
  )
  message('Permutation temperature')
  
  A_B_env$tmax = round(A_B_env$tmax)
  print(
    oneway_test(tmax ~ id,
                data = A_B_env)
  )
  
  
  library(alphahull)
  
  # Find the convex hull of the points being plotted
  hull_state <- A_B %>% filter(id == 'state') %>% drop_na(precip, tmax) %>% 
    slice(chull(precip, tmax))
  
  hull_gbif <- A_B %>% filter(id == 'gbif') %>% drop_na(precip, tmax) %>% 
    slice(chull(precip, tmax))
  
  
  # col_pal = wes_palette("Darjeeling2", 5)[1:2]
  # 
  # 
  # col_pal = wes_palette("Darjeeling2", 5)[1:2]
  # 
  col_pal = c('#046C9A', 'white')
  col_pal = c('#046C9A', 'bisque3')
  col_pal = wes_palette("Darjeeling2", 5)[1:2]
  col_pal = c('#046C9A', 'bisque3')
  # state could be D69C4E
  
  env_space_plot = ggplot(data = A_B, aes(x = precip, y = tmax, group = id, colour = id, fill = id)) +
    geom_point(aes(alpha = 0.5, color = id)) +  # Set color according to the id variable
    geom_polygon(data = hull_state, aes(fill = "state", alpha = 0.5), color = "black") +  # Add black contour for hull_state
    geom_polygon(data = hull_gbif, aes(fill = "gbif", alpha = 0.5), color = "black") +  # Add black contour for hull_gbif
    # scale_fill_manual(values = c("state" = "azure4", "gbif" = "#ECCBAE")) +  # Assign fill colors
    scale_fill_manual(values = c("state" = "bisque3", "gbif" = "#046C9A")) +  # Assign fill colors
    scale_color_manual(values = c("state" = "bisque3", "gbif" = "#046C9A")) +  # Assign point colors
    ggtitle(paste0(unique(species_sf$states_abbrev), ' environmental space')) +
    theme(panel.background = element_rect(fill = "white"),  # Set the background color to white
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.line = element_line(color = "black"),  # Set the axis line color to black
          axis.text = element_text(size = 16, color = "black"),  # Set the axis text size and color
          axis.title = element_text(size = 16, color = "black"),  # Set the axis title size and color
          plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Set the plot margin
          aspect.ratio = 1,  # Set the aspect ratio
          legend.key.size = unit(2, "lines")) +  # Set the legend key size
    xlab('Monthly average \n precipitation in mm') +
    ylab('Monthly average \n temperature in ° C') +
    stat_ellipse(type = "norm", linetype = 2, lwd = 1.1) +  # Add black contour for the stat_ellipse elements
    # stat_ellipse(type = "t", lwd = 1.1) + # Add black contour for the stat_ellipse elements
    theme(axis.text.x = element_text(face = "bold", size = 16 ,color='black'),
          axis.title.x = element_text(face = "bold", size = 16 ,color='black'),
          axis.text.y = element_text(face = "bold", size = 16 ,color='black'),
          axis.title.y = element_text(face = "bold", size = 16 ,color='black'))+
    theme(legend.position="none") # Remove legend
  
  ggsave(env_space_plot, file = paste0(
    '/Users/diegoellis/projects/Manuscripts_collab/EnvironmentalJustice/NationalNatureAssesment/TREE_ScienceSociety/Resubmission/Figures/'
    , unique(species_sf$species),'_env_scatter.pdf') )
  # --- --- ---
  # Future
  # --- --- ---
  
  
  future_tmax_rcp_45 = getMACA(count_centroid,
                               varname     = c('tasmax'), # 	> tasmax [K] (Daily Maximum Near-Surface Air Temperature)
                               timeRes='month',
                               model     = 'CCSM4',# 2,
                               scenario  = c("rcp45"),
                               startDate = "2080-1-1", 
                               endDate   = "2080-12-31")
  
  future_tmin_rcp_45 = getMACA(count_centroid,
                               varname     = c('tasmin'), # 	> tasmax [K] (Daily Maximum Near-Surface Air Temperature)
                               timeRes='month',
                               model     = 'CCSM4',# 2,
                               scenario  = c("rcp45"),
                               startDate = "2080-1-1", 
                               endDate   = "2080-12-31")
  
  tx_tmax_future <- raster(terra::app(future_tmax_rcp_45$air_temperature, mean))
  tx_tmin_future <- raster(terra::app(future_tmin_rcp_45$air_temperature, mean))
  
  
  texas_tmax_future = raster::extract(tx_tmax, count_centroid_sp)
  texas_tmin_future = raster::extract(tx_tmin, count_centroid_sp)
  
  mosquito_tmax_future = raster::extract(tx_tmax, mosquito_count_centroid_sp)
  mosquito_tmin_future = raster::extract(tx_tmin, mosquito_count_centroid_sp)
  
  
  tx_fut = data.frame(tmax_fut = texas_tmax_future, tmin_fut = texas_tmin_future, id = 'state')
  gbif_fut = data.frame(tmax_fut = mosquito_tmax_future, tmin_fut = mosquito_tmin_future, id = 'gbif')
  
  A_B_future <- rbind(tx_fut, gbif_fut)
  
  # Find the convex hull of the points being plotted
  hull_state <- A_B_future %>% filter(id == 'state') %>% drop_na(tmin_fut, tmax_fut) %>% 
    slice(chull(tmin_fut, tmax_fut))
  
  hull_gbif <- A_B_future %>% filter(id == 'gbif') %>% drop_na(tmin_fut, tmax_fut) %>% 
    slice(chull(tmin_fut, tmax_fut))
  
  env_space_plot_future = ggplot(data = A_B_future, aes(x = tmin_fut, y = (tmax_fut - 273.15), group = id, colour = id, fill = id)) +
    geom_point(aes(alpha = 0.5, color = id)) +  # Set color according to the id variable
    geom_polygon(data = hull_state, aes(fill = "state", alpha = 0.5), color = "black") +  # Add black contour for hull_state
    geom_polygon(data = hull_gbif, aes(fill = "gbif", alpha = 0.5), color = "black") +  # Add black contour for hull_gbif
    # scale_fill_manual(values = c("state" = "azure4", "gbif" = "#ECCBAE")) +  # Assign fill colors
    scale_fill_manual(values = c("state" = "azure3", "gbif" = "#046C9A")) +  # Assign fill colors
    scale_color_manual(values = c("state" = "azure3", "gbif" = "#046C9A")) +  # Assign point colors
    ggtitle(paste0(unique(species_sf$states_abbrev), ' environmental space')) +
    theme(panel.background = element_rect(fill = "white"),  # Set the background color to white
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.line = element_line(color = "black"),  # Set the axis line color to black
          axis.text = element_text(size = 16, color = "black"),  # Set the axis text size and color
          axis.title = element_text(size = 16, color = "black"),  # Set the axis title size and color
          plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Set the plot margin
          aspect.ratio = 1,  # Set the aspect ratio
          legend.key.size = unit(2, "lines")) +  # Set the legend key size
    xlab('Future average \n minimum temperature in C') +
    ylab('Future average \n maximum temperature in ° C') +
    stat_ellipse(type = "norm", linetype = 2, lwd = 1.1) +  # Add black contour for the stat_ellipse elements
    stat_ellipse(type = "t", lwd = 1.1) + # Add black contour for the stat_ellipse elements
    theme(axis.text.x = element_text(face = "bold", size = 16 ,color='black'),
          axis.title.x = element_text(face = "bold", size = 16 ,color='black'),
          axis.text.y = element_text(face = "bold", size = 16 ,color='black'),
          axis.title.y = element_text(face = "bold", size = 16 ,color='black'))+
    theme(legend.position="none") # Remove legend
  
  
  
  # Save Figures
  
  require(patchwork)
  require(gridExtra)
  patchwork_mosquito = density_hist + plot_racial_makeup_white + env_space_plot
  grid.arrange(density_hist, plot_racial_makeup_white
               ,env_space_plot, ncol = 3)
  
  # (density_hist | plot_racial_makeup_white | env_space_plot)
  
  ggsave(
    filename = paste0('/Users/diegoellis/projects/Manuscripts_collab/EnvironmentalJustice/NationalNatureAssesment/PNAS_Fig/',unique(species$species), '_multipanel_terraclimate_v2.png'),
    patchwork_mosquito, height=18, width=18, units="in")
  
  ggsave(paste0(
    outdir,unique(species_sf$species),"_income_map_v2.jpeg"),
    species_state_income_map,
    "jpeg", height=12, width=12, units="in")
  
  ggsave(paste0(
    outdir,unique(species_sf$species),"income_density_v2.jpeg"),
    density_hist,
    "jpeg", height=12, width=12, units="in")
  
  ggsave(paste0(
    outdir,unique(species_sf$species),"racial_density_v2.jpeg"),
    plot_racial_makeup_white,
    "jpeg", height=12, width=12, units="in")
  
  ggsave(paste0(
    outdir,unique(species_sf$species),"current_temp_precip_v2.jpeg"),
    env_space_plot_future,
    "jpeg", height=12, width=12, units="in")
  
  ggsave(paste0(
    outdir,unique(species_sf$species),"future_temp_precip_v2.jpeg"),
    env_space_plot_future,
    "jpeg", height=12, width=12, units="in")
  
  # --- --- --- --- --- --- --- --- --- --- --- ---
  # Stats
  # --- --- --- --- --- --- --- --- --- --- --- ---
  
  # --- --- ---
  # Income
  # --- --- ---
  
  # Is this different:
  # chisq.test(round(gbif_state_USCENSUS$estimate_10k), round(state_income$estimate/10000), correct=FALSE)
  
  
  # Remove NA values
  income_species <- gbif_state_USCENSUS[!is.na(gbif_state_USCENSUS$estimate),]
  income_all <- state_income[!is.na(state_income$estimate),]
  
  # Compute empirical cumulative distribution function (ECDF)
  ecdf_mosquitoes <- ecdf(round(income_species$estimate/1000))
  ecdf_all <- ecdf(round(income_all$estimate/1000) )
  
  # Create a sequence of percentiles
  percentiles <- seq(0, 1, by = 0.01)
  
  # Calculate the income values corresponding to each percentile
  income_percentiles_species <- quantile(
    round( income_species$estimate / 1000), percentiles, na.rm = TRUE)
  
  income_percentiles_all <- quantile(
    round( income_all$estimate / 1000), percentiles, na.rm = TRUE)
  
  # Display the results
  income_result_df <- data.frame(
    Percentile = percentiles * 100,
    Income_Species = income_percentiles_species,
    Income_All = income_percentiles_all
  )
  
  # print(result_df)
  # They are significantly different
  KLD(income_result_df$Income_Species,income_result_df$Income_All)$mean.sum.KLD
  KLD(income_result_df$Income_Species,income_result_df$Income_All)$intrinsic.discrepancy
  # Also significantly different: Income:
  
  print(paste0(unique(species$species), ' state vs GBIF income distribution'))
  ks.test(round(gbif_state_USCENSUS$estimate), round(state_income$estimate)) # Is normal distributed
  ks.test(round(income_result_df$Income_Species), round(income_result_df$Income_All)) # Is normal distributed
  
  print(chisq.test(income_result_df$Income_Species,income_result_df$Income_All))
  wilcox.test(round(income_result_df$Income_Species), round(income_result_df$Income_All))
  wilcox.test(round(gbif_state_USCENSUS$estimate), round(state_income$estimate))
  t.test(round(gbif_state_USCENSUS$estimate), round(state_income$estimate))
  
  
  print('Environment')
  # --- --- --- --- --- --- --- --- --- --- --- ---  
  # Environment current temperature:
  # --- --- --- --- --- --- --- --- --- --- --- ---
  ecdf_mosquitoes <- ecdf(round(gbif$tmax))
  ecdf_all <- ecdf(round(tx$tmax) )
  # Create a sequence of percentiles
  percentiles <- seq(0, 1, by = 0.01)
  # Calculate the income values corresponding to each percentile
  tmax_percentiles_species <- quantile(
    round( gbif$tmax, 1), percentiles, na.rm = TRUE)
  
  tmax_percentiles_all <- quantile(
    round( tx$tmax,1), percentiles, na.rm = TRUE)
  # Display the results
  tmax_current_result_df <- data.frame(
    Percentile = percentiles * 100,
    tmax_Species = tmax_percentiles_species,
    tmax_All = tmax_percentiles_all
  )
  # They are significantly different
  KLD(tmax_current_result_df$tmax_Species,tmax_current_result_df$tmax_All)
  ks.test(tmax_current_result_df$tmax_Species,tmax_current_result_df$tmax_All) # Normal
  t.test(tmax_current_result_df$tmax_Species,tmax_current_result_df$tmax_All) # Significantly difference
  print(chisq.test(tmax_current_result_df$tmax_Species,tmax_current_result_df$tmax_All))
  wilcox.test((tmax_current_result_df$tmax_Species), (tmax_current_result_df$tmax_All))
  
  
  # --- --- --- --- --- --- --- --- --- --- --- ---  
  # Environment current precipitation:
  # --- --- --- --- --- --- --- --- --- --- --- ---
  percentiles <- seq(0, 1, by = 0.01)
  # Calculate the income values corresponding to each percentile
  precip_percentiles_species <- quantile(
    round( gbif$precip), percentiles, na.rm = TRUE)
  
  precip_percentiles_all <- quantile(
    round( tx$precip), percentiles, na.rm = TRUE)
  # Display the results
  precip_result_df <- data.frame(
    Percentile = percentiles * 100,
    precip_Species = precip_percentiles_species,
    precip_All = precip_percentiles_all
  )
  # They are significantly different
  ks.test(precip_result_df$precip_Species,precip_result_df$precip_All) # Normal
  KLD(precip_result_df$precip_Species,precip_result_df$precip_All)
  t.test(precip_result_df$precip_Species,precip_result_df$precip_All)
  print(chisq.test(precip_result_df$precip_Species,precip_result_df$precip_All))
  wilcox.test((precip_result_df$precip_Species), (precip_result_df$precip_All))
  
  # --- --- --- --- --- --- --- --- --- --- --- ---  
  # Environment current temperature:
  # --- --- --- --- --- --- --- --- --- --- --- ---
  
  # --- --- --- --- --- --- --- --- --- --- --- ---
  # Environment future
  # --- --- --- --- --- --- --- --- --- --- --- ---
  
  percentiles <- seq(0, 1, by = 0.01)
  # Calculate the income values corresponding to each percentile
  tmax_future_percentiles_species <- quantile(
    round( gbif_fut$tmax_fut,1), percentiles, na.rm = TRUE)
  
  tmax_future_percentiles_all <- quantile(
    round( tx_fut$tmax_fut,1), percentiles, na.rm = TRUE)
  # Display the results
  tmax_future_result_df <- data.frame(
    Percentile = percentiles * 100,
    tmax_future_Species = tmax_future_percentiles_species,
    tmax_future_All = tmax_future_percentiles_all
  )
  # They are significantly different
  ks.test(tmax_future_result_df$tmax_future_Species,tmax_future_result_df$tmax_future_All)
  KLD(tmax_future_result_df$tmax_future_Species,tmax_future_result_df$tmax_future_All)$intrinsic.discrepancy
  t.test(tmax_future_result_df$tmax_future_Species,tmax_future_result_df$tmax_future_All)
  print(chisq.test(tmax_future_result_df$tmax_future_Species,tmax_future_result_df$tmax_future_All))
  wilcox.test((tmax_future_result_df$tmax_future_Species), (tmax_future_result_df$tmax_future_All))
  # Now look bellow poverty line and mean household income:
  
  print(paste0(unique(income_species$species), ' income information '))
  summary(income_species$estimate)
  
  # Calculate the percentage of values below the poverty line for gbif
  below_poverty_line <- subset(income_species, estimate < poverty)
  percentage_below_poverty_line <- (nrow(below_poverty_line) / nrow(income_species)) * 100
  cat("Percentage of values below the poverty line for species  :",unique(income_species$species),' :' ,round(percentage_below_poverty_line, 2), "%\n")
  
  # Calculate the percentage of census tract below the poverty line for state
  below_poverty_line <- subset(state_income, estimate < poverty)
  percentage_below_poverty_line <- (nrow(below_poverty_line) / nrow(state_income)) * 100
  cat("Percentage of values below the poverty line for state :",unique(state_income$states_abbrev),' :' ,round(percentage_below_poverty_line, 2), "%\n")
  
  
  # --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  # Racial composition:
  
  race_tmp = rbind(gbif_racial_comp_comp_tmp, state_race_tmp)
  
  
  percentiles <- seq(0, 1, by = 0.01)
  # Calculate the income values corresponding to each percentile
  race_percentiles_species <- quantile(
    round(     race_tmp[race_tmp$type == 'gbif' & race_tmp$variable %in% 'White',]$percent
               ,1), percentiles, na.rm = TRUE)
  
  race_percentiles_all <- quantile(
    round(     race_tmp[race_tmp$type == 'state' & race_tmp$variable %in% 'White',]$percent
               ,1), percentiles, na.rm = TRUE)
  # Display the results
  race_result_df <- data.frame(
    Percentile = percentiles * 100,
    race_Species = race_percentiles_species,
    race_All = race_percentiles_all
  )
  
  ks.test((race_result_df$race_Species), (race_result_df$race_All))
  wilcox.test((race_result_df$race_Species), (race_result_df$race_All))
  t.test((race_result_df$race_Species), (race_result_df$race_All))
  
  require('ggpubr')
  
  # --- --- --- --- --- --- --- ---
  # Mean:                         #
  # --- --- --- --- --- --- --- ---
  
  
  ggqqplot(income_result_df$Income_Species, title = paste(
    'Income - '))
  ggqqplot(income_result_df$Income_All)
  
  ggqqplot(result_df$race_Species)
  ggqqplot(result_dfrace_All)
  
  ggqqplot(tmax_future_result_df$tmax_future_Species)
  ggqqplot(tmax_future_result_df$tmax_future_All)
  
  ggqqplot(precip_result_df$precip_Species)
  ggqqplot(precip_result_df$precip_All)
  
  ggqqplot(tmax_current_result_df$tmax_Species)
  ggqqplot(tmax_current_result_df$tmax_All)
  
  # Raw data income:
  shapiro.test(round(income_species$estimate/1000))
  ks.test(round(income_all$estimate/1000), 'pnorm')
  
  ggqqplot(income_species$estimate)
  ggqqplot(income_all$estimate)
  
  ggqqplot(round(income_species$estimate/1000))
  ggqqplot(round(income_all$estimate/1000))
  
  # Raw data tmax:
  
  shapiro.test(round(gbif$tmax,1))
  ks.test(round(tx$tmax,1), 'pnorm')
  
  ggqqplot(round(gbif$tmax,1))
  ggqqplot(round(tx$tmax,1))
  
  # Raw data precip:
  shapiro.test(round(gbif$precip))
  ks.test(round(tx$precip), 'pnorm')
  
  ggqqplot(round(gbif$precip))
  ggqqplot(round(tx$precip))
  
  # Raw data temperature future:
  shapiro.test(round(gbif_fut$tmax_fut))
  ks.test(round(tx_fut$tmax_fut), 'pnorm')
  
  ggqqplot(round(gbif_fut$tmax_fut,1))
  ggqqplot(round(tx_fut$tmax_fut,1))
  
  # Raw data race:
  ggqqplot(race_result_df$race_Species)
  ggqqplot(race_result_df$race_All)
  
  shapiro.test(round(race_tmp[race_tmp$type == 'gbif' & race_tmp$variable %in% 'White',]$percent))
  ks.test(round(race_tmp[race_tmp$type == 'state' & race_tmp$variable %in% 'White',]$percent), 'pnorm')
  
  
  # Literature + AvnoNet + Urbanization score + 
  # ecography functional diversity analysis =
  # density complete counts
  
  # Data dumps:
  # eBird
  # gbif
  # avonet
  
  
  
}
# bird biodiversity homogenizing or urban species or traits functional traits + OR density effect cuuz less diversity in low income environment
# avonet
# diet data or predation events
# high income probably 

plot_social_economics_terraclimate_species(mosquito)
plot_social_economics_terraclimate_species(species = lanternfly)
plot_social_economics_terraclimate_species(monarch)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure 2 of the manuscript assess unsampled niche space
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# species = lanternfly


# Remove one of the ellypses
# Add map but without background 

# --- --- --- --- --- --- --- --- --- --- --- --- ---
# Stats revision 1 ######
# --- --- --- --- --- --- --- --- --- --- --- --- ---
# Exclude the GeoIDs in which biodiversity records are present in the state:

# For income:
state_only_income = state_income[!state_income$GEOID %in% unique(gbif_state_USCENSUS$GEOID), ]
# For race:
tmp_race = tmp


# Perform a ks.test now

ks.test(round(gbif_state_USCENSUS$estimate), round(state_only_income$estimate)) # Is normal distributed



# Income:
wilcox.test(
  round(gbif_state_USCENSUS$estimate), 
  round(state_only_income$estimate)
)

wilcox.test(
  round(gbif_state_USCENSUS$estimate), 
  round(state_only_income$estimate)
)

# Perform Mann-Whitney U-test for income
wilcox.test(
  round(gbif_state_USCENSUS$estimate), 
  round(state_only_income$estimate), 
  alternative = "two.sided",
  exact = FALSE, 
  conf.int = TRUE
)
# Perform Mann-Whitney U-test for environment

head(tmp)



# Additionally perform a Mann-Whitney U-test

# 
