# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure for test cases for the Social, Political Consequences Manuscript
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# gbif_species_name = 'Lycorma delicatula'
# Load packages & functions:
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
library(tidyr)

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


col_pal = wes_palette("Darjeeling2", 5)[1:2]

gbif_species_name = "Danaus plexippus" # Monarch
gbif_species_name = 'Lycorma delicatula' # Lanternfly
gbif_species_name = 'Ailanthus altissima' # Tree of heaven
gbif_species_name = 'Agrilus planipennis' # Esmerald Ash bohrer
gbif_species_name = 'Rhinella marina'
gbif_species_name = 'Halyomorpha halys'
gbif_species_name = 'Aedes aegypti'


download_plot_socioeco_gbif_data <- function(gbif_species_name = NULL){
# outdir outdir folder
  # species name
oc_d = occ_download(pred("taxonKey", name_backbone(gbif_species_name)$usageKey),format = "SIMPLE_CSV", user = 'diego_ellis_soto', pwd = 'Atelopus1!', email = 'diego.ellissoto@yale.edu')

# Add has coordinate here ! -> and more ! 
stat <- "PREPARING"
while(stat %in% c('PREPARING', 'RUNNING')) {
  met <- occ_download_meta(oc_d)
  stat <- met$status
  Sys.sleep(5)
}


gbif_data_d <- occ_download_get(oc_d, overwrite = TRUE) %>%
  occ_download_import() %>% filter(countryCode == 'US')  %>%
  drop_na(decimalLatitude, decimalLongitude) %>% 
  mutate(states = lonlat_to_state(
    data.frame(
      decimalLongitude,
      decimalLatitude
    )
  ),
  states_abbrev = state.abb[match(states,state.name)]
  )


save(gbif_data_d, file = paste0('/Users/diegoellis/projects/Proposals_funding/Yale_internal_grants/Redlining/', unique(gbif_data_d$scientificName)[1],'.Rdata'))


# Identify the state with th emost obervations
state_of_interest = which.max(table(gbif_data_d$states_abbrev))

# tail(sort(table(gbif_data_d$states_abbrev), order='desc'),1)
gbif_data_d_backup = gbif_data_d
# Subset to only the state with the most  records
 gbif_data_d = gbif_data_d %>% dplyr::filter(states_abbrev == paste0(names(state_of_interest))) %>%
#  gbif_data_d = gbif_data_d %>% dplyr::filter(states_abbrev == 'CA') %>%
  drop_na(decimalLatitude, decimalLongitude) %>%
  dplyr::select(species, genus, scientificName, countryCode, decimalLongitude, decimalLatitude, eventDate, day, month, year, institutionCode, collectionCode, basisOfRecord, dateIdentified,recordedBy, states_abbrev)

# 
state_income <- get_acs(
  state = paste0(names(state_of_interest)),
 # state = 'CA',
  # county = "Orange",
  geography = "block group",
  variables = "B19013_001",
  geometry = TRUE,
  year = 2020
)

state_income %>%
  ggplot(aes(fill = estimate)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c(option = "magma") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16)) +
  theme(legend.key.size = unit(2, "lines")) # You can adjust the value to make it larger or smaller


gbif_data_d_sf = gbif_data_d %>% st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = st_crs(4326))  %>% st_transform(st_crs(state_income))
# St intersection with records + plot biodiversity records on top
# save(gbif_data_d_sf,file = '/Users/diegoellis/Desktop/mosquito_texas_income.Rdata')


# St_join new jersey income with the species:
gbif_state_USCENSUS <- st_par(  gbif_data_d_sf
                                        , state_income 
                                        , n_cores = 4
                                        , sf_func = st_join
                                        , join = st_within)

gbif_state_USCENSUS = st_join(gbif_data_d_sf, state_income)

pal_estimate_state <- colorNumeric(cm.cols1(100), domain=state_income$estimate)

# save(state_income,file = '/Users/diegoellis/Desktop/texas_income.Rdata')

species_state_income_map = state_income %>%
  ggplot(aes(fill = estimate)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c(option = "magma") +
  geom_sf(data = gbif_state_USCENSUS, aes(fill = estimate), size= 0.01, color = 'white') +
  theme_bw() + ggtitle(paste0(unique(gbif_data_d_sf$states_abbrev), ' species records \n across income of census tracts'))

# 
# Histogram of values 
density_hist =   ggplot() +
  geom_density(aes(estimate,
                   fill = "GBIF records"),
               alpha = .2,
               data = gbif_state_USCENSUS)  +
  geom_density(aes(estimate, fill = "State background income"), alpha = .2, data = state_income) +
  # scale_fill_manual(
  #   name = "datasets",
  #   values = c(data1 = "Lanternfly collections", data2 = "New Jersey")
  #   ) +
  # ggtitle('New Jersey - Lanternfly') +
  ggtitle(paste0(unique(gbif_data_d_sf$states_abbrev), ' species records \n across income of census tracts ', unique(gbif_state_USCENSUS$species))) +
# TO DO ADD THE MEDIAN HOUSEHOLD INCOME FOR KANSAS and for each state ! 
    # geom_vline(xintercept = 89703, linetype = 'dashed') + # Median household income
  # geom_vline(xintercept = 50004, linetype = 'dashed') + # Average salary
  # geom_vline(xintercept = 35801, linetype = 'dashed') + # Poverty
  # annotate("text", x=102000, y=0.0000093, label="Household \n income", angle=20, size=2, color="grey49") +
  # annotate("text", x=53000, y=0.0000093, label="Average \n salary", angle=20, size=2, color="grey49") +
  # annotate("text", x=34000, y=0.0000093, label="Poverty", angle=20, size=2, color="grey49") +
  scale_fill_manual(values = col_pal) + theme_classic() + ylab('Frequency') + xlab('Median household income in $')

ggsave(paste0("/Users/diegoellis/Desktop/income_", unique(gbif_state_USCENSUS$species) ,".jpeg"), (species_state_income_map / density_hist), "jpeg", height=12, width=12, units="in")

# Llegue hasta Race Var: 


racevars <- c(White = "P2_005N", 
              Black = "P2_006N", 
              Asian = "P2_008N", 
              Hispanic = "P2_002N")

state_race <- get_decennial(
  geography = "block group",
  variables = racevars,
    state = paste0(names(state_of_interest)),
  #  state = 'CA',
  # county = "Harris County",
  geometry = TRUE,
  summary_var = "P2_001N",
  year = 2020
) 

# St_join new jersey income with lanternfly
gbif_state_race <- st_par(  gbif_data_d_sf
                                    , state_race 
                                    , n_cores = 4
                                    , sf_func = st_join
                                    , join = st_within)


gbif_racial_comp = gbif_state_race  %>%
  mutate(percent = 100 * (value / summary_value)) %>% 
  drop_na(variable) 

state_race = state_race  %>%
  mutate(percent = 100 * (value / summary_value))

gbif_racial_comp_comp_tmp = gbif_racial_comp %>% as_tibble() %>% dplyr::select(variable, percent) %>% mutate(type='gbif')
# save(gbif_racial_comp_comp_tmp,file ='/Users/diegoellis/Desktop/mosquito_texas_race.Rdata')
state_race_tmp = state_race %>% as_tibble() %>% dplyr::select(variable, percent)  %>% mutate(type='state')

tmp = rbind(gbif_racial_comp_comp_tmp, state_race_tmp)



plot_racial_makeup= ggplot(tmp, aes(percent)) +
  geom_density(data = tmp, aes(fill = type), alpha = 0.4) +
  facet_wrap(~ variable, scales = "free") +
  # ggtitle('New Jersey- Lanternfly') +
  ggtitle(paste0(unique(gbif_data_d_sf$states_abbrev), ' species records \n across income of census tracts ', unique(gbif_state_USCENSUS$species))) +
  scale_fill_manual(values = col_pal) + theme_classic() + ylab('Frequency') + xlab('Percentage in %')



plot_racial_makeup_white= ggplot(tmp[tmp$variable %in% 'White',], aes(percent)) +
  geom_density(data = tmp[tmp$variable %in% 'White',], aes(fill = type), alpha = 0.4, linewidth = 0.8) +
  ggtitle(paste0(unique(gbif_data_d_sf$states_abbrev), ' species records \n across income of census tracts ', unique(gbif_state_USCENSUS$species))) +
  scale_fill_manual(values = col_pal) + theme_classic() + ylab('Frequency') + xlab('Percentage (in %) of census block group \n identifying as white') +
  theme(axis.text.x = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16))

plot_racial_makeup_white

# save(state_race,file = '/Users/diegoellis/projects/texas_race.Rdata')

plot_race = state_race %>%
  mutate(percent = 100 * (value / summary_value)) %>%
  ggplot(aes(fill = percent)) +
  facet_wrap(~variable) +
  geom_sf(color = NA) +
  geom_sf(data = gbif_racial_comp, aes(fill = percent), size= 0.01, col = 'white') +
  # geom_sf(data = gbif_state_USCENSUS, aes(fill = estimate), size= 0.01, color = 'white') +
  theme_void() + 
  scale_fill_viridis_c() + 
  labs(fill = "% of population\n(2020 Census)")


plot_race = state_race %>%
  mutate(percent = 100 * (value / summary_value)) %>%
  ggplot(aes(fill = percent)) +
  facet_wrap(~variable) +
  geom_sf(color = NA) +
  geom_sf(data = gbif_racial_comp, aes(fill = percent), size= 0.01, col = 'white') +
  # geom_sf(data = gbif_state_USCENSUS, aes(fill = estimate), size= 0.01, color = 'white') +
  theme_void() + 
  scale_fill_viridis_c() + 
  labs(fill = "% of population\n(2020 Census)")




ggsave(paste0("/Users/diegoellis/Desktop/race_",unique(gbif_state_USCENSUS$species) ,"_density.jpeg"), plot_racial_makeup, "jpeg", height=12, width=12, units="in")


ggsave(paste0("/Users/diegoellis/Desktop/race_",unique(gbif_state_USCENSUS$species),".jpeg"), (plot_race / plot_racial_makeup), "jpeg", height=12, width=12, units="in")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Add education and mean age
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

state_edu_age <- get_acs(
  state = paste0(names(state_of_interest)),
  geography = 'block group',
  variables = c(
    median_age = "B01002_001",
    edu = 'B15002_001',
    bachelor_degree = 'B06009_005'    
  ),
  geometry = TRUE,
  year = 2020
)


gbif_edu_age <- st_par(  gbif_data_d_sf
                         , state_edu_age 
                         , n_cores = 4
                         , sf_func = st_join
                         , join = st_within)

density_hist_edu =   ggplot() +
  geom_density(aes(estimate,
                   fill = "GBIF records"),
               alpha = .2,
               data = gbif_edu_age[gbif_edu_age$variable =='edu',])  +
  geom_density(aes(estimate, fill = paste0("State ", unique(gbif_data_d_sf$states_abbrev))), alpha = .2, data = state_edu_age) +
  # ggtitle('Education New Jersey - Lanternfly') +
  ggtitle(paste0(unique(gbif_data_d_sf$states_abbrev), ' species records \n across education ', unique(gbif_state_USCENSUS$species))) +  
  scale_fill_manual(values = col_pal) + theme_classic() + ylab('Frequency') + xlab('Educational level')

density_hist_age =   ggplot() +
  geom_density(aes(estimate,
                   fill = "GBIF records"),
               alpha = .2,
               data = gbif_edu_age[gbif_edu_age$variable =='median_age',])  +
  geom_density(aes(estimate, fill = paste0("State ", unique(gbif_data_d_sf$states_abbrev))),
               alpha = .2,
               data = state_edu_age[state_edu_age$variable =='median_age',]) +
  # ggtitle('Median age New Jersey - Lanternfly') +
  ggtitle(paste0(unique(gbif_data_d_sf$states_abbrev), ' species records \n across median age ', unique(gbif_state_USCENSUS$species))) +  
  scale_fill_manual(values = col_pal) + theme_classic() + ylab('Frequency') + xlab('Median age')

# 
# density_hist_bachelor =   ggplot() +
#   geom_density(aes(estimate,
#                    fill = "GBIF records"),
#                alpha = .2,
#                data = gbif_edu_age[gbif_edu_age$variable =='bachelor_degree',])  +
#   geom_density(aes(estimate, fill = paste0("State ", unique(gbif_data_d_sf$states_abbrev))), alpha = .2, data = state_edu_age[state_edu_age$variable =='bachelor_degree',]) +
#   ggtitle(paste0('Education ',paste0("State ", unique(gbif_data_d_sf$states_abbrev)),
#                  ' - ', unique(gbif_data_d_sf$species) )) +
#   scale_fill_manual(values = col_pal) + theme_classic() + ylab('Frequency') + xlab('Bachelor degree level')
# 

ggsave(paste0(
  "/Users/diegoellis/Desktop/age_education_",unique(gbif_state_USCENSUS$species),".jpeg"),
  (density_hist_edu / density_hist_age),
  "jpeg", height=12, width=12, units="in")


} # End of for looP:


download_plot_socioeco_gbif_data(gbif_species_name = 'Ailanthus altissima')
download_plot_socioeco_gbif_data(gbif_species_name = 'Agrilus planipennis')
download_plot_socioeco_gbif_data(gbif_species_name = 'Rhinella marina')
download_plot_socioeco_gbif_data(gbif_species_name = 'Halyomorpha halys')
download_plot_socioeco_gbif_data(gbif_species_name = 'Ixodes scapularis')
download_plot_socioeco_gbif_data(gbif_species_name = 'Aedes aegypti')




