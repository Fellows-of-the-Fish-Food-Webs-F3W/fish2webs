###############
## FUNCTIONS ##
###############

fish2size <- function(){
  
  ## 1. Initiate
  ref_protocol <- get_ref_protocol_operation_aspe() |> as_tibble()
  op_objective <- get_objective_operation_aspe() |> as_tibble()
  sampling_point <- get_sampling_point_aspe() |> as_tibble()
  op_description <- get_description_operation_aspe() |> as_tibble()
  ref_isolation <- get_ref_isolation_operation_aspe() |> as_tibble()
  raw_ref_species <- get_species_aspe() |> as_tibble()
  
  detail_sampling <- get_elementary_sampling_aspe() |> as_tibble()
  ref_detail_sampling <- get_ref_elementary_sampling_aspe() |> as_tibble()
  ref_prospection <- get_ref_prospection_method_operation_aspe() |> as_tibble()
  ref_passage <- get_ref_passage_aspe() |> as_tibble()
  
  point_group <- get_point_group_aspe() |> as_tibble()
  ref_point_group <- get_ref_point_group_aspe() |> as_tibble()
  
  species <- cleaning_species_ref_aspe() |> as_tibble()
  
  ref_batch_type <- get_ref_type_batch_aspe() |> as_tibble()
  raw_fish_batch <- get_fish_batch_aspe() |> as_tibble()
  raw_individual_measurement <- get_individual_measurement_aspe() |> as_tibble()
  raw_ref_length <- get_ref_type_length_aspe() |> as_tibble()
  ref_length <- cleaning_ref_length_type_aspe() |> as_tibble()
  
  op <- clean_operation_aspe() %>%
    as_tibble() %>%
    filter(
      protocol %in% c("complete", "partial_by_point", "partial_over_bank")
    )
  
  ele_sampling <- cleaning_elementary_sampling() |> as_tibble()
  point_group <- cleaning_point_group()
  
  fish_batch <- clean_fish_batch(fish_batch = raw_fish_batch)
  ind_measure <- clean_individual_measurement_aspe(
    ind_measure = raw_individual_measurement
  )
  
  ## 2. Filter out operations
  filtered_sampling <- filter_operation_batch_measure(
    operation = op,
    op_protocol_to_keep = c("complete", "partial_by_point",  "partial_over_bank"),
    op_objective_to_exclude = vec_op_objective_to_exclude(),
    oldest_sampling_date = "1995-01-01",
    omit_na_site = TRUE,
    point_group = point_group,
    min_prop_point_group_on_bank = .999,
    ele_sampling = ele_sampling,
    max_passage_number = 1,
    fish_batch = fish_batch,
    ind_measure = ind_measure
  )
  
  ## TODO: adding fishbase validated names in cleaning_species_ref_aspe()
  filtered_sampling[c("fish_batch", "ind_measure")] <-
    filtered_sampling[c("fish_batch", "ind_measure")] |>
    purrr::map(\(x) sanitize_species_code(data = x))
  summary(as.factor(filtered_sampling$fish_batch$length_type))

  ## 3. Harmonise length measurements
  species_fork_length <- filtered_sampling$fish_batch |>
    filter(length_type == "fork") |>
    distinct(species_code) |>
    left_join(species,
              by = join_by(species_code)
    )
  ll <- rfishbase::length_length(species_list = species_fork_length$latin_name)
  
  sanitized_length_type <- convert_fork_to_total(
    fish_batch = filtered_sampling$fish_batch,
    ind_measure = filtered_sampling$ind_measure,
    species_ref = cleaning_species_ref_aspe(),
    fishbase_length_length = ll,
    conversion_vector = coefficients_fork2total(),
    convert_intercept_cm2mm = TRUE,
    manual_priority = FALSE,
    verbose = TRUE
  )
  
  sanitized_length <- remove_impossible_lengths(
    ind_measure = sanitized_length_type$ind_measure,
    fish_batch = sanitized_length_type$fish_batch,
    species_ref = species,
    remove_outliers = TRUE
  )
  
  sanitized_length$outlier_summary |>
    arrange(desc(n_outliers))
  sanitized_length$outliers |>
    filter(species_code == "HOT") |>
    select(-c(operation_id:measure_id)) |>
    arrange(desc(size)) |>
    print(n = 200)
  
  sanitized <- sanitize_batch_data(
    fish_batch = sanitized_length$fish_batch,
    ind_measure = sanitized_length$ind_measure,
    species_ref = species,
    min_individuals_G = 2,
    min_individuals_SL = 6
  )
  
  sanitized$filtering_log
  sanitized$validation_issues
  sanitized$fish_batch |>
    filter(min_length == max_length) |>
    select(batch_type, min_length, max_length, number)
  sanitized$fish_batch |>
    filter(batch_type == "G") |>
    filter(is.na(min_length) | is.na(max_length))
  

  ## 4. Final file
  fish_individual_size <- generate_individual_sizes(sanitized, verbose = TRUE)
  if (write_output == T){
    write.csv(fish_individual_size, "tab_fish_individual_size.csv", quote = F, row.names = F)
  }
  
  ## End
  return(fish_individual_size)
  
}

#
###

##########
## MAIN ##
##########

##############
## INITIATE ##
##############

## Goal: Install fishdatabuilder and foodwebdatabuilder

library(pkgfilecache)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(devtools)
devtools::load_all()
document()

## General packages
library(dplyr)

## F3W packages
library(fishdatabuilder)
library(foodwebbuilder)
#
###

###################
## DATA DOWNLOAD ##
###################

## Goal: Download and list data files

## Download the data
download_optional_data()

## Extract the data
extract_optional_data()

library(devtools)
library(pkgfilecache)
load_all()
#> ℹ Loading fishdatabuilder
document()
#> ℹ Updating fishdatabuilder documentation
#> ℹ Loading fishdatabuilder

#
###

####################
## MAIN FISH2SIZE ##
####################

fish_individual_size <- fish2size(write_output = T)

#
###

# #####
# ## ##
# #####
# 
# list_optional_files()
# 
# # op <- get_raw_operation_aspe() |> as_tibble() # OLD
# # ref_objective <- get_ref_objective_operation_aspe() |> as_tibble() # OLD?
# ref_protocol <- get_ref_protocol_operation_aspe() |> as_tibble()
# op_objective <- get_objective_operation_aspe() |> as_tibble()
# sampling_point <- get_sampling_point_aspe() |> as_tibble()
# op_description <- get_description_operation_aspe() |> as_tibble()
# ref_isolation <- get_ref_isolation_operation_aspe() |> as_tibble()
# raw_ref_species <- get_species_aspe() |> as_tibble()
# 
# detail_sampling <- get_elementary_sampling_aspe() |> as_tibble()
# ref_detail_sampling <- get_ref_elementary_sampling_aspe() |> as_tibble()
# ref_prospection <- get_ref_prospection_method_operation_aspe() |> as_tibble()
# ref_passage <- get_ref_passage_aspe() |> as_tibble()
# 
# point_group <- get_point_group_aspe() |> as_tibble()
# ref_point_group <- get_ref_point_group_aspe() |> as_tibble()
# 
# species <- cleaning_species_ref_aspe() |> as_tibble()
# 
# ref_batch_type <- get_ref_type_batch_aspe() |> as_tibble()
# raw_fish_batch <- get_fish_batch_aspe() |> as_tibble()
# raw_individual_measurement <- get_individual_measurement_aspe() |> as_tibble()
# raw_ref_length <- get_ref_type_length_aspe() |> as_tibble()
# ref_length <- cleaning_ref_length_type_aspe() |> as_tibble()
# 
# op <- clean_operation_aspe() %>%
#   as_tibble() %>%
#   filter(
#     protocol %in% c("complete", "partial_by_point", "partial_over_bank")
#   )
# head(op)
# 
# ele_sampling <- cleaning_elementary_sampling() |> as_tibble()
# point_group <- cleaning_point_group()
# 
# fish_batch <- clean_fish_batch(fish_batch = raw_fish_batch)
# ind_measure <- clean_individual_measurement_aspe(
#   ind_measure = raw_individual_measurement
# )
# 
# ## Check
# head(fish_batch)
# head(ind_measure)
# 
# #
# ###
# 
# #####
# ## ##
# #####
# 
# operation_id_sample <- sample(op$operation_id, 100)
# head(operation_id_sample)
# 
# batch_test <- fish_batch %>%
#   filter(operation_id %in% operation_id_sample) %>%
#   select(-c(
#     weight, estimated_weight, mep_id,
#     prelevement_id, operation_id
#   ))
# 
# measure_test <- ind_measure %>%
#   filter(operation_id %in% operation_id_sample) %>%
#   select(-c(site_id:prelevement_id, measure_id))
# 
# #
# ###
# 
# #####
# ## ##
# #####
# 
# filtered_sampling <- filter_operation_batch_measure(
#   operation = op,
#   op_protocol_to_keep = c("complete", "partial_by_point",  "partial_over_bank"),
#   op_objective_to_exclude = vec_op_objective_to_exclude(),
#   oldest_sampling_date = "1995-01-01",
#   omit_na_site = TRUE,
#   point_group = point_group,
#   min_prop_point_group_on_bank = .999,
#   ele_sampling = ele_sampling,
#   max_passage_number = 1,
#   fish_batch = fish_batch,
#   ind_measure = ind_measure
# )
# unique(filtered_sampling$operation$site_id)
# 
# # Sanitize
# ## TODO: adding fishbase validated names in cleaning_species_ref_aspe()
# filtered_sampling[c("fish_batch", "ind_measure")] <-
#   filtered_sampling[c("fish_batch", "ind_measure")] |>
#   purrr::map(\(x) sanitize_species_code(data = x))
# summary(as.factor(filtered_sampling$fish_batch$length_type))
# # Check
# ## TODO: check on validated names, i.e. output different than entry
# filtered_sampling[c("fish_batch", "ind_measure")] |>
#   purrr::map(\(x) check_aspe_fish_species(species = x$species_code)
#   )
# 
# #
# ###
# 
# #####
# ## ##
# #####
# 
# species_fork_length <- filtered_sampling$fish_batch |>
#   filter(length_type == "fork") |>
#   distinct(species_code) |>
#   left_join(species,
#             by = join_by(species_code)
#   )
# ll <- rfishbase::length_length(species_list = species_fork_length$latin_name)
# 
# sanitized_length_type <- convert_fork_to_total(
#   fish_batch = filtered_sampling$fish_batch,
#   ind_measure = filtered_sampling$ind_measure,
#   species_ref = cleaning_species_ref_aspe(),
#   fishbase_length_length = ll,
#   conversion_vector = coefficients_fork2total(),
#   convert_intercept_cm2mm = TRUE,
#   manual_priority = FALSE,
#   verbose = TRUE
# )
# 
# sanitized_length <- remove_impossible_lengths(
#   ind_measure = sanitized_length_type$ind_measure,
#   fish_batch = sanitized_length_type$fish_batch,
#   species_ref = species,
#   remove_outliers = TRUE
# )
# 
# sanitized_length$outlier_summary |>
#   arrange(desc(n_outliers))
# sanitized_length$outliers |>
#   filter(species_code == "HOT") |>
#   select(-c(operation_id:measure_id)) |>
#   arrange(desc(size)) |>
#   print(n = 200)
# 
# sanitized <- sanitize_batch_data(
#   fish_batch = sanitized_length$fish_batch,
#   ind_measure = sanitized_length$ind_measure,
#   species_ref = species,
#   min_individuals_G = 2,
#   min_individuals_SL = 6
# )
# 
# sanitized$filtering_log
# sanitized$validation_issues
# sanitized$fish_batch |>
#   filter(min_length == max_length) |>
#   select(batch_type, min_length, max_length, number)
# sanitized$fish_batch |>
#   filter(batch_type == "G") |>
#   filter(is.na(min_length) | is.na(max_length))
# 
# #
# ###
# 
# #####
# ## ##
# #####
# 
# fish_individual_size <- generate_individual_sizes(sanitized, verbose = TRUE)
# head(fish_individual_size)
# 
# fish_individual_size |>
#   filter(species_code == "CHP")
# # write_csv(
# #   x = fish_individual_size,
# #   file = here::here("data-raw", "fish_individual_size.csv")
# # )
# 
# head(sanitized)
# 
# #
# ###

###############################################################################
##                                                                           ##
##                            FOOD WEB BUILDER                               ##
##                                                                           ##
###############################################################################

## inputs
num_classes <- 9
ind_measure <- fish_individual_size

size2webs <- function(num_classes, ind_measure, resource_diet_shift, fish_diet_shift, pred_win, write_output=T){
  
  ## Filter out missing species
  ind_clean <- remove_missing_species(
    ind_measure     = ind_measure, # INPUT
    fish_diet_shift = fish_diet_shift, # INPUT
    pred_win        = pred_win # INPUT
  )
  
  ## Clean column names
  colnames(ind_clean) <- c("operation_id", "batch_id", "species_code", "size")
  
  ## Compute size classes
  message("Computing size classes...")
  size_classes <- compute_size_classes(
    ind_measure = ind_clean,
    num_classes = NUM_CLASSES
  )
  
  ## Compute metaweb
  message("Building metaweb...")
  metaweb <- build_metaweb(
    tab_size_classes    = size_classes, 
    pred_win            = pred_win, # INPUT
    fish_diet_shift     = fish_diet_shift, # INPUT
    resource_diet_shift = resource_diet_shift, # INPUT
    num_classes         = NUM_CLASSES, # INPUT
    selected_resources  = c("det", "biof", "phytob", "macroph", "phytopl", "zoopl", "zoob") # WB: Make sure by default it's all resources
  )
  
  ## Only for debug purposes sample operations 
  selected_op <- unique(ind_clean$operation_id)[1:1000] # DEBUG
  s <- NULL; for (selected_op_ in selected_op) s <- c(s, which(ind_clean$operation_id == selected_op_))
  ind_clean <- ind_clean[s,]
  
  ## Build local foodwebs
  message("Extracting local food webs...")
  local_foodwebs <- build_local_foodweb(
    ind_measure       = ind_clean,
    local_id          = "operation_id",         # column in ind_measure
    metaweb           = metaweb,
    tab_size_classes  = size_classes,
    selected_resources  = c("det", "biof", "phytob", "macroph", "phytopl", "zoopl", "zoob")
  )
  
  ## Compute local food web metrics
  message("Compute local food web metrics...")
  tab_local_foodweb_summary_metrics <- NULL
  for (local_foodweb in local_foodwebs){
    metrics <- compute_metrics_summary(local_foodweb)
    tab_local_foodweb_summary_metrics <- rbind(tab_local_foodweb_summary_metrics, metrics)
  }
  tab_local_foodweb_summary_metrics <- data.frame(tab_local_foodweb_summary_metrics)
  tab_local_foodweb_summary_metrics$operation_id <- names(local_foodwebs)
  
  ## Flatten food webs and store
  message("Flattening local food webs for storage...")
  tab_metaweb_flattened <- flatten_foodweb(metaweb)
  tab_local_foodwebs_flattened <- flatten_foodweb_list(local_foodwebs)
  
  ## Write output
  if (write_output == T){
    message("Writing outputs...")
    write.csv(size_classes, "tab_size_classes.csv", quote=F, row.names=F)
    write.csv(tab_metaweb_flattened, "tab_metaweb.csv", quote=F, row.names=F)
    write.csv(tab_local_foodwebs_flattened, "tab_local_foodwebs.csv", quote=F, row.names=F)
    write.csv(tab_local_foodweb_summary_metrics, "tab_local_foodwebs_summary_metrics.csv", quote=F, row.names=F)
  }
  
  ## End
  return(list("tab_size_classes" = size_classes, 
              "tab_metaweb_flattened" = tab_metaweb_flattened,
              "tab_local_foodwebs_flattened" = tab_local_foodwebs_flattened,
              "tab_local_foodwebs_summary_metrics" = tab_local_foodweb_summary_metrics
              )
         )
}

## Test function
outputs <- size2webs(num_classes, ind_measure, resource_diet_shift, fish_diet_shift, pred_win, write_output=F)

##############################
## EXAMPLE ANALYSIS METAWEB ##
##############################

## Rebuild metaweb
metaweb <- unflatten_foodweb(outputs$tab_metaweb_flattened)
# local_foodwebs <- unflatten_foodweb(outputs$tab_metaweb)

## Structural metrics
basal_nodes <- get_basal_nodes(metaweb)
leaf_nodes  <- get_leaf_nodes(metaweb)

degree_in  <- compute_inward_degree(metaweb)
degree_out <- compute_outward_degree(metaweb)

TL <- compute_trophic_level(metaweb)
TB <- compute_trophic_breadth(metaweb, TL)

## Fluxes (matrix of accumulated bottom-up fluxes)
fluxes <- compute_bottom_up_fluxes(metaweb)

incoming_flux <- colSums(fluxes)  # received by consumers (columns)
outgoing_flux <- rowSums(fluxes)  # emitted from prey/resources (rows)

##
plot_network(
  metaweb,
  x = TB,
  y = TL,
  labels = colnames(metaweb),
  xlab = "Trophic breadth (sd of resource TL)",
  ylab = "Trophic level",
  add_legend = "bottomright",
  label_space_x = 0.10
)

plot_network_radial(
  fluxes,
  x = TL,
  y = TL,
  labels = colnames(metaweb),
  rotations = 1.5,
  scale = 1.25,
  line_width_max = 1
)

plot_network_radial(
  fluxes,
  x = degree_in,
  labels = colnames(metaweb),
  rotations = 0.5,
  scale = 1.25,
  line_width_max = 1
)

#
###

####################################
## EXAMPLE ANALYSIS LOCAL FOODWEB ##
####################################

## Select a local food web
local_foodwebs <- outputs$tab_local_foodwebs_flattened
s <- which(local_foodwebs$operation_id == unique(local_foodwebs$operation_id)[2])
local_foodweb <- unflatten_foodweb(tab_local_foodwebs_flattened[s,])

## Structural metrics
basal_nodes <- get_basal_nodes(local_foodweb)
leaf_nodes  <- get_leaf_nodes(local_foodweb)

degree_in  <- compute_inward_degree(local_foodweb)
degree_out <- compute_outward_degree(local_foodweb)

TL <- compute_trophic_level(local_foodweb)
TB <- compute_trophic_breadth(local_foodweb, TL)

## Fluxes (matrix of accumulated bottom-up fluxes)
fluxes <- compute_bottom_up_fluxes(local_foodweb)

incoming_flux <- colSums(fluxes)  # received by consumers (columns)
outgoing_flux <- rowSums(fluxes)  # emitted from prey/resources (rows)

##
plot_network(
  local_foodweb,
  x = TB,
  y = TL,
  labels = colnames(local_foodweb),
  xlab = "Trophic breadth (sd of resource TL)",
  ylab = "Trophic level",
  add_legend = "bottomright",
  label_space_x = 0.10
)

plot_network_radial(
  fluxes,
  x = TL,
  y = TL,
  labels = colnames(local_foodweb),
  rotations = 1.5,
  scale = 1.25,
  line_width_max = 1
)

plot_network_radial(
  fluxes,
  x = degree_in,
  labels = colnames(local_foodweb),
  rotations = 0.5,
  scale = 1.25,
  line_width_max = 1
)

#
###