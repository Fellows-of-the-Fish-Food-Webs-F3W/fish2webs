############
## main.R ##
############

## Goal: Define main functions to generate the F2W dataset containing curated individual
## fish length measurements and associated food webs.

#' Generate validated individual fish size measurements from ASPE surveys
#'
#' @description
#' Extracts, filters, cleans, and harmonises fish survey data from the ASPE
#' database to generate the main fish-related tables of the F3W dataset.
#'
#' The workflow:
#' 1. cleans and filters sampling operations and fish records;
#' 2. harmonises species codes and fish length measurements;
#' 3. converts fork length to total length when required;
#' 4. removes biologically implausible measurements and outliers;
#' 5. validates batch-level records;
#' 6. reconstructs individual fish sizes;
#' 7. converts individual lengths to body weights;
#' 8. computes community- and species-level metrics;
#' 9. generates site, operation, individual, and community tables.
#'
#' @param write_output Logical. If `TRUE`, writes the generated tables as CSV
#'   files in the working directory.
#'
#' @return
#' A named list containing five data frames:
#' \itemize{
#'   \item `site_information`
#'   \item `operation_information`
#'   \item `fish_individual_size_weight`
#'   \item `community_metrics`
#'   \item `species_level_metrics`
#' }
#'
#' @seealso [size2webs()]
#'
#' @export
fish2size <- function(write_output = FALSE) {

  ## 1. Prepare ASPE data
  message("Preparing ASPE data...")

  species <- cleaning_species_ref_aspe() |>
    tibble::as_tibble()

  ref_length <- cleaning_ref_length_type_aspe() |>
    tibble::as_tibble()

  operation <- clean_operation_aspe() |>
    tibble::as_tibble()

  ele_sampling <- cleaning_elementary_sampling() |>
    tibble::as_tibble()

  point_group <- cleaning_point_group()

  fish_batch <- clean_fish_batch()

  ind_measure <- clean_individual_measurement_aspe()

  station <- clean_station_aspe(
    station = get_raw_station_aspe(),
    ref_coordinates = get_raw_ref_coordinates_station_aspe(),
    crs_to = 4326
  )


  ## 2. Filter sampling operations
  message("Filtering sampling operations...")

  filtered_sampling <- filter_operation_batch_measure(
    operation = operation,
    op_protocol_to_keep = c(
      "complete",
      "partial_by_point",
      "partial_over_bank"
    ),
    op_objective_to_exclude = vec_op_objective_to_exclude(),
    oldest_sampling_date = "1995-01-01",
    omit_na_site = TRUE,
    point_group = point_group,
    min_prop_point_group_on_bank = 0.999,
    ele_sampling = ele_sampling,
    max_passage_number = 1,
    fish_batch = fish_batch,
    ind_measure = ind_measure
  )

  # Harmonise species codes in batch and individual-measurement data.
  filtered_sampling[c("fish_batch", "ind_measure")] <-
    filtered_sampling[c("fish_batch", "ind_measure")] |>
    purrr::map(\(x) sanitize_species_code(data = x))


  ## 3. Harmonise length measurements
  message("Cleaning and harmonising fish length measurements...")
  species_fork_length <- filtered_sampling$fish_batch |>
    dplyr::filter(length_type == "fork") |>
    dplyr::distinct(species_code) |>
    dplyr::left_join(
      species,
      by = dplyr::join_by(species_code)
    )

  length_length <- rfishbase::length_length(
    species_list = species_fork_length$latin_name
  )

  sanitized_length_type <- convert_fork_to_total(
    fish_batch = filtered_sampling$fish_batch,
    ind_measure = filtered_sampling$ind_measure,
    species_ref = cleaning_species_ref_aspe(),
    fishbase_length_length = length_length,
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

  sanitized <- sanitize_batch_data(
    fish_batch = sanitized_length$fish_batch,
    ind_measure = sanitized_length$ind_measure,
    species_ref = species,
    min_individuals_G = 2,
    min_individuals_SL = 6
  )


  ## 4. Generate individual fish data
  message("Generating individual fish sizes and weights...")

  fish_individual_size <- generate_individual_sizes(
    sanitized,
    verbose = TRUE
  )

  fish_individual_weight <- convert_length_to_weight(
    fish_data = fish_individual_size,
    species_ref = species,
    verbose = TRUE
  )


  ## 5. Compute community metrics
  message("Computing community metrics...")

  community_metrics <- compute_community_metrics(
    fish_data = fish_individual_weight$fish_data,
    operation = filtered_sampling$operation
  )


  ## 6. Generate output tables
  retained_operations <- unique(fish_individual_weight$fish_data$operation_id)

  tab3 <- filtered_sampling$operation |>
    dplyr::filter(operation_id %in% retained_operations) |>
    dplyr::distinct(site_id) |>
    dplyr::inner_join(
      station |>
        dplyr::select(site_id, x, y),
      by = "site_id"
    ) |>
    dplyr::arrange(site_id)

  tab4 <- filtered_sampling$operation |>
    dplyr::filter(operation_id %in% retained_operations) |>
    dplyr::mutate(
      time = format(date_time, "%H:%M:%S")
    ) |>
    dplyr::select(
      operation_id,
      site_id,
      date,
      time,
      protocol,
      computed_surface
    ) |>
    dplyr::arrange(operation_id)

  tab5 <- fish_individual_weight$fish_data |>
    dplyr::select(
      operation_id,
      species_code,
      size_mm,
      measured,
      weight_g
    )

  tab6 <- community_metrics |>
    dplyr::select(
      operation_id,
      total_richness,
      total_abundance,
      total_biomass_g,
      richness_per_m2,
      abundance_per_m2,
      biomass_g_per_m2
    )

  tab7_abundance <- community_metrics |>
    dplyr::select(
      operation_id,
      abundance_by_species
    ) |>
    tidyr::unnest(abundance_by_species)

  tab7_biomass <- community_metrics |>
    dplyr::select(
      operation_id,
      biomass_by_species
    ) |>
    tidyr::unnest(biomass_by_species)

  tab7 <- tab7_abundance |>
    dplyr::left_join(
      tab7_biomass,
      by = c("operation_id", "species_code")
    ) |>
    dplyr::select(
      operation_id,
      species_code,
      total_abundance,
      abundance_per_m2,
      total_biomass_g,
      biomass_g_per_m2
    ) |>
    dplyr::arrange(
      operation_id,
      species_code
    )


  ## 7. Write outputs
  if (write_output) {

    message("Writing outputs...")

    write.csv(
      tab3,
      "3_tab_site_information.csv",
      quote = FALSE,
      row.names = FALSE
    )

    write.csv(
      tab4,
      "4_tab_operation_information.csv",
      quote = FALSE,
      row.names = FALSE
    )

    write.csv(
      tab5,
      "5_tab_fish_individual_size_weight.csv",
      quote = FALSE,
      row.names = FALSE
    )

    write.csv(
      tab6,
      "6_tab_community_metrics.csv",
      quote = FALSE,
      row.names = FALSE
    )

    write.csv(
      tab7,
      "7_tab_species_level_metrics.csv",
      quote = FALSE,
      row.names = FALSE
    )
  }


  ## 8. Return
  message("Done!")

  invisible(
    list(
      site_information = tab3,
      operation_information = tab4,
      fish_individual_size_weight = tab5,
      community_metrics = tab6,
      species_level_metrics = tab7
    )
  )
}

#' Build size-structured food webs from individual fish sizes
#'
#' @description
#' Constructs a size-based global metaweb and operation-specific local food
#' webs from a table of individual fish size measurements.
#'
#' The full `ind_measure` dataset is used to build the metaweb. Local food webs
#' are then extracted for each sampling operation represented in
#' `ind_measure`, or for the subset specified by `selected_operations`.
#'
#' Individuals are assigned to size classes, trophic species are inferred from
#' species-specific diet shifts and predator-prey size constraints, and network
#' metrics are computed for each local food web.
#'
#' @param num_classes Integer. Number of size classes to use when constructing
#'   the size-structured web.
#' @param ind_measure A data frame of individual fish measurements. Must
#'   contain at least `operation_id`, `batch_id`, `species_code`, and `size`.
#' @param resource_diet_shift A data frame or table describing resource-use
#'   shifts used to parameterize resource nodes in the metaweb.
#' @param fish_diet_shift A data frame or table describing fish diet shifts used
#'   to parameterize trophic interactions in the metaweb.
#' @param pred_win Numeric or table defining the predator-prey size window used
#'   to constrain trophic links.
#' @param write_output Logical. If `TRUE`, writes the generated outputs to CSV
#'   files in the working directory.
#' @param selected_operations Optional character vector of `operation_id`
#'   values. If provided, local food webs are built only for these operations.
#'
#' @return
#' A named list with four elements:
#' \describe{
#'   \item{tab_size_classes}{A data frame of size-class boundaries for each
#'   species.}
#'   \item{tab_metaweb_flattened}{A flattened representation of the global
#'   metaweb.}
#'   \item{tab_local_foodwebs_flattened}{A flattened representation of all local
#'   food webs.}
#'   \item{tab_local_foodwebs_summary_metrics}{A data frame of summary network
#'   metrics for each local food web.}
#' }
#'
#' @details
#' The metaweb is always constructed from the full cleaned dataset, regardless
#' of any `selected_operations` subset. The `selected_operations` argument only
#' affects which local food webs are extracted and summarized.
#'
#' The function proceeds in four main stages:
#' 1. removes individuals with missing species information;
#' 2. computes size classes for all species present in the dataset;
#' 3. builds the global metaweb using diet-shift and size-constraint rules;
#' 4. extracts local food webs by operation and computes summary metrics.
#'
#' If `write_output = TRUE`, the following files are written:
#' `tab_size_classes.csv`, `tab_metaweb.csv`, `tab_local_foodwebs.csv`, and
#' `tab_local_foodwebs_summary_metrics.csv`.
#'
#' @seealso [compute_size_classes()], [build_metaweb()],
#'   [build_local_foodweb()]
#'
#' @export
size2webs <- function(num_classes, ind_measure, resource_diet_shift, fish_diet_shift, pred_win, write_output=T, selected_operations=NULL){

  ## Filter out missing species
  ind_clean <- remove_missing_species(
    ind_measure     = ind_measure, # INPUT
    fish_diet_shift = fish_diet_shift, # INPUT
    pred_win        = pred_win # INPUT
  )

  ## Clean column names
  #colnames(ind_clean) <- c("operation_id", "batch_id", "species_code", "size")

  ## Compute size classes
  message("Computing size classes...")
  size_classes <- compute_size_classes(
    ind_measure = ind_clean,
    num_classes = num_classes
  )

  ## Compute metaweb
  message("Building metaweb...")
  metaweb <- build_metaweb(
    tab_size_classes    = size_classes,
    pred_win            = pred_win, # INPUT
    fish_diet_shift     = fish_diet_shift, # INPUT
    resource_diet_shift = resource_diet_shift, # INPUT
    num_classes         = num_classes, # INPUT
    selected_resources  = c("det", "biof", "phytob", "macroph", "phytopl", "zoopl", "zoob") # WB: Make sure by default it's all resources
  )

  ## Subset operations
  if (!is.null(selected_operations)){
    s <- NULL; for (selected_operations_ in selected_operations) s <- c(s, which(ind_clean$operation_id == selected_operations_))
    ind_clean <- ind_clean[s,]
  }

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
    write.csv(size_classes, "8_tab_trophic_species_size_classes.csv", quote=F, row.names=F)
    write.csv(tab_metaweb_flattened, "9_tab_metaweb.csv", quote=F, row.names=F)
    write.csv(tab_local_foodwebs_flattened, "10_tab_local_foodwebs.csv", quote=F, row.names=F)
    write.csv(tab_local_foodweb_summary_metrics, "11_tab_local_foodweb_metrics.csv", quote=F, row.names=F)
  }

  ## End
  message("Done!")
  return(list("tab_trophic_species_size_classes" = size_classes,
              "tab_metaweb_flattened" = tab_metaweb_flattened,
              "tab_local_foodwebs_flattened" = tab_local_foodwebs_flattened,
              "tab_local_foodwebs_metrics" = tab_local_foodweb_summary_metrics
    )
  )
}

#
###
