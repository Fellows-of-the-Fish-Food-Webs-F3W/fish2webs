############
## main.R ##
############

## Goal: Define main functions to generate the F2W dataset containing curated individual 
## fish length measurements and associated food webs.

#' Generate validated individual fish size measurements from ASPE surveys
#'
#' @description
#' Extracts, filters, cleans, and harmonises fish survey data from the ASPE
#' database to produce a standardized table of individual fish sizes.
#'
#' The workflow:
#' 1. extracts the required survey metadata and measurement tables;
#' 2. filters operations to retain only eligible sampling events;
#' 3. standardises species identifiers and length types;
#' 4. converts fork length to total length when required;
#' 5. removes implausible measurements and outliers;
#' 6. validates batch-level records;
#' 7. reconstructs individual size measurements from the cleaned data.
#'
#' The resulting table is suitable for downstream size-based ecological and
#' food-web analyses.
#'
#' @param write_output Logical. If `TRUE`, writes the final individual-size
#'   table to `tab_fish_individual_size.csv` in the working directory.
#'
#' @return
#' A data frame of validated individual fish size measurements, typically
#' containing one row per individual.
#'
#' @details
#' This function is designed as a cleaning and harmonisation pipeline for ASPE
#' fish data. It retains only sampling operations matching the targeted
#' protocols, applies species and length-type sanitation, and enforces
#' biologically and structurally plausible measurements before returning a
#' consolidated individual-size table.
#'
#' @seealso [size2webs()]
#'
#' @export
fish2size <- function(write_output = F){
  
  ## 1. Initiate
  message("Extracting metadata...")
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
  message("Filtering sampling operations...")
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
  message("Harmonising fish length measurements...")
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
  message("Generating final length measurements...")
  fish_individual_size <- generate_individual_sizes(sanitized, verbose = TRUE)
  if (write_output == T){
    message("Writing outputs...")
    write.csv(fish_individual_size, "tab_fish_individual_size.csv", quote = F, row.names = F)
  }
  
  ## End
  message("Done!")
  return(fish_individual_size)
  
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
  colnames(ind_clean) <- c("operation_id", "batch_id", "species_code", "size")
  
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
    write.csv(size_classes, "tab_size_classes.csv", quote=F, row.names=F)
    write.csv(tab_metaweb_flattened, "tab_metaweb.csv", quote=F, row.names=F)
    write.csv(tab_local_foodwebs_flattened, "tab_local_foodwebs.csv", quote=F, row.names=F)
    write.csv(tab_local_foodweb_summary_metrics, "tab_local_foodwebs_summary_metrics.csv", quote=F, row.names=F)
  }
  
  ## End
  message("Done!")
  return(list("tab_size_classes" = size_classes, 
              "tab_metaweb_flattened" = tab_metaweb_flattened,
              "tab_local_foodwebs_flattened" = tab_local_foodwebs_flattened,
              "tab_local_foodwebs_summary_metrics" = tab_local_foodweb_summary_metrics
    )
  )
}

#
###
