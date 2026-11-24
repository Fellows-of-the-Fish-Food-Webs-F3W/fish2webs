##########
## MAIN ##
##########

## Goal: Install fishdatabuilder and foodwebdatabuilder

###############
## FUNCTIONS ##
###############

plot_network <- function (M, x = NULL, y = NULL, labels = NULL, xlab = "", ylab = "", 
                  line_width_max = 1, label_space_x = 0.05, add_legend = "topright",
                  cex.dots=2, cex.text=1.5) 
{
  d = dim(M)[1]
  if (is.null(labels) == T) 
    labels = 1:d
  if (is.null(x) == T) 
    x = cos((1:d)/0.25)
  if (is.null(y) == T) 
    y = sin((1:d)/0.25)
  delta_x = (max(x) - min(x)) * label_space_x
  plot(x, y, xlim = c(min(x) - delta_x, max(x) + delta_x * 
                        1.5), cex = 0, bty = "n", xlab = xlab, ylab = ylab, cex.lab = 1.5, 
       bty = "l")
  for (j in 1:ncol(M)) {
    color_ = rainbow(ncol(M))[j]
    for (i in 1:nrow(M)) {
      line_width = abs(M[i, j])/max(abs(M)) * line_width_max
      lines(x = c(x[i], x[j] - delta_x/2), y = c(y[i], 
                                                 y[j]), col = color_, lwd = line_width)
      # arrows(x0 = x[i], x1 = x[j] - delta_x/2, y0 = y[i], 
      #        y1 = y[j], col = color_, lwd = line_width)
    }
  }
  points(x - delta_x/2, y, pch = 1, cex = cex.dots)
  points(x, y, pch = 16, cex = cex.dots)
  text(x + delta_x, y, labels = labels, cex = cex.text)
  if (add_legend != "off") {
    legend(add_legend, legend = c("Ingoing effects", "Outgoing effects"), 
           pch = c(1, 16), col = c("black", "black"), cex = 1.5, 
           bg = "white", box.col = "white")
  }
}

#
###

#########
## DEV ##
#########

remove.packages("fish2webs")
remove.packages("fishdatabuilder")
remove.packages("foodwebbuilder")
.rs.restartR()
devtools::document()
devtools::build()
devtools::install()
devtools::load_all()
.rs.restartR()
library(fish2webs)

#
###

###########
## SETUP ##
###########

install.packages("pak")
pak::pak("Fellows-of-the-Fish-Food-Webs-F3W/fish2webs")
library(fish2webs)
library(fishdatabuilder)

## Load data
data(resource_diet_shift, package = "foodwebbuilder")
data(fish_diet_shift, package = "foodwebbuilder")
data(pred_win, package = "foodwebbuilder")

#
###

###############
## LOAD DATA ##
###############

tab_fish_individual_size <- read.csv("tab_fish_individual_size.csv")
head(tab_fish_individual_size)


#
###

####################
## MAIN FISH2SIZE ##
####################

library(dplyr)

write_output <- F

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

#
###

############
## CHECKS ##
############

## To do:
## Implement retaining of site_id in ind_clean

## Checks
colnames(fish_individual_size)
colnames(op_description)
colnames(ind_measure)

## operation id --> site id --> coordinates
## Note: this should start from ind_clean.
## So either site_id needs to be kept in ind_clean or 
## the mapping operation id to site id needs
## to be unique. 
## Question: is the operation id is only valid within a site?
op_id_ <- op$operation_id[300]
s <-  which(op$operation_id == op_id_)
site_id_ <- unique(op$site_id[s])
if (length(site_id_)>1){
  message("Problem: op_id maps to multiple site_id")
} else {
  s <- which(sampling_point$pop_sta_id == site_id_)
  coords_x_ <- sampling_point$pop_coordonnees_x[s]
  coords_y_ <- sampling_point$pop_coordonnees_y[s]
  plot(coords_x_,coords_y_)
}
## Problem: this maps to different coordinates, so it cannot be the same site.
## Solution: switch to pop_id instead of pop_sta_id
op_id_ <- op$operation_id[300]
s <-  which(op$operation_id == op_id_)
site_id_ <- unique(op$site_id[s]) # Problem op_id maps to multiple site_id
if (length(site_id_)>1){
  message("Problem: op_id maps to multiple site_id")
} else {
  s <- which(sampling_point$pop_id == site_id_)
  coords_x_ <- sampling_point$pop_coordonnees_x[s]
  coords_y_ <- sampling_point$pop_coordonnees_y[s]
  plot(coords_x_,coords_y_)
}

## site_id --> coordinates
## If operation id is only valid within a site, we start instead with the site.
site_id_ <- ind_measure$site_id[100]
s <- which(sampling_point$pop_sta_id == site_id_)
coords_x_ <- sampling_point$pop_coordonnees_x[s]
coords_y_ <- sampling_point$pop_coordonnees_y[s]
plot(coords_x_,coords_y_)
## This maps to different coordinates, so it cannot be the same site.
site_id_ <- ind_measure$site_id[100]
s <- which(sampling_point$pop_id == site_id_)
coords_x_ <- sampling_point$pop_coordonnees_x[s]
coords_y_ <- sampling_point$pop_coordonnees_y[s]
plot(coords_x_,coords_y_)
## This maps to a single pair of coordinates.

## Efficient acquisition of site coordinates
unique_site_ids <- unique(sampling_point$pop_id)
s <- match(unique_site_ids, sampling_point$pop_id)
unique_site_ids_and_coords <- sampling_point[s,]
length(unique_site_ids)
nrow(unique_site_ids_and_coords)
get_site_coords <- function(site_id, unique_site_ids_and_coords=unique_site_ids_and_coords){
  s <- which(unique_site_ids_and_coords$pop_id == site_id)
  coords_x_ <- unique_site_ids_and_coords$pop_coordonnees_x[s]
  coords_y_ <- unique_site_ids_and_coords$pop_coordonnees_y[s]
  return(c(coords_x_, coords_y_))
}

## Efficient acquisition of times
unique_op_ids <- unique(op$operation_id)
s <- match(unique_op_ids, op$operation_id)
unique_op_ids_and_times <- op[s,]
get_op_time <- function(op_id, op_ids_and_times=unique_op_ids_and_times){
  s <- which(op_ids_and_times$operation_id == op_id)
  return(op_ids_and_times$date[s])
}
get_op_time(unique_op_ids[5])

#
###

######################################################
## EXAMPLE ANALYSIS: SINGLE-SITE COMMUNITY PATTERNS ##
######################################################

## All species
site_id_ <- ind_measure$site_id[9000]
site_ind_measure <- ind_measure[which(ind_measure$site_id == site_id_),]
res <- table(site_ind_measure$species_code)
colors <- rev(rainbow(length(res), start = 0.0, end = 0.8, alpha=0.75))
barplot(rev(res), horiz=T, las=1, col=colors, family = "mono")

## Selected species
selected_species <- c("BRO", "SIL", "PER", "TRF")
s <- NULL
for (selected_species_ in selected_species){
  s <- c(s, which(site_ind_measure$species_code ==  selected_species_))
}
subset_site_ind_measure <- site_ind_measure[s,]
res <- table(subset_site_ind_measure$species_code)
colors <- rainbow(length(res))
barplot(res, horiz=T, las=1, col=colors)

## Species abundances through operation
selected_species <- sort(unique(site_ind_measure$species_code)) # c("BRO", "SIL", "PER", "TRF", "GOU")
unique_op_ids <- unique(site_ind_measure$operation_id)
res <- NULL
for (unique_op_ids_ in unique_op_ids){
  s <- which(site_ind_measure$operation_id == unique_op_ids_)
  site_ind_measure_ <- site_ind_measure[s,]
  res_ <- NULL
  for (selected_species_ in selected_species){
    individual_count <- length(which(site_ind_measure_$species_code ==  selected_species_))
    res_ <- c(res_, individual_count)
  }
  res <- rbind(res, res_)
}
tab_species_abundances <- data.frame(res)
colnames(tab_species_abundances) <- selected_species
barplot(tab_species_abundances[,"PER"], col="lightgreen", names.arg = as.factor(unique_op_ids))
barplot(tab_species_abundances[,"BRO"], col="gold", add=T)
barplot(tab_species_abundances[,"TRF"], col="salmon", add=T)

## Sepecies abundances time series
times <- NULL 
for (unique_op_ids_ in unique_op_ids){
  times <- c(times, get_op_time(unique_op_ids_))
}
times <- as.Date(times)
times <- as.numeric(times - min(times)) # Time in days
s <- order(times)
x <- times[s]/365
y <- tab_species_abundances[s,]
plot(x, y$BRO , type="l")
plot(x, y$PER , type="l")
plot(x, y$TRF , type="l")
plot(x, y$SIL , type="l")
plot(x, y$GOU , type="l")

pdf(file="plot-single-site-abundances.pdf", width=5, height=4)
par(mar = c(5, 5, 1, 1), xpd = TRUE, mfrow=c(1,1), bty="l")
#
## Minimal stacked area plot
sp <- selected_species# c("BRO", "PER", "TRF", "SIL", "GOU")
## Convert abundances to proportions
p <- as.matrix(y[, sp])
p <- p / rowSums(p)
## Cumulative proportions for stacked areas
cp <- t(apply(p, 1, cumsum))
plot(x, x,
     type = "n",
     ylim = c(0, 1),
     xlab = "Time (years)",
     ylab = "Proportion",
     axes = FALSE)
# cols <- gray.colors(length(sp), start = 0.2, end = 0.85)
cols <- rainbow(length(sp), start = 0.0, end = 0.8, alpha=0.75)
for (i in seq_along(sp)) {
  lower <- if (i == 1) rep(0, length(x)) else cp[, i - 1]
  upper <- cp[, i]
  
  polygon(c(x, rev(x)),
          c(lower, rev(upper)),
          col = cols[i],
          border = NA
          )
}
axis(1)
axis(2, las = 1)
box()
# legend("right",
#        inset = c(-0.2, -0.0),
#        legend = sp,
#        fill = cols,
#        border = NA,
#        cex=0.7,
#        bty = "n",
#        xpd = T)
#
dev.off()
#
pdf(file="plot-single-site-counts.pdf", width=2.5, height=4)
par(mar = c(5, 5, 1, 1), xpd = TRUE, mfrow=c(1,1), bty="l")
#
## Counts
res <- table(site_ind_measure$species_code)
colors <- rainbow(length(res), start = 0.0, end = 0.8, alpha=0.75)
barplot(res, horiz=T, las=1, col=colors, xlab="Total count", cex.names=0.7)
#
dev.off()

## Visualise networks
k <- 1
order_ <- order(times)
for (k in 1:length(unique_op_ids)){
  
  pdf(file=paste("plot-single-site-network-",k,".pdf",sep=""), width=5, height=4)
  par(mar = c(5, 5, 1, 1))
  
  ## Check time
  unique_op_ids_ <- unique_op_ids[order_[k]]
  print(times[order_[k]]/365)
  
  ## Build local foodwebs
  subset_ind_clean <- ind_clean[which(ind_clean$operation_id == unique_op_ids_),]
  message("Extracting local food webs...")
  local_foodwebs <- build_local_foodweb(
    ind_measure       = subset_ind_clean,
    local_id          = "operation_id",         # column in ind_measure
    metaweb           = metaweb,
    tab_size_classes  = size_classes,
    selected_resources  = c("det", "biof", "phytob", "macroph", "phytopl", "zoopl", "zoob")
  )
  
  ##
  local_foodweb <- local_foodwebs[[1]]
  
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
  noise_x <- rnorm(length(TB),0,0.1)
  noise_y <- rnorm(length(TL),0,0.1)
  plot_network(
    local_foodweb,
    x = TB+noise_x,
    y = TL+noise_y,
    labels = colnames(local_foodweb),
    xlab = "Trophic breadth",
    ylab = "Trophic level",
    add_legend = "off",
    label_space_x = 0.075,
    line_width_max = 2,
    cex.dots=1,
    cex.text=0.5
  )
  # plot_network_radial(
  #   fluxes,
  #   x = degree_in,
  #   labels = colnames(local_foodweb),
  #   rotations = 0.5,
  #   scale = 1.25,
  #   line_width_max = 2
  # )
  
  ##
  legend("bottomright", legend = paste("t =",round(times[order_[k]]/365),"years"), bty="n", cex=1.5)
  # legend(-0.55, -0.25, legend = paste("t =",round(times[order_[k]]/365),"years"), bty="n", cex=2)
  
  ##
  dev.off()
  
}


# plot_network_radial(
#   fluxes,
#   x = TL,
#   y = TL,
#   labels = colnames(local_foodweb),
#   rotations = 1.5,
#   scale = 1.25,
#   line_width_max = 2
# )
# 
# plot_network_radial(
#   fluxes,
#   x = degree_in,
#   labels = colnames(local_foodweb),
#   rotations = 0.5,
#   scale = 1.25,
#   line_width_max = 2
# )

#
###

########################
## SPECIES MULTI-SITE ##
########################

pdf(file="plot-multi-site.pdf", width=8, height=8)
par(mar = c(4, 4, 4, 4), xpd = TRUE)
selected_species <- c("TRF", "BRO", "SIL")
colors <- rainbow(length(selected_species))
# selected_species_ <- selected_species[1]
k <- 1
res_range <- list()
for (selected_species_ in selected_species){
  species_ind_measure <-ind_measure[which(ind_measure$species_code == selected_species_),]
  nrow(species_ind_measure)
  res <- table(species_ind_measure$site_id)
  res_range[[selected_species_]] <- range(res)
  length(res)
  s <- match(names(res), unique_site_ids_and_coords$pop_id)
  coords <- unique_site_ids_and_coords[s,]
  head(coords)
  s <- which((coords$pop_coordonnees_x >= 4e+05) & 
               # (coords$pop_coordonnees_x < 8e+06) & 
               (coords$pop_coordonnees_y >= 6.65e+06) & 
               (coords$pop_coordonnees_y < 8e+06)
             )
  coords <- coords[s,]
  res <- res[s]
  u <- res / max(res) # u <- (res-min(res))/(max(res)-min(res))
  # colors <- rainbow(100)[round(u*100)]
  if (k == 1){
    plot(coords$pop_coordonnees_y, -coords$pop_coordonnees_x, col=colors[k], 
         pch=16, cex=2*u^0.25, 
         xaxt="n", yaxt="n", xlab="", ylab="")
  } else {
    points(coords$pop_coordonnees_y, -coords$pop_coordonnees_x, col=colors[k], 
           pch=16, cex=2*u^0.25)
  }
  ## Specific site
  s <- which(sampling_point$pop_id == site_id_)
  # points(sampling_point$pop_coordonnees_y[s], -sampling_point$pop_coordonnees_x[s], col="black", cex=3, pch=16)
  points(sampling_point$pop_coordonnees_y[s], -sampling_point$pop_coordonnees_x[s], col="black", cex=3, pch=1, lwd=8)
  points(sampling_point$pop_coordonnees_y[s], -sampling_point$pop_coordonnees_x[s], col="white", cex=3, pch=1, lwd=3)
  #
  k <- k + 1
}
legend_prop <- c(0.01, 0.1, 0.5, 1)
legend("bottomright",
       inset = c(-0.0, 0.0),
       legend = c("Species", selected_species,
                  "", "Prop. of max abundance",
                  paste0(100 * legend_prop, "%")),
       col = c(NA, colors, NA, NA, rep("black", length(legend_prop))),
       pch = c(NA, rep(16, length(selected_species)),
               NA, NA, rep(16, length(legend_prop))),
       pt.cex = c(NA, rep(1.5, length(selected_species)),
                  NA, NA, 2 * legend_prop^0.25),
       bty = "n",
       xpd = TRUE)
dev.off()

#
###

####################
## MAIN SIZE2WEBS ##
####################

library(foodwebbuilder)

num_classes = 9
write_output=F

## Filter out missing species
ind_clean <- remove_missing_species(
  ind_measure     = ind_measure, # INPUT
  fish_diet_shift = fish_diet_shift, # INPUT
  pred_win        = pred_win # INPUT
)

## Clean column names
colnames(ind_clean) <- c("operation_id", "species_code", "size")

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
outputs <- list("tab_size_classes" = size_classes, 
            "tab_metaweb_flattened" = tab_metaweb_flattened,
            "tab_local_foodwebs_flattened" = tab_local_foodwebs_flattened,
            "tab_local_foodwebs_summary_metrics" = tab_local_foodweb_summary_metrics
            )

#
###

# ####################
# ## MAIN FISH2WEBS ##
# ####################
# 
# ## Format individual fish sizes
# fish_individual_size <- fish2size(write_output = F)
# 
# ## Generate food webs from size for a few selected operations
# ## For demo purposes
# selected_op <- unique(fish_individual_size$operation_id)[1:100]
# outputs <- size2webs(num_classes = 3, fish_individual_size, resource_diet_shift, fish_diet_shift, pred_win, write_output=F, selected_op)
# 
# ## Generate food webs from size for all operations
# ## Takes a while so not for demo
# outputs <- size2webs(num_classes = 9, fish_individual_size, resource_diet_shift, fish_diet_shift, pred_win, write_output=F)
# 
# #
# ###

#########################################
## EXAMPLE ANALYSIS SIZE DISTRIBUTIONS ##
#########################################

## Number of panels in figure rows and columns
d <- 4

## Get species codes
unique_species_codes <- unique(fish_individual_size$species_code)

## Common x-axis limits for comparison
xlim <- range(fish_individual_size$size_mm, na.rm = TRUE)

## Layout settings
par(
  mfrow = c(d, d),
  mar = c(3, 4, 2, 1),
  oma = c(2, 0, 2, 0)
)

for (i in seq_len(d^2)) {
  
  sizes <- fish_individual_size[
    fish_individual_size$species_code == unique_species_codes[i],
    "size_mm"
  ]
  
  hist(
    sizes,
    breaks = 20,
    xlim = xlim,
    col = "lightblue",
    border = "white",
    main = paste("Species:", unique_species_codes[i]),
    xlab = if (i == d) "Size (mm)" else "",
    ylab = "Count"
  )
}

mtext("Fish size distributions by species", outer = TRUE, cex = 1.4)
par(mfrow = c(1, 1))

#
###

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
s <- which(local_foodwebs$operation_id == unique(local_foodwebs$operation_id)[100])
local_foodweb <- unflatten_foodweb(outputs$tab_local_foodwebs_flattened[s,])

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