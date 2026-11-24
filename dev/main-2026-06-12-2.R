##########
## MAIN ##
##########

## Goal: Install fishdatabuilder and foodwebdatabuilder

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

## Load data
data(resource_diet_shift, package = "foodwebbuilder")
data(fish_diet_shift, package = "foodwebbuilder")
data(pred_win, package = "foodwebbuilder")

#
###

####################
## MAIN FISH2SIZE ##
####################

## Format individual fish sizes
fish_individual_size <- fish2size(write_output = F)

## Generate food webs from size for a few selected operations
## For demo purposes
selected_op <- unique(fish_individual_size$operation_id)[1:100]
outputs <- size2webs(num_classes = 3, fish_individual_size, resource_diet_shift, fish_diet_shift, pred_win, write_output=F, selected_op)

## Generate food webs from size for all operations
## Takes a while so not for demo
outputs <- size2webs(num_classes = 9, fish_individual_size, resource_diet_shift, fish_diet_shift, pred_win, write_output=F)

#
###

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
