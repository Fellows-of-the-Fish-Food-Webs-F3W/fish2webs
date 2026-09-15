######################
## EXAMPLE ANALYSIS ##
######################

## Goal: 

######################
## GLOBAL VARIABLES ##
######################

## User defined variables
INSTALL <- FALSE # {TRUE, FALSE} to install all packages in dependencies before running
DEBUG_MODE <- TRUE # {TRUE, FALSE} to run whole pipeline on smaller operation set
OUTPUT_FOLDER <- "./outputs" # Folder where the output tables will be stored
system("mkdir ./outputs") # TODO: implement saving files in output folder

#
###

###########
## SETUP ##
###########

## Install library
if (INSTALL == TRUE){
  install.packages("pak")
  pak::pak("Fellows-of-the-Fish-Food-Webs-F3W/fish2webs")
}

## Load library
library(fish2webs)
par(family = "sans")

## Load data
data(resource_diet_shift, package = "foodwebbuilder")
data(fish_diet_shift, package = "foodwebbuilder")
data(pred_win, package = "foodwebbuilder")

#
###

################
## LOAD FILES ##
################

## Move to output dir
setwd(OUTPUT_FOLDER)

tab_site_information <- read.csv("1_tab_site_information.csv")
tab_operation_information <- read.csv("2_tab_operation_information.csv")
tab_species_information <- read.csv("3_tab_species_information.csv")
tab_fish_individual_size_weight <- read.csv("4_tab_fish_individual_size_weight.csv")
tab_community_metrics <- read.csv("5_tab_community_metrics.csv")
tab_species_level_metrics <- read.csv("6_tab_species_level_metrics.csv")
tab_trophic_species_size_classes <- read.csv("9_tab_trophic_species_size_classes.csv")
tab_metaweb <- read.csv("10_tab_metaweb.csv")
tab_local_foodwebs <- read.csv("11_tab_local_foodwebs.csv")
tab_local_foodweb_metrics <- read.csv("12_tab_local_foodweb_metrics.csv")

#
###

#########################################
## EXAMPLE ANALYSIS SIZE DISTRIBUTIONS ##
#########################################

## Number of panels in figure rows and columns
d <- 4

## Get species codes
unique_species_codes <- unique(tab_fish_individual_size_weight$species_code)

## Common x-axis limits for comparison
xlim <- range(tab_fish_individual_size_weight$size_mm, na.rm = TRUE)

## Layout settings
par(
  mfrow = c(d, d),
  mar = c(3, 4, 2, 1),
  oma = c(2, 0, 2, 0)
)

for (i in seq_len(d^2)) {
  
  sizes <- tab_fish_individual_size_weight[
    tab_fish_individual_size_weight$species_code == unique_species_codes[i],
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

library(foodwebbuilder)

## Rebuild metaweb
metaweb <- unflatten_foodweb(tab_metaweb)
dim(metaweb)

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

## plot TL vs TB 
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

## Plot radial 
plot_network_radial(
  fluxes,
  x = TL,
  y = TL,
  labels = colnames(metaweb),
  rotations = 1.5,
  scale = 1.25,
  line_width_max = 1
)

## Plot radial 2
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
local_foodwebs <- tab_local_foodwebs
s <- which(local_foodwebs$operation_id == unique(local_foodwebs$operation_id)[100])
local_foodweb <- unflatten_foodweb(local_foodwebs[s,])

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

######################
## HELPER FUNCTIONS ##
######################

## Efficient acquisition of times
get_op_time <- function(op_id, op_ids_and_times=tab_operation_information){
  s <- which(op_ids_and_times$operation_id == op_id)
  return(op_ids_and_times$date[s])
}
head(tab_operation_information)
get_op_time(39)

## TODO: update code of function in foodwebbuilder repository
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
                        1.5), cex = 0, bty = "n", xlab = xlab, ylab = ylab, 
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

######################################################
## EXAMPLE ANALYSIS: SINGLE-SITE COMMUNITY PATTERNS ##
######################################################

## Select site and obtain all operations
site_id_ <- unique(tab_operation_information$site_id)[310]
s <- which(tab_operation_information$site_id == site_id_)
operations_ <- tab_operation_information$operation_id[s]
print(site_id_)

## Collect individual fish measurements
tab_fish_individual_size_weight_single_site <- NULL
for (operation in operations_){
  s <- which(tab_fish_individual_size_weight$operation_id == operation)
  tab_fish_individual_size_weight_single_site_ <- tab_fish_individual_size_weight[s,]
  tab_fish_individual_size_weight_single_site <- rbind(
    tab_fish_individual_size_weight_single_site, 
    tab_fish_individual_size_weight_single_site_
    )
}

## Format table
tab_fish_individual_size_weight_single_site <- data.frame(tab_fish_individual_size_weight_single_site)
colnames(tab_fish_individual_size_weight_single_site) <- colnames(tab_fish_individual_size_weight) 

## Checks
print(site_id_)
head(tab_fish_individual_size_weight_single_site)
nrow(tab_fish_individual_size_weight_single_site)

## Shorten name
site_ind_measure <- tab_fish_individual_size_weight_single_site

## Figure: Single site species counts
pdf(file="plot-single-site-counts.pdf", width=2.5, height=4)
par(mar = c(5, 5, 2, 2), xpd = TRUE, mfrow=c(1,1), bty="l", family="mono")
#
## Counts
res <- table(site_ind_measure$species_code)
colors <- rainbow(length(res), start = 0.0, end = 0.8, alpha=1.0)
bar_positions <- barplot(res, horiz=T, las=1, col=colors, xlab="Total count", cex.names=0.7)
# axis(1, family="sans")
# axis(2, at=bar_positions, labels=names(res), las=1, family="mono")
#
dev.off()

## Figure: Single site species biomass
pdf(file="plot-single-site-biomasses.pdf", width=2.5, height=4)
par(mar = c(5, 5, 2, 2), xpd = TRUE, mfrow=c(1,1), bty="l", family="mono")
#
## Biomass
res <- NULL
site_unique_species_codes <- sort(unique(site_ind_measure$species_code))
for (species_code_ in site_unique_species_codes){
  tab <- site_ind_measure[which(site_ind_measure$species_code == species_code_),]
  res <- c(res, sum(tab$weight_g))
}
names(res) <- site_unique_species_codes
#
## Plot
colors <- rainbow(length(res), start = 0.0, end = 0.8, alpha=1.0)
bar_positions <- barplot(res, horiz=T, las=1, col=colors, xlab="Total biomass (g)", cex.names=0.7)
# axis(1, family="sans")
# axis(2, at=bar_positions, labels=names(res), las=1, family="mono")
#
dev.off()

## Figure: Single site species size distributions
#
## Three most abundant species
res <- table(site_ind_measure$species_code)
selected_species_codes <- names(sort(res, decreasing = TRUE))[1:3]
#
k <- 1
for (species_code_ in selected_species_codes){
  #
  ## Size distribution
  pdf(file=paste("plot-single-site-sizes-", k, ".pdf", sep=""), width=2.5, height=4)
  par(mar = c(5, 5, 2, 2), xpd = TRUE, mfrow=c(1,1), bty="l", family="mono")
  tab <- site_ind_measure[which(site_ind_measure$species_code == species_code_),]
  hist(tab$size_mm, las=1, xlab="Body size (mm)", main=paste(species_code_))
  dev.off()
  k <- k + 1
}
#

## Selected species counts
selected_species <- c("VAN", "GOU", "ABL")
#
## Subset
s <- NULL
for (selected_species_ in selected_species){
  s <- c(s, which(site_ind_measure$species_code ==  selected_species_))
}
subset_site_ind_measure <- site_ind_measure[s,]
#
## Counts
res <- table(subset_site_ind_measure$species_code)
#
## Plot
colors <- rainbow(length(res))
barplot(res, horiz=T, las=1, col=colors)

## Species abundances through operations
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
#
## Convert to date to continuous time
times <- NULL 
for (unique_op_ids_ in unique_op_ids){
  times <- c(times, get_op_time(unique_op_ids_))
}
times <- as.Date(times)
times <- as.numeric(times - min(times)) # Time in days
s <- order(times)
x <- times[s]/365
y <- tab_species_abundances[s,]

## Figure: Time series of proportional abundances
pdf(file="plot-single-site-proportional-abundances.pdf", width=5, height=4)
par(mar = c(5, 5, 2, 2), xpd = TRUE, mfrow=c(1,1), bty="l", family="mono")
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
cols <- rainbow(length(sp), start = 0.0, end = 0.8, alpha=1.0)
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

## Figure: Time series of local food webs
k <- 1
order_ <- order(times)
for (k in 1:length(unique_op_ids)){
  
  pdf(file=paste("plot-single-site-network-",k,".pdf",sep=""), width=5, height=4)
  par(mar = c(5, 5, 2, 2), family="mono")
  
  ## Check time
  unique_op_ids_ <- unique_op_ids[order_[k]]
  print(times[order_[k]]/365)
  
  ## Select a local food web
  local_foodwebs <- tab_local_foodwebs
  s <- which(local_foodwebs$operation_id == unique_op_ids_)
  local_foodweb <- unflatten_foodweb(local_foodwebs[s,])
  
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
  
  ## Visualise network
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
  legend("bottomright", legend = c(
    paste("t =", round(times[order_[k]]/365), "years"),
    paste("mean TL =", round(mean(TL), 2)),
    paste("max TL =", round(max(TL), 2))
  ), bty="n")
  # legend(-0.55, -0.25, legend = paste("t =",round(times[order_[k]]/365),"years"), bty="n", cex=2)
  
  ##
  dev.off()
  
}

#
###

########################
## SAMPLES MULTI-SITE ##
########################

## Figure
pdf(file="plot-multi-site-num-samples.pdf", width=8, height=8)

## Graphical parameters
par(mar = c(4, 4, 4, 4), xpd = TRUE, mfrow=c(1,1), family="mono")

## Keep only one observation per operation
s <- match(unique(tab_fish_individual_size_weight$operation_id), tab_fish_individual_size_weight$operation_id) 
tab <- tab_fish_individual_size_weight[s,]
nrow(tab)
length(unique(tab_fish_individual_size_weight$operation_id))

## Get site identity
s <- match(tab$operation_id, tab_operation_information$operation_id)
site_ids <- tab_operation_information$site_id[s]
head(site_ids)
length(site_ids)
length(unique(site_ids))

## Add column
tab$site_id <- site_ids

## Get operation count per site
res <- table(tab$site_id)
res_range <- range(res)
length(names(res))
head(names(res))

## Add column
s <- match(tab$site_id, names(res))
tab$res <- res[s]
tab$res_site_check <- names(res)[s]
head(tab)

## Get site coordinates
s <- match(tab$site_id, tab_site_information$site_id)
site_coordinates <- tab_site_information[s,]

## Add columns
tab$site_x <- site_coordinates$x
tab$site_y <- site_coordinates$y

## Filter coordinates
s <- which(tab$site_y > 0) # TODO: Fix bug with coordinate
tab <- tab[s,] #site_coordinates <- site_coordinates[s,]

## Plot
u <- tab$res / max(tab$res) # u <- (res-min(res))/(max(res)-min(res))
plot(tab$site_x, tab$site_y, col="black", 
     pch=16, cex=2*u^0.25, 
     xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
#
## Visualise specific site
s <- which(tab$site_id == 8092)
points(tab$site_x[s], tab$site_y[s], col="black", cex=3, pch=1, lwd=8)
points(tab$site_x[s], tab$site_y[s], col="cyan", cex=3, pch=1, lwd=3)
#
legend_labels <- round(c(min(u), 0.1, 0.5, 1) * max(res))
legend_prop <- c(min(u), 0.1, 0.5, 1)
legend("topleft",
       inset = c(-0.0, 0.0),
       legend = c("N sampling events", legend_labels),
       # legend = c("N sampling events", paste0(100 * legend_prop, "%")),
       col = c(rep("black", length(legend_prop))),
       pch = c(16, rep(16, length(legend_prop))),
       pt.cex = c(0, 2*legend_prop^0.25),
       bty = "n",
       xpd = TRUE)

## End 
dev.off()

#
###

########################
## SPECIES MULTI-SITE ##
########################

## User-defined entries
selected_species <- c("TRF", "BRO", "SIL")

## Figure
pdf(file="plot-multi-site-species-distributions.pdf", width=8, height=8)
#
## Graphical parameters
par(mar = c(4, 4, 4, 4), xpd = TRUE, family="mono")
colors <- rainbow(length(selected_species))
#
## For each species
k <- 1
res_range <- list()
for (selected_species_ in selected_species){
# selected_species_ <- selected_species[1]
  
  ## Subset only selected species
  s <- which(tab_fish_individual_size_weight$species_code == selected_species_)
  species_ind_measure <- tab_fish_individual_size_weight[s,]
  nrow(species_ind_measure)
  
  ## Get site identity
  s <- match(species_ind_measure$operation_id, tab_operation_information$operation_id)
  site_ids <- tab_operation_information$site_id[s]
  head(site_ids)
  length(site_ids)
  length(unique(site_ids))
  
  ## Get coordinates
  s <- match(site_ids, tab_site_information$site_id)
  site_coordinates <- tab_site_information[s,]
  head(site_coordinates)
  tail(site_coordinates)

  ## Get fish count per site
  res <- table(site_ids)
  res_range[[selected_species_]] <- range(res)
  length(res)
  head(res)
  head(site_ids)
  
  ## Match
  s <- match(site_ids, names(res))
  res <- res[s]
  head(res)
  head(site_ids)
  
  ## Plot
  u <- res / max(res) # u <- (res-min(res))/(max(res)-min(res))
  if (k == 1){
    plot(site_coordinates$x, site_coordinates$y, col=colors[k], 
         pch=16, cex=2*u^0.25, 
         xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  } else {
    points(site_coordinates$x, site_coordinates$y, col=colors[k], 
           pch=16, cex=2*u^0.25)
  }
  #
  k <- k + 1
}
#
legend_prop <- c(0.01, 0.1, 0.5, 1)
legend("bottomleft",
       inset = c(-0.0, 0.0),
       legend = c("Species", 
                  selected_species,
                  "", 
                  "% max abundance",
                  paste0(100 * legend_prop, "%")),
       col = c(NA, 
               colors, 
               NA, 
               NA, 
               rep("black", length(legend_prop))),
       pch = c(NA, 
               rep(16, length(selected_species)),
               NA, 
               NA, 
               rep(16, length(legend_prop))),
       pt.cex = c(NA, 
                  rep(1.5, length(selected_species)),
                  NA, 
                  NA, 
                  2 * legend_prop^0.25),
       bty = "n",
       xpd = TRUE)
#
## End
dev.off()

#
###
