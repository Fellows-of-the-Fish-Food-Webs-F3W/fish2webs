# fish2webs

`fish2webs` is an R package for constructing size-structured aquatic food webs from fish survey data. The package provides a workflow to:

1. Extract and harmonise individual fish size measurements from ASPE surveys.
2. Build size-structured trophic metawebs based on species identity, ontogenetic diet shifts, and predator-prey size constraints.
3. Generate local food webs for individual sampling operations.
4. Compute structural and trophic network metrics for both global and local food webs.

The package builds on the companion packages:

- **fishdatabuilder**: tools for extracting and cleaning fish survey data.
- **foodwebbuilder**: tools for constructing and analysing food webs.

---

# Installation

Install the package directly from GitHub using `pak`:

```r
install.packages("pak")
pak::pak("Fellows-of-the-Fish-Food-Webs-F3W/fish2webs")

library(fish2webs)
library(dplyr) # This is the only way to make the fishdatabuilder package work...
```

The required companion packages are installed automatically.

---

# Quick start

## Load datasets

```r
data(resource_diet_shift, package = "foodwebbuilder")
data(fish_diet_shift, package = "foodwebbuilder")
data(pred_win, package = "foodwebbuilder")
```

## Generate individual fish size measurements

The first step is to extract, clean, and harmonise fish size measurements from the ASPE database.

```r
fish_individual_size <- fish2size(write_output = FALSE)
```

The resulting table contains validated individual fish size measurements suitable for downstream ecological analyses.

---

## Build food webs

### Example using a subset of operations

For demonstration purposes, food webs can be generated for a limited number of sampling operations.

```r
selected_operations <- unique(
  fish_individual_size$operation_id
)[1:100]

outputs <- size2webs(
  num_classes = 3,
  ind_measure = fish_individual_size,
  resource_diet_shift = resource_diet_shift,
  fish_diet_shift = fish_diet_shift,
  pred_win = pred_win,
  write_output = FALSE,
  selected_operations = selected_operations
)
```

### Full dataset

Generating food webs for all operations may require substantial computation time.

```r
outputs <- size2webs(
  num_classes = 9,
  ind_measure = fish_individual_size,
  resource_diet_shift = resource_diet_shift,
  fish_diet_shift = fish_diet_shift,
  pred_win = pred_win,
  write_output = FALSE
)
```

---

# Output objects

`size2webs()` returns a list containing:

| Object | Description |
|----------|----------|
| `tab_size_classes` | Species-specific size class boundaries |
| `tab_metaweb_flattened` | Flattened representation of the global metaweb |
| `tab_local_foodwebs_flattened` | Flattened representation of local food webs |
| `tab_local_foodwebs_summary_metrics` | Summary metrics for each local food web |

---

# Example: fish size distributions

Visualise size distributions for selected species.

```r
d <- 4

unique_species_codes <- unique(
  fish_individual_size$species_code
)

xlim <- range(
  fish_individual_size$size_mm,
  na.rm = TRUE
)

par(
  mfrow = c(d, d),
  mar = c(3, 4, 2, 1),
  oma = c(2, 0, 2, 0)
)

for (i in seq_len(d^2)) {

  sizes <- fish_individual_size[
    fish_individual_size$species_code ==
      unique_species_codes[i],
    "size_mm"
  ]

  hist(
    sizes,
    breaks = 20,
    xlim = xlim,
    col = "lightblue",
    border = "white",
    main = paste(
      "Species:",
      unique_species_codes[i]
    ),
    xlab = "",
    ylab = "Count"
  )
}

mtext(
  "Fish size distributions by species",
  outer = TRUE,
  cex = 1.4
)

par(mfrow = c(1, 1))
```

---

# Example: analyse the global metaweb

Reconstruct the metaweb and compute structural metrics.

```r
metaweb <- unflatten_foodweb(
  outputs$tab_metaweb_flattened
)

basal_nodes <- get_basal_nodes(metaweb)
leaf_nodes  <- get_leaf_nodes(metaweb)

degree_in  <- compute_inward_degree(metaweb)
degree_out <- compute_outward_degree(metaweb)

TL <- compute_trophic_level(metaweb)

TB <- compute_trophic_breadth(
  metaweb,
  TL
)

fluxes <- compute_bottom_up_fluxes(metaweb)
```

Visualise trophic structure:

```r
plot_network(
  metaweb,
  x = TB,
  y = TL,
  labels = colnames(metaweb),
  xlab = "Trophic breadth",
  ylab = "Trophic level"
)
```

---

# Example: analyse a local food web

Select one operation and reconstruct its local food web.

```r
local_foodwebs <-
  outputs$tab_local_foodwebs_flattened

s <- which(
  local_foodwebs$operation_id ==
    unique(local_foodwebs$operation_id)[100]
)

local_foodweb <- unflatten_foodweb(
  local_foodwebs[s, ]
)
```

Compute structural metrics:

```r
basal_nodes <- get_basal_nodes(local_foodweb)
leaf_nodes  <- get_leaf_nodes(local_foodweb)

degree_in  <- compute_inward_degree(local_foodweb)
degree_out <- compute_outward_degree(local_foodweb)

TL <- compute_trophic_level(local_foodweb)

TB <- compute_trophic_breadth(
  local_foodweb,
  TL
)

fluxes <- compute_bottom_up_fluxes(local_foodweb)
```

Visualise the local food web:

```r
plot_network(
  local_foodweb,
  x = TB,
  y = TL,
  labels = colnames(local_foodweb),
  xlab = "Trophic breadth",
  ylab = "Trophic level"
)
```

---

# Workflow summary

```text
ASPE surveys
      │
      ▼
fish2size()
      │
      ▼
Individual fish sizes
      │
      ▼
size2webs()
      │
      ├── Size classes
      ├── Global metaweb
      ├── Local food webs
      └── Food-web metrics
```

---

# Citation

If you use `fish2webs` in published work, please cite the package and the associated methodological publication(s).
