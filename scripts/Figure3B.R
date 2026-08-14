#!/usr/bin/env Rscript
#" Pie chart from Excel: column 1 = names, column 2 = numbers.

options(echo = TRUE)
Sys.time()

library(openxlsx)
library(ggplot2)

INPUT_XLSX <- "../data/data.Figure3B.xlsx" ## path relative to working directory or absolute path

## ----- Visual tuning -----
PIE_LABEL_SIZE_PT <- 12 ## slice numbers (larger text on the chart)
LEGEND_TEXT_SIZE_PT <- 20 ## legend wording (small)
LEGEND_KEY_MM <- 6 ## coloured square (~height/width of swatch); lower = tighter legend

## ----- Colours -----
## Names MUST match Excel column‑1 slice names exactly after trim().
## Omit names() or leave NULL only if Excel labels differ from pathways below —
## then unnamed colours recycle in plotted order instead.

select_pathway <- c(
  "Biological oxidations",
  "Metabolism of lipids",
  "Metabolism of vitamins and cofactors",
  "Metabolism of porphyrins",
  "Metabolism of amino acids and derivatives",
  "Metabolism of carbohydrates",
  "Metabolism of nucleotides",
  "TCA cycle and respiratory electron transport",
  "Transport of small molecules",
  "Integration of energy metabolism"
)

mycolors <- c(
  "#C53B32", # red
  "#E06B3C", # orange
  "#ECB14A", # yellow
  "#CCB871", # dark yellow green
  "#72BA76", # green
  "#78ADAA", # blue green
  "#5CBEE6", # light blue
  "#3D3D77", # dark purple
  "#992E75", # purple
  "#F5D6A7" # pale yellow
)

COLORS_CUSTOM <- mycolors
names(COLORS_CUSTOM) <- select_pathway

## Uncomment to revert to ggplot2 default hues instead:
## COLORS_CUSTOM <- NULL

## ggplot2 text size is in mm; ~pt / (72.27/72)
GGPLOT_PT <- 72.27 / 72
size_from_pt <- function(pt) pt / GGPLOT_PT

make_fill_scale <- function(level_names, COLORS_CUSTOM) {
  n <- length(level_names)
  if (!n) stop("No factors to colour.")
  
  hues <- scales::hue_pal()(n)
  vals <- stats::setNames(hues, level_names)
  
  if (!is.null(COLORS_CUSTOM) && length(COLORS_CUSTOM) > 0L) {
    if (is.null(names(COLORS_CUSTOM))) {
      vals <- stats::setNames(
        rep_len(as.character(unname(COLORS_CUSTOM)), n),
        level_names
      )
    } else {
      hit <- intersect(level_names, names(COLORS_CUSTOM))
      vals[hit] <- as.character(unname(COLORS_CUSTOM[hit]))
    }
  }
  
  ggplot2::scale_fill_manual(
    breaks = level_names,
    limits = level_names,
    drop = FALSE,
    values = vals[level_names]
  )
}

#" ----------
#" 
#" # Plot UC

prefix <- "UC"
out_file <- paste0("Figure3B.", prefix, ".png")

input_file <- INPUT_XLSX
if (!file.exists(input_file)) stop("File not found: ", input_file)

raw <- openxlsx::read.xlsx(input_file, sheet = prefix, colNames = TRUE)
if (ncol(raw) < 2L) stop("Excel file needs at least 2 columns.")

d <- data.frame(
  name = raw[[1]],
  value = suppressWarnings(as.numeric(raw[[2]])),
  stringsAsFactors = FALSE
)

d$name <- gsub("The citric acid \\(TCA\\) cycle", "TCA cycle", d$name)

d <- d[!is.na(d$name) & !is.na(d$value), , drop = FALSE]
d$name <- trimws(as.character(d$name))

if (!nrow(d)) stop("No rows left after dropping missing values.")
if (any(d$value < 0, na.rm = TRUE)) warning("Negative values present; pie chart sums may look odd.")

# Drop zero slices (optional — remove next line to include zeros)
d <- d[d$value != 0, , drop = FALSE]

if (!nrow(d)) stop("No non-zero numeric values to plot.")

d$name <- factor(d$name, levels = names(COLORS_CUSTOM))

## legend wraps across columns when there are many categories
leg_ncol <- 2

p <- ggplot2::ggplot(d, ggplot2::aes(x = "", y = value, fill = name)) +
  ggplot2::geom_col(aes(colour = name), width = 5) +
  ggplot2::coord_polar(theta = "y") +
  make_fill_scale(levels(d$name), alpha(COLORS_CUSTOM, 0.7)) +
  scale_colour_manual(breaks = levels(d$name), values = alpha(COLORS_CUSTOM, 1), guide = "none") +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.spacing.x = ggplot2::unit(12, "pt"),
    legend.spacing.y = ggplot2::unit(4, "pt"),
    legend.key.size = ggplot2::unit(LEGEND_KEY_MM, "mm"),
    legend.text = ggplot2::element_text(size = LEGEND_TEXT_SIZE_PT),
    legend.margin = ggplot2::margin(t = 2, r = 0, b = 0, l = 0),
    legend.box.margin = ggplot2::margin(t = -4, r = 0, b = 0, l = 0),
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      direction = "horizontal",
      ncol = leg_ncol,
      byrow = TRUE,
      reverse = TRUE
    )
  ) +
  ggplot2::labs(fill = NULL, title = NULL) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        abs(value - round(value)) < 1e-9,
        scales::comma(round(value)),
        scales::comma(value, accuracy = 0.01)
      )
    ),
    position = ggplot2::position_stack(vjust = 0.5),
    size = size_from_pt(PIE_LABEL_SIZE_PT)
  )

ggplot2::ggsave(out_file, plot = p, width = 14, height = 16, dpi = 300)
message("Saved: ", normalizePath(out_file, mustWork = FALSE))

#" ----------
#" 
#" # Plot CD

prefix <- "CD"
out_file <- paste0("Figure3B.", prefix, ".png")

input_file <- INPUT_XLSX
if (!file.exists(input_file)) stop("File not found: ", input_file)

raw <- openxlsx::read.xlsx(input_file, sheet = prefix, colNames = TRUE)
if (ncol(raw) < 2L) stop("Excel file needs at least 2 columns.")

d <- data.frame(
  name = raw[[1]],
  value = suppressWarnings(as.numeric(raw[[2]])),
  stringsAsFactors = FALSE
)

d$name <- gsub("The citric acid \\(TCA\\) cycle", "TCA cycle", d$name)

d <- d[!is.na(d$name) & !is.na(d$value), , drop = FALSE]
d$name <- trimws(as.character(d$name))

if (!nrow(d)) stop("No rows left after dropping missing values.")
if (any(d$value < 0, na.rm = TRUE)) warning("Negative values present; pie chart sums may look odd.")

# Drop zero slices (optional — remove next line to include zeros)
d <- d[d$value != 0, , drop = FALSE]

if (!nrow(d)) stop("No non-zero numeric values to plot.")

d$name <- factor(d$name, levels = names(COLORS_CUSTOM))

## legend wraps across columns when there are many categories
leg_ncol <- 2

p <- ggplot2::ggplot(d, ggplot2::aes(x = "", y = value, fill = name)) +
  ggplot2::geom_col(aes(colour = name), width = 5) +
  ggplot2::coord_polar(theta = "y") +
  make_fill_scale(levels(d$name), alpha(COLORS_CUSTOM, 0.7)) +
  scale_colour_manual(breaks = levels(d$name), values = alpha(COLORS_CUSTOM, 1), guide = "none") +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.spacing.x = ggplot2::unit(12, "pt"),
    legend.spacing.y = ggplot2::unit(4, "pt"),
    legend.key.size = ggplot2::unit(LEGEND_KEY_MM, "mm"),
    legend.text = ggplot2::element_text(size = LEGEND_TEXT_SIZE_PT),
    legend.margin = ggplot2::margin(t = 2, r = 0, b = 0, l = 0),
    legend.box.margin = ggplot2::margin(t = -4, r = 0, b = 0, l = 0),
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      direction = "horizontal",
      ncol = leg_ncol,
      byrow = TRUE,
      reverse = TRUE
    )
  ) +
  ggplot2::labs(fill = NULL, title = NULL) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        abs(value - round(value)) < 1e-9,
        scales::comma(round(value)),
        scales::comma(value, accuracy = 0.01)
      )
    ),
    position = ggplot2::position_stack(vjust = 0.5),
    size = size_from_pt(PIE_LABEL_SIZE_PT)
  )

ggplot2::ggsave(out_file, plot = p, width = 14, height = 16, dpi = 300)
message("Saved: ", normalizePath(out_file, mustWork = FALSE))

sessionInfo()

#" ----------
#"
#" # Done.
