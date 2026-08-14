library(ggplot2)

cell_colors <- c(
  "B.memory"="#8B0000", "B.naive"="#FF6347", "Plasma.cells"="#FF69B4",
  "T.CD4.naive"="#6A1B9A", "T.CD4.memory"="#AB47BC", "T.regs"="#CE93D8",
  "T.CD8.naive"="#0A3A85", "T.CD8.em"="#1565C0", "T.CD8.cm"="#42A5F5",
  "NK"="#311B92",
  "Macrophages"="#BF360C", "DC"="#B8960F", "Mon.Classical"="#FF6F00", "Mon.NonClassical"="#FFAB91",
  "Neutrophils"="#E8D44D",
  "Epithelial"="#1B5E20", "Fibroblasts"="#43A047", "Endothelial"="#76C87E", "NonImmune"="#4E342E",
  "EMT"="#808080"
)


stagger_labels <- function(centroids, y_step = NULL, x_band = NULL, y_band = NULL) {
  span <- max(
    diff(range(centroids$coord_X, na.rm = TRUE)),
    diff(range(centroids$coord_Y, na.rm = TRUE)),
    1
  )
  if (is.null(y_step)) y_step <- span * 0.04
  if (is.null(x_band)) x_band <- span * 0.12
  if (is.null(y_band)) y_band <- span * 0.085

  label_X <- centroids$coord_X
  label_Y <- centroids$coord_Y
  n <- nrow(centroids)
  parent <- seq_len(n)

  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      if (abs(centroids$coord_X[i] - centroids$coord_X[j]) <= x_band &&
          abs(centroids$coord_Y[i] - centroids$coord_Y[j]) <= y_band) {
        ri <- i
        while (parent[ri] != ri) ri <- parent[ri]
        rj <- j
        while (parent[rj] != rj) rj <- parent[rj]
        if (ri != rj) parent[rj] <- ri
      }
    }
  }

  roots <- vapply(seq_len(n), function(i) {
    r <- i
    while (parent[r] != r) r <- parent[r]
    r
  }, integer(1))

  for (cl in unique(roots)) {
    idx <- which(roots == cl)
    if (length(idx) == 1L) next
    idx <- idx[order(centroids$coord_Y[idx])]
    k <- length(idx)
    offsets <- (seq_len(k) - (k + 1) / 2) * y_step
    label_Y[idx] <- centroids$coord_Y[idx] + offsets
  }

  centroids$label_X <- label_X
  centroids$label_Y <- label_Y
  centroids
}

make_plot <- function(file, title, zoom = TRUE, patient_col = NULL, label_nudges = NULL) {
  d <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)

  n_cells <- nrow(d)
  if (!is.null(patient_col) && patient_col %in% names(d)) {
    n_patients <- length(unique(d[[patient_col]]))
  } else {
    n_patients <- length(unique(d$sample))
  }

  set.seed(42)
  d <- d[sample(nrow(d)), ]

  centroids <- aggregate(cbind(coord_X, coord_Y) ~ cellstates, data = d, FUN = median)

  cell_counts <- table(d$cellstates)
  centroids$label <- paste0(
    centroids$cellstates, " (",
    trimws(format(as.integer(cell_counts[centroids$cellstates]), big.mark = ",")), ")"
  )

  xlim <- NULL
  ylim <- NULL
  if (zoom) {
    xr <- quantile(d$coord_X, c(0.001, 0.999))
    yr <- quantile(d$coord_Y, c(0.001, 0.999))
    pad <- 300
    xlim <- c(xr[1] - pad, xr[2] + pad)
    ylim <- c(yr[1] - pad, yr[2] + pad)
  }

  centroids <- stagger_labels(centroids)

  if (!is.null(label_nudges)) {
    for (nm in names(label_nudges)) {
      idx <- match(nm, centroids$cellstates)
      if (!is.na(idx)) {
        centroids$label_Y[idx] <- centroids$label_Y[idx] + label_nudges[[nm]]
      }
    }
  }

  p <- ggplot(d, aes(x = coord_X, y = coord_Y, color = cellstates)) +
    geom_point(size = 0.8, alpha = 0.7, stroke = 0, shape = 16) +
    geom_segment(
      data = centroids,
      aes(x = coord_X, y = coord_Y, xend = label_X, yend = label_Y, color = cellstates),
      linewidth = 0.25,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_label(
      data = centroids,
      aes(x = label_X, y = label_Y, label = label, color = cellstates),
      size = 4.5,
      fontface = "bold",
      fill = alpha("white", 0.85),
      label.padding = unit(0.1, "lines"),
      linewidth = 0,
      show.legend = FALSE
    ) +
    scale_color_manual(values = cell_colors, drop = TRUE) +
    theme_classic(base_size = 18) +
    theme(
      plot.title     = element_text(hjust = 0.5, face = "bold", size = 40),
      axis.title     = element_text(size = 20, face = "bold"),
      axis.text      = element_blank(),
      axis.ticks     = element_blank(),
      axis.line      = element_line(linewidth = 0.6),
      legend.position = "none",
      plot.margin    = margin(15, 15, 15, 15)
    ) +
    labs(title = title, x = "SPRING-1", y = "SPRING-2")

  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }

  annotation <- paste0("Cells = ", format(n_cells, big.mark = ","),
                        "\nPatients = ", n_patients)
  p <- p + annotate("text", x = Inf, y = -Inf, label = annotation,
                     hjust = 1.1, vjust = -0.5, size = 10, fontface = "italic")

  return(p)
}

# --- CD Ileum ---
cat("Plotting CD Ileum...\n")
p1 <- make_plot(
  "../data/metadata_coords_ileum.tsv",
  "CD Ileum",
  zoom = TRUE,
  patient_col = "sample",
  label_nudges = c("NK" = 550)
)
ggsave("cellstates_ileum.pdf", p1, width = 12, height = 9, device = "pdf", bg = "white")
cat("  Saved cellstates_ileum.pdf\n")

# --- UC Colon ---
cat("Plotting UC Colon...\n")
p2 <- make_plot(
  "../data/metadata_coords_colon.tsv",
  "UC Colon",
  zoom = FALSE,
  patient_col = "subject"
)
ggsave("cellstates_colon.pdf", p2, width = 12, height = 9, device = "pdf", bg = "white")
cat("  Saved cellstates_colon.pdf\n")

cat("Done.\n")
