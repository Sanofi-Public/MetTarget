# Network figure from TSV: edges + DE node stats (ggplot2).
#
# Install:
#   install.packages("ggplot2")

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Missing package ggplot2. Install: install.packages(\"ggplot2\")", call. = FALSE)
  }
})

library(ggplot2)

pick_col <- function(df, synonyms) {
  cn <- colnames(df)
  canon <- function(x) tolower(trimws(gsub("\\s+", " ", x)))
  for (syn in synonyms) {
    w <- cn[canon(cn) == canon(syn)]
    if (length(w)) return(w[[1]])
  }
  NA_character_
}

read_tsv <- function(path) {
  df <- utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "\"",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    na.strings = c("", "NA", "NaN", "NULL"),
    strip.white = TRUE
  )
  if (!ncol(df) || ncol(df) == 1L) {
    stop(
      "Read 0 usable columns — is the file comma-separated?\n ",
      path,
      call. = FALSE
    )
  }
  df
}

read_edges_tsv <- function(path) {
  d <- read_tsv(path)
  from_c <- pick_col(d, c("from", "source", "tail"))
  to_c <- pick_col(d, c("to", "target", "head"))
  if (is.na(from_c) || is.na(to_c)) {
    stop("Edge TSV needs recognizable 'from' and 'to' columns: ", path)
  }
  type_c <- pick_col(d, "type")
  colour_c <- pick_col(d, c("color", "colour", "edge_color", "edge_colour"))
  line_c <- pick_col(d, c("linetype", "line_type", "linestyle", "line_style"))
  out <- data.frame(
    from = as.character(d[[from_c]]),
    to = as.character(d[[to_c]]),
    stringsAsFactors = FALSE
  )
  out$type <- if (!is.na(type_c)) as.character(d[[type_c]]) else rep(NA_character_, nrow(out))
  out$linestyle <- if (!is.na(line_c)) as.character(d[[line_c]]) else rep(NA_character_, nrow(out))
  out$colour <- if (!is.na(colour_c)) as.character(d[[colour_c]]) else rep("grey", nrow(out))
  out <- out[stats::complete.cases(out[, c("from", "to"), drop = FALSE]), , drop = FALSE]
  out <- out[!is.na(out$from) & !is.na(out$to), , drop = FALSE]
  out <- out[nzchar(out$from) & nzchar(out$to), , drop = FALSE]
  if (!nrow(out)) stop("No valid edges found in ", path)
  out
}

read_nodes_tsv <- function(path) {
  d <- read_tsv(path)
  name_c <- pick_col(d, c("name", "gene", "id", "symbol"))
  lfc_c <- pick_col(d, c("log2FC", "log2fc", "log2_fc", "logFC", "logfc"))
  p_c <- pick_col(d, c("p-value", "p.value", "p_value", "pvalue", "p"))
  nlp_c <- pick_col(d, c("neg.log10p", "neg_log10p", "neg_log10_p", "-log10p"))
  hub_c <- pick_col(d, c("hub", "is_hub"))
  shape_c <- pick_col(d, c("shape", "node_shape"))
  x_c <- pick_col(d, c("x", "layout_x"))
  y_c <- pick_col(d, c("y", "layout_y"))
  if (is.na(name_c)) stop("Node TSV needs 'name' (or gene/id): ", path)
  nm <- trimws(as.character(d[[name_c]]))
  out <- data.frame(
    name = nm,
    log2FC = if (!is.na(lfc_c)) suppressWarnings(as.numeric(d[[lfc_c]])) else rep(NA_real_, nrow(d)),
    p_value = if (!is.na(p_c)) suppressWarnings(as.numeric(d[[p_c]])) else rep(NA_real_, nrow(d)),
    neg_log10P = if (!is.na(nlp_c)) suppressWarnings(as.numeric(d[[nlp_c]])) else rep(NA_real_, nrow(d)),
    x = if (!is.na(x_c)) suppressWarnings(as.numeric(d[[x_c]])) else rep(NA_real_, nrow(d)),
    y = if (!is.na(y_c)) suppressWarnings(as.numeric(d[[y_c]])) else rep(NA_real_, nrow(d)),
    stringsAsFactors = FALSE
  )
  if (!is.na(hub_c)) {
    hv <- d[[hub_c]]
    hub_vec <- rep(FALSE, nrow(out))
    if (is.logical(hv)) {
      hub_vec <- hv
      if (length(hub_vec) == 1L) hub_vec <- rep(hub_vec, nrow(out))
    } else {
      hv2 <-tolower(trimws(as.character(hv)))
      hub_vec <- hv2 %in% c("true", "yes", "1", "hub", "y") | hv == 1L
      if (length(hub_vec) == 1L) hub_vec <- rep(hub_vec, nrow(out))
    }
    if (length(hub_vec) != nrow(out)) hub_vec <- rep(as.logical(hub_vec[1]), nrow(out))
    out$hub <- hub_vec
  } else {
    out$hub <- rep(FALSE, nrow(out))
  }
  if (!is.na(shape_c)) {
    sch <- trimws(as.character(d[[shape_c]]))
    if (length(sch) != nrow(out)) sch <- rep(sch, length.out = nrow(out))
    sch[is.na(sch) | !nzchar(sch)] <- "circle"
    out$shape <- sch
  } else {
    out$shape <- rep("circle", nrow(out))
  }
  out <- out[!is.na(out$name) & nzchar(out$name), , drop = FALSE]
  if (!nrow(out)) stop("No usable node rows in ", path)
  out
}

add_missing_nodes <- function(nodes, edges) {
  need <- unique(c(edges$from, edges$to))
  miss <- setdiff(need, nodes$name)
  if (!length(miss)) return(nodes)
  pad <- data.frame(
    name = miss,
    log2FC = NA_real_,
    p_value = NA_real_,
    neg_log10P = NA_real_,
    x = NA_real_,
    y = NA_real_,
    hub = FALSE,
    shape = "circle",
    stringsAsFactors = FALSE
  )
  rbind(nodes, pad)
}

layout_norm <- function(xy, pad = 0.08) {
  xr <- range(xy[, 1], na.rm = TRUE)
  yr <- range(xy[, 2], na.rm = TRUE)
  span <- max(diff(xr), diff(yr), 1e-6)
  mx <- mean(xr)
  my <- mean(yr)
  x <- (xy[, 1] - mx) / span * (1 - 2 * pad) + 0.5
  y <- (xy[, 2] - my) / span * (1 - 2 * pad) + 0.5
  cbind(pmin(pmax(x, 0), 1), pmin(pmax(y, 0), 1))
}

fr_layout_from_edges <- function(from, to, nodes, niter = 2000, seed = 42) {
  nodes <- unique(as.character(nodes))
  n <- length(nodes)
  idx <- stats::setNames(seq_len(n), nodes)
  ei <- unname(idx[from])
  ej <- unname(idx[to])

  set.seed(seed)
  pos <- matrix(stats::rnorm(n * 2), ncol = 2)
  k <- sqrt(1 / n)
  area <- 1

  for (iter in seq_len(niter)) {
    disp <- matrix(0, n, 2)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i == j) next
        d <- pos[i, ] - pos[j, ]
        dist <- sqrt(sum(d * d))
        if (dist < 1e-9) dist <- 1e-9
        repulsive <- (k * k / dist) * (d / dist)
        disp[i, ] <- disp[i, ] + repulsive
      }
    }
    for (e in seq_along(ei)) {
      i <- ei[e]
      j <- ej[e]
      d <- pos[j, ] - pos[i, ]
      dist <- sqrt(sum(d * d))
      if (dist < 1e-9) dist <- 1e-9
      attractive <- (dist * dist / k) * (d / dist)
      disp[i, ] <- disp[i, ] + attractive
      disp[j, ] <- disp[j, ] - attractive
    }
    temp <- area * (1 - (iter - 1) / max(niter - 1, 1))
    for (i in seq_len(n)) {
      disp_len <- sqrt(sum(disp[i, ] * disp[i, ]))
      if (disp_len > 1e-9) {
        max_disp <- min(disp_len, temp)
        pos[i, ] <- pos[i, ] + disp[i, ] * max_disp / disp_len
      }
    }
  }

  co <- layout_norm(pos)
  rownames(co) <- nodes
  co
}

assign_node_layout <- function(nodes_tbl, edges_tbl) {
  if (all(is.finite(nodes_tbl$x)) && all(is.finite(nodes_tbl$y))) {
    return(nodes_tbl)
  }
  co <- fr_layout_from_edges(edges_tbl$from, edges_tbl$to, nodes_tbl$name)
  ix <- match(nodes_tbl$name, rownames(co))
  if (any(is.na(ix))) stop("Some node names missing from computed layout.")
  nodes_tbl$x <- co[ix, 1]
  nodes_tbl$y <- co[ix, 2]
  nodes_tbl
}

shrink_segment <- function(x0, y0, x1, y1, r0 = 0.022, r1 = 0.02) {
  dx <- x1 - x0
  dy <- y1 - y0
  len <- sqrt(dx * dx + dy * dy)
  if (!is.finite(len) || len < 1e-9) {
    return(list(x = x0, y = y0, xend = x1, yend = y1))
  }
  ux <- dx / len
  uy <- dy / len
  list(x = x0 + r0 * ux, y = y0 + r0 * uy, xend = x1 - r1 * ux, yend = y1 - r1 * uy)
}

resolve_linetype <- function(types) {
  valid <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  t_raw <- ifelse(is.na(types), NA_character_, as.character(types))
  t <- ifelse(is.na(t_raw), "dashed", tolower(trimws(t_raw)))
  t[!nzchar(t)] <- "dashed"
  unk <- !(t %in% valid)
  if (!any(unk)) return(t)
  uniq <- sort(unique(types[unk]))
  pal <- rep(c("dotdash", "longdash", "twodash"), length.out = max(1L, length(uniq)))
  mp <- stats::setNames(pal[seq_along(uniq)], as.character(uniq))
  t[unk] <- unname(mp[as.character(types[unk])])
  t[!t %in% valid] <- "dashed"
  t
}

# Map biological edge `type` column to arrowhead style targets.
interaction_effect_from_type <- function(t) {
  xl <- ifelse(is.na(t), "", tolower(trimws(as.character(t))))
  inh_hit <- xl %in% c("inhibit", "inhibition", "inh") |
    grepl("^inhibit", xl, perl = TRUE) |
    grepl("repress", xl, perl = TRUE) |
    grepl("^suppress", xl, perl = TRUE) |
    xl %in% c("negative", "down", "downregulate", "down-regulate")
  act_hit <- xl %in% c(
    "activate", "activation", "activating",
    "+", "positive",
    "up", "upregulate", "up-regulate"
  ) |
    grepl("^activ", xl, perl = TRUE) |
    grepl("promot", xl, perl = TRUE) |
    grepl("stimul", xl, perl = TRUE)
  
  conflict <- inh_hit & act_hit
  if (any(conflict)) {
    message("Some edges matched both activate and inhibit; using inhibit for those rows.")
    act_hit <- act_hit & !conflict
  }
  
  eff <- rep(NA_character_, length(xl))
  eff[inh_hit] <- "inhibit"
  eff[!inh_hit & act_hit] <- "activate"
  
  unresolved <- nzchar(xl) & is.na(eff)
  if (any(unresolved)) {
    bad <- paste(unique(utils::head(xl[unresolved], 20)), collapse = ", ")
    if (nzchar(bad)) {
      message(sprintf("Unknown interaction type (%s); using activate arrows.", bad))
    }
  }
  miss <- !nzchar(xl) | is.na(eff) | unresolved
  eff[miss] <- "activate"
  eff
}

# Kinase hubs: fixed hue per project brief.
KINASE_IDS <- toupper(c("SPHK1", "SPHK2"))
KINASE_FILL <- "#C19364"

shape_key <- function(s) {
  k <- trimws(as.character(s))
  k <- tolower(gsub("\\s+", "_", k))
  k[is.na(k) | !nzchar(k) | k == "na"] <- "circle"
  k
}

shape_pch_vectorized <- function(labels) {
  k <- shape_key(labels)
  pch <- rep(21L, length(k))
  circle_like <- k %in% c("circle", "round", "dot", "point", "o")
  pch[circle_like] <- 21L
  pch[k %in% c("square", "rect", "rectangle", "box")] <- 22L
  pch[k %in% c("diamond", "diam")] <- 23L
  pch[k %in% c("triangle", "triangle_up", "tri", "triangleup", "^")] <- 24L
  pch[k %in% c(
    "triangle_down", "triangle_dn", "tri_dn", "tri_down",
    "invtriangle", "inv_triangle", "v"
  )] <- 25L
  pch
}

build_shape_scale <- function(shape_labels) {
  lev <- unique(as.character(shape_labels))
  canon <- shape_key(lev)
  known <- canon %in% c(
    "circle", "round", "dot", "point", "o",
    "square", "rect", "rectangle", "box",
    "diamond", "diam",
    "triangle", "triangle_up", "tri", "triangleup", "^",
    "triangle_down", "triangle_dn", "tri_dn", "tri_down",
    "invtriangle", "inv_triangle", "v"
  )
  if (any(!known)) {
    unk <- lev[!known]
    message(sprintf(
      "Unknown shape label(s) %s — drawn as circle (pch 21).",
      paste(unique(unk), collapse = ", ")
    ))
  }
  list(levels = lev, pch = stats::setNames(shape_pch_vectorized(lev), lev))
}

EDGE_TOKEN_COLOURS <- stats::setNames(
  c("#E8923C", "#E8C547", "#8A8A8A", "#8A8A8A"),
  tolower(c("orange", "yellow", "grey", "gray"))
)

edge_colour_identity <- function(col_vec) {
  raw <- ifelse(is.na(col_vec), "grey", trimws(as.character(col_vec)))
  raw[!nzchar(raw)] <- "grey"
  
  canon_key <- tolower(gsub("\\s+", "", raw))
  mapped <- EDGE_TOKEN_COLOURS[canon_key]
  v <- ifelse(!is.na(mapped), mapped, raw)
  
  is_hex <- grepl("^#[[:xdigit:]]{6}([[:xdigit:]]{2})?$", v) |
    grepl("^#[[:xdigit:]]{3}$", v)
  is_named <- v %in% grDevices::colors()
  ok <- is_hex | is_named
  if (any(!ok)) {
    bad <- unique(v[!ok])
    message(sprintf(
      "Non-colour edge value(s) %s — forcing neutral grey (#8A8A8A).",
      paste(utils::head(bad, 8), collapse = ", ")
    ))
    v[!ok] <- "#8A8A8A"
  }
  v
}

# --- Plot --------------------------------------------------------------------

for (ID in c("4E", "4G")) {
  edges_path <- paste0("../data/data.Figure", ID, ".edges.txt")
  nodes_path <- paste0("../data/data.Figure", ID, ".nodes.txt")
  out_path <- paste0("Figure", ID, ".png")

  edges_tbl <- read_edges_tsv(edges_path)
nodes_tbl <- read_nodes_tsv(nodes_path)
nodes_tbl <- add_missing_nodes(nodes_tbl, edges_tbl)
nodes_tbl <- nodes_tbl[!duplicated(nodes_tbl$name), , drop = FALSE]

nodes_tbl <- assign_node_layout(nodes_tbl, edges_tbl)
nodes_tbl$shape[toupper(nodes_tbl$name) %in% KINASE_IDS] <- "square"

edges <- edges_tbl
edges$x <- nodes_tbl$x[match(edges$from, nodes_tbl$name)]
edges$y <- nodes_tbl$y[match(edges$from, nodes_tbl$name)]
edges$ex <- nodes_tbl$x[match(edges$to, nodes_tbl$name)]
edges$ey <- nodes_tbl$y[match(edges$to, nodes_tbl$name)]
if (any(is.na(edges$x) | is.na(edges$y) | is.na(edges$ex) | is.na(edges$ey))) {
  stop("Edge endpoints must exist in merged node table.")
}

edges$linetype <- resolve_linetype(ifelse(is.na(edges$linestyle), "dashed", edges$linestyle))
edges$effect <- interaction_effect_from_type(edges$type)
edges$edge_id <- seq_len(nrow(edges))

for (j in seq_len(nrow(edges))) {
  s <- shrink_segment(edges$x[j], edges$y[j], edges$ex[j], edges$ey[j])
  edges$x[j] <- s$x
  edges$y[j] <- s$y
  edges$ex[j] <- s$xend
  edges$ey[j] <- s$yend
}

edges$colour_char <- edge_colour_identity(edges$colour)

edges_act <- edges[edges$effect == "activate", , drop = FALSE]
edges_inh <- edges[edges$effect == "inhibit", , drop = FALSE]

nlp <- nodes_tbl$neg_log10P
nlp_fin <- nlp[is.finite(nlp)]
if (!length(nlp_fin)) {
  nodes_tbl$nlp_plot <- rep(1, nrow(nodes_tbl))
  message("All neg_log10P values NA — node area uses uniform size (=1).")
} else {
  ref <- suppressWarnings(stats::median(nlp_fin))
  if (!is.finite(ref)) ref <- mean(nlp_fin)
  tiny <- max(stats::median(nlp_fin, na.rm = TRUE) * 0.08, suppressWarnings(min(nlp_fin)) * 0.5, 1e-9)
  nodes_tbl$nlp_plot <- ifelse(is.finite(nlp), nlp, ref)
  nodes_tbl$nlp_plot <- pmax(nodes_tbl$nlp_plot, tiny)
}

nodes_tbl$log2FC_plot <- ifelse(
  is.finite(nodes_tbl$log2FC),
  pmax(pmin(nodes_tbl$log2FC, 2), -2),
  NA_real_
)

shp_sc <- build_shape_scale(nodes_tbl$shape)
nodes_tbl$shape_fac <- factor(nodes_tbl$shape, levels = shp_sc$levels)

kin_mask <- toupper(nodes_tbl$name) %in% KINASE_IDS
nodes_main <- nodes_tbl[!kin_mask, , drop = FALSE]
nodes_kin <- nodes_tbl[kin_mask, , drop = FALSE]

gradient_scale <- scale_fill_gradient2(
  name = expression(log[2]~FC),
  low = "#699BF8",
  mid = "#F7F7F7",
  high = "#EB4C42",
  midpoint = 0,
  limits = c(-2, 2),
  na.value = "#d9d9d9",
  aesthetics = "fill",
  guide = "colourbar"
)

arrow_activate <- ggplot2::arrow(length = grid::unit(2.65, "mm"), type = "closed", angle = 22)
arrow_inhibit <- ggplot2::arrow(length = grid::unit(3.0, "mm"), type = "closed", angle = 90)

p <- ggplot() +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.box = "vertical",
    legend.position = "right"
  )

if (nrow(edges_act) > 0L) {
  p <- p +
    geom_segment(
      data = edges_act,
      aes(
        x = x,
        y = y,
        xend = ex,
        yend = ey,
        colour = colour_char,
        linetype = linetype,
        group = edge_id
      ),
      linewidth = 0.45,
      lineend = "round",
      arrow = arrow_activate,
      inherit.aes = FALSE
    )
}

if (nrow(edges_inh) > 0L) {
  p <- p +
    geom_segment(
      data = edges_inh,
      aes(
        x = x,
        y = y,
        xend = ex,
        yend = ey,
        colour = colour_char,
        linetype = linetype,
        group = edge_id
      ),
      linewidth = 0.45,
      lineend = "round",
      arrow = arrow_inhibit,
      inherit.aes = FALSE
    )
}

p <- p +
  scale_colour_identity(guide = "none") +
  scale_linetype_identity(guide = "none")

if (nrow(nodes_main) > 0L) {
  p <- p +
    geom_point(
      data = nodes_main,
      aes(
        x = x,
        y = y,
        fill = log2FC_plot,
        size = nlp_plot,
        shape = shape_fac
      ),
      colour = "#333333",
      stroke = 0.35,
      inherit.aes = FALSE
    )
}

if (nrow(nodes_kin) > 0L) {
  p <- p +
    geom_point(
      data = nodes_kin,
      aes(x = x, y = y, size = nlp_plot),
      shape = 22,
      fill = KINASE_FILL,
      colour = "#2b2b2b",
      stroke = 0.4,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
}

if (nrow(nodes_tbl) > 0L) {
  p <- p +
    scale_size_area(max_size = 12, name = expression(-log[10] * italic(P)))
}

if (nrow(nodes_main) > 0L) {
  p <- p +
    gradient_scale +
    scale_shape_manual(values = shp_sc$pch, name = "Shape", guide = "none")
}

p <- p +
  geom_text(
    data = nodes_tbl,
    aes(x = x, y = y, label = name),
    size = 3.2,
    fontface = "plain",
    colour = "#151515",
    inherit.aes = FALSE
  )

suppressMessages({
  ggplot2::ggsave(
    filename = out_path,
    plot = p,
    width = 9.5,
    height = 7,
    dpi = 300,
    bg = "white"
  )
})

  message(sprintf("Wrote %s", out_path))
}

# Done.