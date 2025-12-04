
stat_stars <- function(mapping = NULL,
                       data = NULL,
                       geom = "text",
                       position = "identity",
                       ...,
                       stat_df,
                       time_col   = "time",
                       p_col      = "p.value",
                       cutpoints  = c(0, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                       symbols    = c("****", "***", "**", "*", "ns"),
                       npc_x      = 0.5,
                       npc_y      = 0.95,
                       size       = 4,
                       vjust      = 0,
                       na.rm = TRUE,
                       show.legend = FALSE,
                       inherit.aes = TRUE) {
  ggplot2::layer(
    stat = StatStars2,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      stat_df = stat_df,
      time_col = time_col,
      p_col = p_col,
      cutpoints = cutpoints,
      symbols = symbols,
      npc_x = npc_x,
      npc_y = npc_y,
      size = size,
      vjust = vjust,
      na.rm = na.rm,
      ...
    )
  )
}
# 
# StatStars2 <- ggplot2::ggproto(
#   "StatStars2", ggplot2::Stat,
#   # We need: x, y (to get ranges), and a dummy aesthetic 'time_panel' that carries the facet value
#   required_aes = c("x", "y", "time_panel"),
#   compute_panel = function(data, scales,
#                            stat_df, time_col, p_col,
#                            cutpoints, symbols, npc_x, npc_y,
#                            size, vjust, na.rm) {
#     
#     # Extract this panel's time value (carried via aes(time_panel = time))
#     time_vals <- unique(data$time_panel)
#     time_vals <- time_vals[!is.na(time_vals)]
#     if (length(time_vals) == 0L) return(data.frame())
#     time_val <- time_vals[1]
#     
#     # Pull the (adjusted) p-value for this facet from emmeans pairs table
#     if (!time_col %in% names(stat_df))
#       stop(sprintf("stat_stars: '%s' not found in emm_pairs_df.", time_col))
#     if (!p_col %in% names(stat_df))
#       stop(sprintf("stat_stars: '%s' not found in emm_pairs_df.", p_col))
#     
#     df_sub <- stat_df[stat_df[[time_col]] == time_val, , drop = FALSE]
#     if (nrow(df_sub) == 0L) return(data.frame())
#     pval <- df_sub[[p_col]][1]
#     
#     # Convert p to stars
#     stars <- as.character(cut(pval, breaks = cutpoints, labels = symbols, include.lowest = TRUE))
#     
#     # Panel data ranges (x may be discrete but is numeric at this stage)
#     x_min <- min(data$x, na.rm = TRUE)
#     x_max <- max(data$x, na.rm = TRUE)
#     y_min <- min(data$y, na.rm = TRUE)
#     y_max <- max(data$y, na.rm = TRUE)
#     
#     # NPC -> data coords (linear interpolation across panel range)
#     x_ann <- x_min + npc_x * (x_max - x_min)
#     y_ann <- y_min + npc_y * (y_max - y_min)
#     
#     data.frame(
#       x = x_ann,
#       y = y_ann,
#       label = stars
#     )
#   }
# )
stat_stars_facet <- function(mapping = NULL,
                       data = NULL,
                       geom = "text",
                       position = "identity",
                       ...,
                       stat_df,
                       time_col   = "time",
                       p_col      = "p.value",
                       cutpoints  = c(0, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                       symbols    = c("****", "***", "**", "*", "ns"),
                       npc_x      = 0.5,
                       npc_y      = 0.95,
                       size       = 4,
                       vjust      = 0,
                       na.rm = TRUE,
                       show.legend = FALSE,
                       inherit.aes = TRUE) {

  ggplot2::layer(
    stat = StatStars2,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      stat_df = stat_df,
      time_col = time_col,
      p_col = p_col,
      cutpoints = cutpoints,
      symbols = symbols,
      npc_x = npc_x,
      npc_y = npc_y,
      size = size,
      vjust = vjust,
      na.rm = na.rm,
      ...
    )
  )
}

StatStars2 <- ggplot2::ggproto(
  "StatStars2", ggplot2::Stat,
  required_aes = c("x", "y", "time_panel"),
  compute_panel = function(data, scales,
                           stat_df, time_col, p_col,
                           cutpoints, symbols, npc_x, npc_y,
                           size, vjust, na.rm) {

    # facet value carried via aes(time_panel = time)
    time_vals <- unique(data$time_panel)
    time_vals <- time_vals[!is.na(time_vals)]
    if (length(time_vals) == 0L) return(data.frame())
    time_val <- time_vals[1]

    # pick p-value
    if (!time_col %in% names(stat_df)) stop(sprintf("stat_stars: '%s' not in stat_df.", time_col))
    if (!p_col %in% names(stat_df))   stop(sprintf("stat_stars: '%s' not in stat_df.", p_col))
    df_sub <- stat_df[stat_df[[time_col]] == time_val, , drop = FALSE]
    if (nrow(df_sub) == 0L) return(data.frame())
    pval <- df_sub[[p_col]][1]
    stars <- as.character(cut(pval, breaks = cutpoints, labels = symbols, include.lowest = TRUE))

    # ----- panel ranges from SCALES (preferred) -----
    # y-range (continuous)
    y_rng <- tryCatch(scales$y$range$range, error = function(e) NULL)
    if (is.null(y_rng) || any(!is.finite(y_rng))) {
      # fallback to layer data range
      y_rng <- range(data$y, na.rm = TRUE)
    }

    # x-range: handle continuous vs discrete
    x_rng <- tryCatch(scales$x$range$range, error = function(e) NULL)

    # If discrete, 'range' may be character levels. Convert to numeric positions 1..k
    if (is.null(x_rng)) {
      # fallback
      xs <- sort(unique(data$x))
      x_min <- min(xs, na.rm = TRUE); x_max <- max(xs, na.rm = TRUE)
    } else if (is.numeric(x_rng)) {
      x_min <- x_rng[1]; x_max <- x_rng[2]
    } else {
      # character / factor levels:
      k <- length(x_rng)
      x_min <- 1
      x_max <- k
    }

    # Convert NPC -> data coords using panel scale ranges
    x_ann <- x_min + npc_x * (x_max - x_min)
    y_ann <- y_rng[1] + npc_y * (y_rng[2] - y_rng[1])

    data.frame(
      x = x_ann,
      y = y_ann,
      # label = pval
      label = stars
    )
  }
)
stat_stars_grid <- function(mapping = NULL,
                       data = NULL,
                       geom = "text",
                       position = "identity",
                       ...,
                       stat_df,
                       time_col     = "time",      # column in stat_df
                       cluster_col  = NULL,        # set to "cluster" if using facet_grid(cluster ~ time)
                       p_col        = "p.value",   # column with p-values/q-values
                       cutpoints    = c(0, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                       symbols      = c("****", "***", "**", "*", "ns"),
                       npc_x        = 0.5,
                       npc_y        = 0.95,
                       size         = 4,
                       vjust        = 0,
                       na.rm = TRUE,
                       show.legend = FALSE,
                       inherit.aes = TRUE) {
  
  ggplot2::layer(
    stat = StatStarsGrid,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      stat_df = stat_df,
      time_col = time_col,
      cluster_col = cluster_col,
      p_col = p_col,
      cutpoints = cutpoints,
      symbols = symbols,
      npc_x = npc_x,
      npc_y = npc_y,
      size = size,
      vjust = vjust,
      na.rm = na.rm,
      ...
    )
  )
}

StatStarsGrid <- ggplot2::ggproto(
  "StatStarsGrid", ggplot2::Stat,
  # x,y needed for panel ranges; we *optionally* accept panel_time/panel_cluster via mapping
  required_aes = c("x", "y", "panel_time", "panel_cluster"),
  compute_panel = function(data, scales,
                           stat_df, time_col, cluster_col, p_col,
                           cutpoints, symbols, npc_x, npc_y,
                           size, vjust, na.rm) {
    
    # --- identify facet keys present in this panel ---
    get_unique_safe <- function(v) {
      u <- unique(v)
      u <- u[!is.na(u)]
      if (length(u) == 0L) return(NULL)
      u[1]
    }
    
    time_val    <- if ("panel_time"    %in% names(data)) get_unique_safe(data$panel_time)    else NULL
    cluster_val <- if ("panel_cluster" %in% names(data)) get_unique_safe(data$panel_cluster) else NULL
    
    # Validate that necessary keys are present
    if (!is.null(cluster_col) && is.null(cluster_val)) {
      stop("stat_stars: facet_grid detected (cluster_col provided) but 'panel_cluster' aesthetic not mapped.")
    }
    if (is.null(time_val)) {
      stop("stat_stars: 'panel_time' aesthetic must be mapped (e.g., aes(panel_time = time)).")
    }
    
    # --- subset stat_df by the facet keys ---
    if (!time_col %in% names(stat_df)) stop(sprintf("stat_stars: '%s' not in stat_df.", time_col))
    if (!p_col   %in% names(stat_df)) stop(sprintf("stat_stars: '%s' not in stat_df.", p_col))
    sdf <- stat_df
    
    sdf <- sdf[sdf[[time_col]] == time_val, , drop = FALSE]
    if (!is.null(cluster_col)) {
      if (!cluster_col %in% names(sdf)) stop(sprintf("stat_stars: '%s' not in stat_df.", cluster_col))
      sdf <- sdf[sdf[[cluster_col]] == cluster_val, , drop = FALSE]
    }
    if (nrow(sdf) == 0L) return(data.frame())
    
    pval  <- sdf[[p_col]][1]
    stars <- as.character(cut(pval, breaks = cutpoints, labels = symbols, include.lowest = TRUE))
    
    # --- panel scale ranges (consistent across layers & facets) ---
    y_rng <- tryCatch(scales$y$range$range, error = function(e) NULL)
    if (is.null(y_rng) || any(!is.finite(y_rng))) y_rng <- range(data$y, na.rm = TRUE)
    
    x_rng <- tryCatch(scales$x$range$range, error = function(e) NULL)
    if (is.null(x_rng)) {
      xs <- sort(unique(data$x))
      x_min <- min(xs, na.rm = TRUE); x_max <- max(xs, na.rm = TRUE)
    } else if (is.numeric(x_rng)) {
      x_min <- x_rng[1]; x_max <- x_rng[2]
    } else { # discrete levels (character)
      k <- length(x_rng)
      x_min <- 1; x_max <- k
    }
    
    # --- NPC -> data coords ---
    x_ann <- x_min + npc_x * (x_max - x_min)
    y_ann <- y_rng[1] + npc_y * (y_rng[2] - y_rng[1])
    
    data.frame(
      x = x_ann,
      y = y_ann,
      label = stars
    )
  }
)
