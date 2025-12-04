# stat_emm_stars(): add FDR stars from emmeans to facetted plots
# Usage:
#   ct <- pairs(emmeans(m, ~ group | time), adjust = "BH") %>% summary()
#   ggplot(df, aes(group, percent)) +
#     geom_boxplot() +
#     facet_wrap(~ time) +
#     stat_emm_stars(
#       aes(time_panel = time),    # <-- carry facet var into the stat
#       emm_pairs_df = ct,
#       p_col   = "p.value",
#       time_col = "time"
#     )

stat_emm_stars <- function(
    mapping = NULL,
    data = NULL,
    geom = "text",
    position = "identity",
    ...,
    emm_pairs_df,
    time_col   = "time",                  # column in emm_pairs_df with facet value
    p_col      = "p.value",               # FDR-adjusted p-values column in emm_pairs_df
    cutpoints  = c(0, 1e-4, 1e-3, 1e-2, 0.05, Inf),
    symbols    = c("****", "***", "**", "*", "ns"),
    x_pos      = "middle",                # "middle" or numeric x
    y_pad      = 0.05,                    # 0..1 -> fraction of y-range; >1 absolute
    size       = 4,
    vjust      = 0,
    contrast_filter = NULL,               # optional function to pre-filter emm_pairs_df
    na.rm = TRUE,
    show.legend = FALSE,
    inherit.aes = TRUE
) {
  stopifnot(length(symbols) == length(cutpoints) - 1)
  ep <- emm_pairs_df
  if (!is.null(contrast_filter)) ep <- contrast_filter(ep)
  
  ggplot2::layer(
    stat = StatEmmStars2, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      ep = ep, time_col = time_col, p_col = p_col,
      cutpoints = cutpoints, symbols = symbols,
      x_pos = x_pos, y_pad = y_pad,
      size = size, vjust = vjust, na.rm = na.rm, ...
    )
  )
}

StatEmmStars2 <- ggplot2::ggproto(
  "StatEmmStars2", ggplot2::Stat,
  # We need: x, y (to get ranges), and a dummy aesthetic 'time_panel' that carries the facet value
  required_aes = c("x", "y", "time_panel"),
  compute_panel = function(data, scales,
                           ep, time_col, p_col,
                           cutpoints, symbols, x_pos, y_pad,
                           size, vjust, na.rm) {
    
    # Extract this panel's time value (carried via aes(time_panel = time))
    time_vals <- unique(data$time_panel)
    time_vals <- time_vals[!is.na(time_vals)]
    if (length(time_vals) == 0L) return(data.frame())
    time_val <- time_vals[1]
    
    # Pull the (adjusted) p-value for this facet from emmeans pairs table
    if (!time_col %in% names(ep))
      stop(sprintf("stat_emm_stars: '%s' not found in emm_pairs_df.", time_col))
    if (!p_col %in% names(ep))
      stop(sprintf("stat_emm_stars: '%s' not found in emm_pairs_df.", p_col))
    
    ep_sub <- ep[ep[[time_col]] == time_val, , drop = FALSE]
    if (nrow(ep_sub) == 0L) return(data.frame())
    pval <- ep_sub[[p_col]][1]
    
    # Convert p to stars
    stars <- as.character(cut(pval, breaks = cutpoints, labels = symbols, include.lowest = TRUE))
    
    # y position (padding above panel max)
    y_max <- max(data$y, na.rm = TRUE)
    y_min <- min(data$y, na.rm = TRUE)
    y_annot <- if (!is.na(y_pad) && y_pad >= 0 && y_pad <= 1) {
      y_max + y_pad * (y_max - y_min)
    } else {
      y_max + if (is.na(y_pad)) 0 else y_pad
    }
    
    # x position (center between groups unless numeric)
    x_center <- if (identical(x_pos, "middle")) {
      rng <- range(unique(data$x))
      mean(rng)
    } else if (is.numeric(x_pos) && length(x_pos) == 1) {
      x_pos
    } else {
      stop("x_pos must be 'middle' or a single numeric value.")
    }
    
    data.frame(
      x = x_center,
      y = y_annot,
      label = stars
    )
  }
)

# model_fun <- function(df) {
#   # Read & prep
#   df <- df |>
#     rename(group = patient_type3,
#            time  = time,
#            percent = percent) |>
#     mutate(
#       patient = factor(patient),
#       group   = factor(group),
#       # levels like "IBR", "IBS"
#       time    = factor(time, levels = c(1, 2, 3))
#     )
#   
#   # Linear mixed model: random intercept for patient, fixed group*time
#   m <- lmer(percent ~ group * time + (1 | patient), data = df)
#   
#   # Model summary and ANOVA table
#   summary(m)
#   anova(m)  # Type III by default in lmerTest
#   
#   # Estimated marginal means and by-time group contrasts (IBR vs IBS at each time)
#   emm <- emmeans(m, ~ group | time)
#   pairs <- pairs(emm, adjust = "none")   # FDR-BH across the three timepoints
#   
#   # If you also want the estimated differences with CIs:
#   # contrast(emm, method = "revpairwise", adjust = "BH")  # IBR - IBS per time
#   # confint(contrast(emm, method = "revpairwise", adjust = "BH"))
#   
#   list(model = m, pairs = pairs)
#   
# }

# model_fun <- function(df) {
#   # --- 1. Prep data ----------------------------------------------------------
#   df <- df |>
#     rename(group = patient_type3,
#            time  = time,
#            percent = percent) |>
#     mutate(
#       patient = factor(patient),
#       group   = factor(group),
#       time    = factor(time, levels = c(1, 2, 3))
#     )
#   
#   # --- 2. Fit model ----------------------------------------------------------
#   m <- lmer(percent ~ group * time + (1 | patient), data = df)
#   
#   # --- 3. Global ANOVA tests -------------------------------------------------
#   anova_table <- anova(m, type = 3) |> 
#     tibble::rownames_to_column("Term")
#   
#   # --- 4. Estimated marginal means ------------------------------------------
#   emm <- emmeans(m, ~ group | time)
#   
#   # --- 5. Pairwise contrasts (IBR - IBS within each time) -------------------
#   ct_unadj <- pairs(emm, adjust = "none") |> as.data.frame()
#   ct_unadj <- ct_unadj |>
#     rename(p_unadjusted = p.value)
#   
#   # Apply FDR (BH) across all timepoints collectively
#   ct_unadj <- ct_unadj |>
#     mutate(p_FDR_BH = p.adjust(p_unadjusted, method = "BH"),
#            sig_raw  = cut(p_unadjusted,
#                           breaks = c(0, 1e-4, 1e-3, 1e-2, 0.05, Inf),
#                           labels = c("****", "***", "**", "*", "ns"),
#                           include.lowest = TRUE),
#            sig_FDR  = cut(p_FDR_BH,
#                           breaks = c(0, 1e-4, 1e-3, 1e-2, 0.05, Inf),
#                           labels = c("****", "***", "**", "*", "ns"),
#                           include.lowest = TRUE))
#   
#   # --- 6. Output -------------------------------------------------------------
#   list(
#     model        = m,
#     anova_table  = anova_table,
#     contrasts    = ct_unadj
#   )
# }