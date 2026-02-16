library(ggplot2)
# Refined function for plotting intercept values with CI & significance stars
draw_intercept <- function(
    diff_stats,
    palette = NULL,
    box_width = 0.3,
    title = "Adjusted HD–SD Difference",
    ylab = "Adjusted Mean Difference (95% CI)"
) {
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  
  # Ensure Antigen is a factor
  diff_stats <- diff_stats %>%
    mutate(Antigen = factor(Antigen))
  
  # Determine palette: use provided or generate
  if (is.null(palette)) {
    n <- nlevels(diff_stats$Antigen)
    cols <- brewer.pal(max(3, min(n, 8)), "Set2")[1:n]
    palette <- setNames(cols, levels(diff_stats$Antigen))
  } else {
    if (is.null(names(palette))) {
      palette <- setNames(rep(palette, length.out = nlevels(diff_stats$Antigen)),
                          levels(diff_stats$Antigen))
    } else {
      palette <- palette[levels(diff_stats$Antigen)]
    }
  }
  
  # Prepare star data only for significant
  star_df <- diff_stats %>%
    filter(star != "") %>%
    mutate(
      x = as.numeric(Antigen),
      y_pos = conf.high * 1.05
    )
  
  # Build plot
  p <- ggplot(diff_stats, aes(x = Antigen, fill = Antigen)) +
    # CI as rectangles
    geom_rect(
      aes(
        xmin = as.numeric(Antigen) - box_width,
        xmax = as.numeric(Antigen) + box_width,
        ymin = conf.low,
        ymax = conf.high
      ),
      alpha = 0.7, color = "black"
    ) +
    # Mean as segment
    geom_segment(
      aes(
        x    = as.numeric(Antigen) - box_width,
        xend = as.numeric(Antigen) + box_width,
        y    = estimate,
        yend = estimate
      ),
      color = "black", size = 1
    ) +
    # Significance stars
    geom_text(
      data = star_df,
      aes(x = x, y = y_pos, label = star),
      family = "serif",
      size   = 6,
      inherit.aes = FALSE
    ) +
    # Scales & labels
    scale_fill_manual(values = palette) +
    scale_x_discrete(labels = levels(diff_stats$Antigen)) +
    labs(
      title = title,
      x     = NULL,
      y     = ylab
    ) +
    # Theme
    theme_bw() +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 14, family = "serif"),
      axis.title.y       = element_text(size = 12, family = "serif"),
      axis.text.x        = element_text(size = 10, angle = 45, hjust = 1, family = "serif"),
      legend.position    = "none",
      panel.grid.major.x = element_blank()
    )
  
  return(p)
}