create_contrast_plot <- function(emmeans_results, title="", box_wide=0.3, 
                                 group_color=c("#4682B4", "#FF6347", "#3CB371", "#9370DB"),
                                 show_x_labels=F) {
  data <- as.data.frame(emmeans_results$emmeans)
  colnames(data)[5:6] <- c("lower.CL","upper.CL")
  contrasts <- as.data.frame(emmeans_results$contrasts)
  contrasts <- contrasts %>%
    separate(contrast, into = c("Group1", "Group2"), sep = " - ") %>%
    mutate(
      x1 = as.numeric(factor(Group1, levels = levels(data$Group))),
      x2 = as.numeric(factor(Group2, levels = levels(data$Group))),
      y = max(data$upper.CL) + 0.02*(max(data$upper.CL)-min(data$lower.CL)), 
      label = ifelse(p.value < 0.001, "***",
                     ifelse(p.value < 0.01, "**",
                            ifelse(p.value < 0.05, "*", "")))
    ) %>%
    filter(label != "") %>% 
    arrange(desc(p.value)) 
  
  p <- ggplot(data, aes(x = Group, fill = Group)) +
    geom_rect(
      aes(
        xmin = as.numeric(Group) - box_wide, 
        xmax = as.numeric(Group) + box_wide,
        ymin = lower.CL,
        ymax = upper.CL
      ),
      alpha = 0.7, color = "black"
    ) +
    geom_segment(
      aes(
        x = as.numeric(Group) - box_wide,
        xend = as.numeric(Group) + box_wide,
        y = emmean,
        yend = emmean
      ),
      color = "black", size = 1
    ) +
    scale_fill_manual(values = group_color) + # Colors for groups
    scale_x_discrete(labels = levels(data$Group)) + # Keep Group as discrete
    labs(
      title = title,
      x = NULL,
      y = "Estimated Marginal Means"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16,family = "Arial"),
      axis.title.y = element_text(size = 14,family = "Arial"),
      axis.text.y = element_text(size = 12), # Tilt x-axis labels
      axis.text.x  = if (show_x_labels)
        ggplot2::element_text(size = 14, angle = 45, hjust = 1, family = "Arial")
      else
        ggplot2::element_blank(),
      legend.position = if (show_x_labels) {"none"} else {"right"}
    )
  
  if (exists("contrasts") && nrow(contrasts) > 0) {
    p <- p +
      geom_segment(
        data = contrasts,
        aes(
          x = x1,
          xend = x1,
          y = y + 0.01*(max(data$upper.CL)-min(data$lower.CL)),
          yend = y + 0.03*(max(data$upper.CL)-min(data$lower.CL))
        ),
        inherit.aes = FALSE, color = "black", size = 0.8
      ) +
      geom_segment(
        data = contrasts,
        aes(
          x = x2,
          xend = x2,
          y = y + 0.01*(max(data$upper.CL)-min(data$lower.CL)),
          yend = y + 0.03*(max(data$upper.CL)-min(data$lower.CL))
        ),
        inherit.aes = FALSE, color = "black", size = 0.8
      ) +
      geom_segment(
        data = contrasts,
        aes(
          x = x1,
          xend = x2,
          y = y + 0.03*(max(data$upper.CL)-min(data$lower.CL)),
          yend = y + 0.03*(max(data$upper.CL)-min(data$lower.CL))
        ),
        inherit.aes = FALSE, color = "black", size = 0.8
      ) +
      geom_text(
        data = contrasts,
        aes(
          x = (x1 + x2) / 2, 
          y = y + 0.04*(max(data$upper.CL)-min(data$lower.CL)), 
          label = label
        ),
        inherit.aes = FALSE,
        size = 5
      )
  }
  p
}