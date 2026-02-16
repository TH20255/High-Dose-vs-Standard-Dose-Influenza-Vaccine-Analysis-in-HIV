library(flextable)
library(gtsummary)
library(tidyverse)
library(readxl)
library(car)

######### Data Preprocessing #########
df_long <- read_excel("./serology_master_dataset.xlsx")
df_long <- df_long %>% distinct()

df_long$Race <- as.factor(df_long$Race)
df_long <- mutate(df_long, Race = fct_collapse(Race,
                                               `Black/African American` = c("Black/African American", "Black/African American,Multi-racial/Mixed","Black/African American,White/Caucasian"),
                                               `White/Caucasian` = c("White/Caucasian", "White/Caucasian,Other"),
                                               Others = c("Asian","Multi-racial/Mixed","Native American","Other")))

colnames(df_long)[7] <- "previous flu vaccine"
colnames(df_long)[15] <- "HAI_Titer"
colnames(df_long)[16] <- "HA_IgG"


df_long$`Age group` <- as.factor(df_long$`Age group`)
df_long$Sex <- as.factor(df_long$Sex)
df_long$`HIV status` <- as.factor(df_long$`HIV status`)
df_long$Race <- as.factor(df_long$Race)
df_long$Ethnicity <- as.factor(df_long$Ethnicity)
df_long$`previous flu vaccine` <- as.factor(df_long$`previous flu vaccine`)

df_long <- df_long %>%
  mutate(Group = case_when(
    Group == "AOC" ~ "Older PWoH",
    Group == "AOH" ~ "Older PWH",
    Group == "AYC" ~ "Young PWoH",
    Group == "AYH" ~ "Young PWH"
  )) %>% 
  mutate(Antigen = case_when(
    Antigen == "H1N1" ~ "A/H1N1",
    Antigen == "H3N2" ~ "A/H3N2",
    Antigen == "B1" ~ "B/Victoria",
    Antigen == "B2" ~ "B/Yamagata"
  ))

df_long_HAI <- df_long[!is.na(df_long$HAI_Titer),]
df_long_HAI <- df_long_HAI[,-c(6,16)]
df_long_IgG <- df_long[!is.na(df_long$HA_IgG),]
df_long_IgG <- df_long_IgG[,-c(6,15)]


######### HAI Titer multivariate analysis #########
library(MASS)
df_wide_HAI <- df_long_HAI %>%
  pivot_wider(names_from = c(Timepoint), 
              values_from = HAI_Titer)
HAI_path <- "./outputs/HAI/"

library(patchwork)
all_contrast_plots <- list()
#### by antigen_type and dose type 
for(antigen_type in c("A/H1N1","B/Victoria","B/Yamagata"))
{
  data_ordinal0 <- df_wide_HAI %>% filter(Antigen == antigen_type)
  for(dose_type in c("SD","HD")) 
  {
    data_ordinal <- data_ordinal0 %>% filter(Dose == dose_type)

    data_ordinal <- data_ordinal[c("Group", "Sex", "Race",  "Ethnicity", "previous flu vaccine", "T0","T3")]
    data_ordinal <- na.omit(data_ordinal)
    
    data_ordinal <- droplevels(data_ordinal)

    data_ordinal$T3 <- as.factor(data_ordinal$T3)
    
    model_ordinal <- polr(T3 ~ ., 
                          data = data_ordinal, 
                          method = "logistic")
    
    ######## contrast plot ########
    library(emmeans)
    source("./contrast_for_ordinalrg.R")
    emmeans_results <- emmeans(model_ordinal, pairwise ~ Group)
    contrast_plot <- create_contrast_plot(emmeans_results, 
                                          title = paste0(antigen_type), 
                                          group_color=c("#9370DB","#3CB371","#FF6347","#4682B4"))
    all_contrast_plots[[which(antigen_type==c("A/H1N1","B/Victoria","B/Yamagata"))+3*(which(dose_type==c("SD","HD"))-1)]] <- contrast_plot
    
    model_summary <- summary(model_ordinal)
    se <- coef(model_summary)[, "Std. Error"]
    
    tbl <- tbl_regression(model_ordinal, 
                          exponentiate = TRUE) %>% 
      add_global_p() %>%
      bold_p(t = 0.05) %>% 
      add_n(location = c('label', 'level')) %>%
      modify_header(label = "**Predictors**") %>%
      modify_spanning_header(everything() ~ paste0("**","HAI Titer ",dose_type," ",antigen_type,"**")) %>% 
      bold_labels() %>%  
      as_flex_table()
    
    library(officer)
    plot_path <- paste0(HAI_path,"contrast_plot.png")
    ggsave(plot_path, plot = contrast_plot, width = 6, height = 4.67, dpi = 300)
    doc <- read_docx()
    
    doc <- doc %>%
      body_add_par("Regression Table", style = "table title")
    
    doc <- doc %>%
      body_add_flextable(tbl)
    
    doc <- doc %>%
      body_add_par("Contrast Plots", style = "graphic title")
    
    doc <- doc %>%
      body_add_img(src = plot_path, width = 4, height = 3, style = "centered")
    
    print(doc, target = paste0(HAI_path,"HAI ",gsub("/", "_",paste0(dose_type,"_",antigen_type)),".docx"))
  }
}

HAI_contrast_plots <- all_contrast_plots
HAI_SD_contrast_plots <- HAI_contrast_plots[1:3]
HAI_HD_contrast_plots <- HAI_contrast_plots[4:6]

combined_plot <- wrap_plots(HAI_contrast_plots, ncol = 2)
print(combined_plot)

######### IgG multivariate analysis #########
library(MASS)
df_wide_IgG <- df_long_IgG %>%
  pivot_wider(names_from = c(Timepoint), 
              values_from = HA_IgG)
IgG_path <- "./outputs/IgG/"

#### by antigen_type and dose type 
all_contrast_plots <- list()
for(antigen_type in c("A/H1N1","A/H3N2","B/Victoria","B/Yamagata"))
{
  data_ordinal0 <- df_wide_IgG %>% filter(Antigen == antigen_type)
  for(dose_type in c("SD","HD")) 
  {
    data_ordinal <- data_ordinal0 %>% filter(Dose == dose_type)
    
    data_ordinal <- data_ordinal[c("Group", "Sex" , "Race" ,  "Ethnicity", "previous flu vaccine", "T0","T3")] 
    
    data_ordinal <- na.omit(data_ordinal)
    
    model_lm <- lm(T3 ~ ., 
                   data = data_ordinal)
    summary(model_lm)
    
    library(emmeans)
    source("./contrast_for_ordinalrg.R")
    emmeans_results <- emmeans(model_lm, pairwise ~ Group)
    contrast_plot <- create_contrast_plot(emmeans_results, 
                                          title = paste0(antigen_type), 
                                          group_color=c("#9370DB","#3CB371","#FF6347","#4682B4"))
    contrast_plot
    all_contrast_plots[[which(antigen_type==c("A/H1N1","A/H3N2","B/Victoria","B/Yamagata"))+4*(which(dose_type==c("SD","HD"))-1)]] <- contrast_plot
    
    tidy_with_power <- function(x,
                                conf.int      = FALSE,
                                conf.level    = 0.95,
                                exponentiate  = FALSE,
                                ...) {
      td <- broom::tidy(x,
                        conf.int     = conf.int,
                        conf.level   = conf.level,
                        exponentiate = exponentiate,
                        ...)
      
      df    <- x$df.residual
      alpha <- 1 - conf.level
      tcrit <- qt(1 - alpha/2, df)
      td <- td %>%
        dplyr::mutate(
          power = pt(-tcrit, df, ncp = statistic) +
            (1 - pt( tcrit, df, ncp = statistic))
        )
      td
    }
    
    tbl <- tbl_regression(model_lm, intercept = F) %>% 
      add_global_p() %>% 
      bold_p(t = 0.05) %>% 
      add_n(location = c('label', 'level')) %>%
      modify_header(label = "**Predictors**") %>%
      modify_spanning_header(everything() ~ paste0("**","IgG ",antigen_type,"**")) %>% 
      bold_labels() %>%  
      as_flex_table()
    
    library(officer)
    plot_path <- paste0(IgG_path,"contrast_plot.png")
    ggsave(plot_path, plot = contrast_plot, width = 6, height = 4.67, dpi = 300)
    doc <- read_docx()
    
    doc <- doc %>%
      body_add_par("Regression Table", style = "table title")
    
    doc <- doc %>%
      body_add_flextable(tbl)
    
    doc <- doc %>%
      body_add_par("Contrast Plots", style = "graphic title")
    
    doc <- doc %>%
      body_add_img(src = plot_path, width = 4, height = 3, style = "centered")
    
    print(doc, target = paste0(IgG_path,"IgG ",gsub("/", "_",paste0(dose_type,"_",antigen_type)),".docx"))
  }
}

IgG_contrast_plots <- all_contrast_plots
IgG_SD_contrast_plots <- IgG_contrast_plots[1:4]
IgG_HD_contrast_plots <- IgG_contrast_plots[5:8]

combined_plot <- wrap_plots(IgG_contrast_plots, ncol = 2)

print(combined_plot)

combined_sd_plot <- wrap_plots(c(IgG_SD_contrast_plots,HAI_SD_contrast_plots), ncol = 2) 
ggsave(paste0("./outputs/","SD_contrast_plot.png"), plot = combined_sd_plot, width = 10, height = 16, dpi = 300)

combined_hd_plot <- wrap_plots(c(IgG_HD_contrast_plots,HAI_HD_contrast_plots), ncol = 2) 
ggsave(paste0("./outputs/","HD_contrast_plot.png"), plot = combined_hd_plot, width = 10, height = 16, dpi = 300)

######### HAI HD vs. SD #########
df_wide_HAI <- df_long_HAI %>%
  dplyr::select(-`Season collected`,-Age) %>% 
  pivot_wider(names_from = c(Dose,Timepoint), 
              values_from = HAI_Titer)
HAI_paired_path <- "./outputs/HAI_paired/"

all_contrast_plots <- list()  
diff_list <- list()

for(antigen_type in c("A/H1N1","B/Victoria","B/Yamagata"))
{
  data_ordinal <- df_wide_HAI %>% 
    filter(Antigen == antigen_type) %>% 
    mutate(diff_T3=HD_T3-SD_T3)
  
    data_ordinal <- data_ordinal[c("Group", "Sex", "Race",  "Ethnicity", "previous flu vaccine", "SD_T0","diff_T3")]
  
  data_ordinal <- na.omit(data_ordinal)
  
  model_lm <- lm(diff_T3 ~ ., 
                 data = data_ordinal)
  summary(model_lm)
  
  td <- broom::tidy(model_lm, conf.int = TRUE)
  interc <- td %>% filter(term == "(Intercept)")
  est   <- interc$estimate
  low   <- interc$conf.low
  high  <- interc$conf.high
  pval  <- interc$p.value
  pval   <- summary(model_lm)$coefficients["(Intercept)", "Pr(>|t|)"]
  star <- if (pval < 0.001) {
    "***"
  } else if (pval < 0.01) {
    "**"
  } else if (pval < 0.05) {
    "*"
  } else {
    ""
  }
  diff_list[[antigen_type]] <- data.frame(
    Antigen   = antigen_type,
    estimate  = est,
    conf.low  = low,
    conf.high = high,
    p.value   = pval,
    star      = star,
    stringsAsFactors = FALSE
  )
  
  tbl <- tbl_regression(model_lm, intercept = TRUE) %>% 
    add_global_p() %>%
    bold_p(t = 0.05) %>% 
    add_n() %>%
    modify_header(label = "**Predictors**") %>%
    modify_spanning_header(everything() ~ paste0("**","HAI Titer ",antigen_type,"**")) %>% 
    bold_labels() %>%  # Bold the labels
    as_flex_table()
  
  library(officer)
  doc <- read_docx()
  
  doc <- doc %>%
    body_add_par("Regression Table", style = "table title")
  
  doc <- doc %>%
    body_add_flextable(tbl)
  
  print(doc, target = paste0(HAI_paired_path,"HAI Titer ",gsub("/", "_", paste0(antigen_type)),".docx"))
}

diff_stats <- do.call(rbind, diff_list)

source("./draw_intercept.R")
effect_plot <- draw_intercept(
  diff_stats,
  palette = c("#66C2A4","#8DA0CB","#E78AC3"),
  title   = "Adjusted HD–SD HAI Titer at 28dpv",
  ylab    = "Adjusted Mean Difference (95% CI)",
)
ggsave(paste0(HAI_paired_path,"HAITiter_adjusted_HD-SD.png"), plot = effect_plot, width = 6, height = 4.67, dpi = 300)

effect_plot

######### IgG HD vs. SD #########
df_wide_IgG <- df_long_IgG %>%
  dplyr::select(-`Season collected`,-Age) %>% 
  pivot_wider(names_from = c(Dose,Timepoint), 
              values_from = HA_IgG)
IgG_paired_path <- "./outputs/IgG_paired/"

all_contrast_plots <- list()  
diff_list <- list()

for(antigen_type in c("A/H1N1","A/H3N2","B/Victoria","B/Yamagata"))
{
  data_ordinal <- df_wide_IgG %>% 
    filter(Antigen == antigen_type) %>% 
    mutate(diff_T3=HD_T3-SD_T3)
  
  data_ordinal <- data_ordinal[c("Group", "Sex", "Race",  "Ethnicity", "previous flu vaccine", "SD_T0","diff_T3")]
  
  data_ordinal <- na.omit(data_ordinal)
  
  model_lm <- lm(diff_T3 ~ ., 
                 data = data_ordinal)
  summary(model_lm)
  
  td <- broom::tidy(model_lm, conf.int = TRUE)
  interc <- td %>% filter(term == "(Intercept)")
  est   <- interc$estimate
  low   <- interc$conf.low
  high  <- interc$conf.high
  pval  <- interc$p.value
  pval   <- summary(model_lm)$coefficients["(Intercept)", "Pr(>|t|)"]
  star <- if (pval < 0.001) {
    "***"
  } else if (pval < 0.01) {
    "**"
  } else if (pval < 0.05) {
    "*"
  } else {
    ""
  }
  diff_list[[antigen_type]] <- data.frame(
    Antigen   = antigen_type,
    estimate  = est,
    conf.low  = low,
    conf.high = high,
    p.value   = pval,
    star      = star,
    stringsAsFactors = FALSE
  )
  
  tbl <- tbl_regression(model_lm, intercept = TRUE) %>% 
    add_global_p() %>% 
    bold_p(t = 0.05) %>% 
    add_n() %>%
    modify_header(label = "**Predictors**") %>%
    modify_spanning_header(everything() ~ paste0("**","IgG ",antigen_type,"**")) %>% 
    bold_labels() %>% 
    as_flex_table()
  
  library(officer)
  doc <- read_docx()
  
  doc <- doc %>%
    body_add_par("Regression Table", style = "table title")
  
  doc <- doc %>%
    body_add_flextable(tbl)
  
  print(doc, target = paste0(IgG_paired_path,"IgG ",gsub("/", "_", paste0(antigen_type)),".docx"))
}

diff_stats <- do.call(rbind, diff_list)

source("./draw_intercept.R")
effect_plot <- draw_intercept(
  diff_stats,
  palette = c("#66C2A4","#FC8D62","#8DA0CB","#E78AC3"),
  title   = "Adjusted HD–SD IgG at 28dpv",
  ylab    = "Adjusted Mean Difference (95% CI)",
)
ggsave(paste0(IgG_paired_path,"IgG_adjusted_HD-SD.png"), plot = effect_plot, width = 6, height = 4.67, dpi = 300)

effect_plot
