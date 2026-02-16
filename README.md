# High-Dose vs Standard-Dose Influenza Vaccine Analysis in HIV

Statistical analysis code for the article: **"Limited effectiveness of high-dose flu vaccine in augmenting influenza A responses in older people with HIV"**

## Citation

Kupritz J*, Davis S*, Singh P, Liu T, Diaz-Pachon D, Rodriguez A, Pahwa R, Pallikkuth S, Pahwa S. Limited effectiveness of high-dose flu vaccine in augmenting influenza A responses in older people with HIV. University of Miami Miller School of Medicine.

**Clinical Trial**: [NCT04487041](https://clinicaltrials.gov/study/NCT04487041)

## Study Overview

This study examines the immunogenicity of high-dose (HD) vs standard-dose (SD) influenza vaccines in people with HIV (PWH) compared to people without HIV (PWoH), across different age groups.

### Study Population (N=235)

| Group | Age Range | PWH | PWoH |
|-------|-----------|-----|------|
| Young | 18-40 years | 42 | 56 |
| Old | ≥60 years | 67 | 71 |

### Vaccine Administration
- **Standard-dose (SD)**: Given during 2020-21, 2021-22, or 2022-23 flu seasons
- **High-dose (HD)**: Given the following season to a subset of participants (HD contains 4x more antigen than SD)

### Measurements
- **HAI titers**: Hemagglutination inhibition titers
- **HA-specific IgG**: Measured using Luminex-based assay
- **Timepoints**: 0, 7, 14, 28, and 180 days post-vaccination (dpv)

### Antigens Analyzed
- A/H1N1 (A/H1)
- A/H3N2 (A/H3)
- B/Victoria (B/Vic)
- B/Yamagata (B/Yam)

## Project Structure

```
multi_codes/
├── multi_analysis.R              # Main analysis script
├── contrast_for_ordinalrg.R      # Contrast plot generation function
├── draw_intercept.R              # HD-SD difference visualization function
├── serology_master_dataset.xlsx  # Source data
├── multi_codes.Rproj             # RStudio project file
└── outputs/
    ├── HAI/                      # HAI titer regression results by dose
    ├── IgG/                      # IgG regression results by dose
    ├── HAI_paired/               # HD vs SD HAI comparison (paired analysis)
    ├── IgG_paired/               # HD vs SD IgG comparison (paired analysis)
    ├── SD_contrast_plot.png      # Combined SD contrast plots
    └── HD_contrast_plot.png      # Combined HD contrast plots
```

## Statistical Methods

### Multivariate Analysis by Dose

**HAI Titer Analysis**
- **Model**: Ordinal logistic regression (`polr` from MASS package)
- **Outcome**: HAI titer at 28 dpv (T3)
- **Covariates**: Group, Sex, Race, Ethnicity, Previous flu vaccination, Baseline titer (T0)

**IgG Analysis**
- **Model**: Linear regression
- **Outcome**: HA-specific IgG at 28 dpv (T3)
- **Covariates**: Group, Sex, Race, Ethnicity, Previous flu vaccination, Baseline level (T0)

### HD vs SD Paired Comparison
- **Model**: Linear regression on paired differences (HD_T3 - SD_T3)
- **Outcome**: Difference in post-vaccination response between HD and SD
- **Interpretation**: Intercept represents the adjusted mean difference (HD - SD), controlling for covariates


## Dependencies

```r
install.packages(c(
  "tidyverse",
  "readxl",
  "MASS",
  "emmeans",
  "flextable",
  "gtsummary",
  "officer",
  "patchwork",
  "car",
  "broom",
  "RColorBrewer"
))
```

## Usage

1. Open the project in RStudio: `multi_codes.Rproj`
2. Ensure `serology_master_dataset.xlsx` is in the project root
3. Run `multi_analysis.R` to execute the full analysis pipeline

## Outputs

- **Word Documents (.docx)**: Regression tables with odds ratios (HAI) or coefficients (IgG), 95% CI, and p-values
- **Contrast Plots (.png)**: Estimated marginal means with 95% CI and significance annotations
- **HD-SD Difference Plots**: Adjusted mean differences between HD and SD with 95% CI

## Statistical Notes

- Significance: * p < 0.05, ** p < 0.01, *** p < 0.001
- HAI results presented as odds ratios (exponentiated coefficients from ordinal regression)
- IgG results presented as beta coefficients from linear regression
- Emmeans used for pairwise group comparisons

## Contact

Correspondence: spallikkuth@med.miami.edu; spahwa@med.miami.edu

Department of Microbiology and Immunology
University of Miami Miller School of Medicine
Miami, Florida, United States of America

