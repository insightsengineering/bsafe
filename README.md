# B-SAFE Shiny Module

<!-- start badges -->

[![Check 🛠](https://github.com/insightsengineering/bsafe/actions/workflows/check.yaml/badge.svg)](https://insightsengineering.github.io/bsafe/)
[![Docs 📚](https://github.com/insightsengineering/bsafe/actions/workflows/docs.yaml/badge.svg)](https://insightsengineering.github.io/bsafe/main-unit-test-report)
[![Code Coverage 📔](https://raw.githubusercontent.com/insightsengineering/bsafe/_xml_coverage_reports/data/main/badge.svg)](https://insightsengineering.github.io/bsafe/main/coverage-report/)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Current Version](https://img.shields.io/github/r-package/v/insightsengineering/bsafe/main?color=purple\&label=package%20version)](https://github.com/insightsengineering/bsafe/tree/main)
<!-- end badges -->

## Installation

``` r
install.packages("remotes")
install.packages("teal")
remotes::install_github("https://github.com/insightsengineering/bsafe")
remotes::install_github("https://github.com/insightsengineering/teal.modules.bsafe")
```

## Usage

``` r
library(teal)
library(bsafe)
library(teal.modules.bsafe)


data <- teal.modules.bsafe::test_data

data <- teal.data::teal_data(
  bsafe_data = data,
  code = expression({
    bsafe_data <- teal.modules.bsafe::test_data
  })
) |>
  teal.data::verify()

app <- teal::init(
  data = data,
  modules = list(
    teal.modules.bsafe:::tm_bsafe(
      label = "teal.modules.bsafe",
      dataset_name = "bsafe_data"
    )
  ),
  header = "DaVinci test of bsafe as teal module"
)
shiny::shinyApp(app$ui, app$server)
```

## Further Information

B-SAFE is an R-Shiny app. The app is an innovate software tool for statistical analysis of adverse event summary data. The app can enhance the descriptive analysis for a current trial with historical information on one or more treatment arms for increased precision. It features a Bayesian Meta-Analytic Predictive (MAP) Prior approach and a robust extension, which incorporates historical information for safety analyses on adverse events into safety analyses for a new trial. The use of historical information has been used for efficacy analyses in the past and now being extended to safety analyses.

[Manual](https://github.com/insightsengineering/teal.modules.bsafe/blob/main/inst/www/manual.pdf)
[Stats Plan](https://github.com/insightsengineering/teal.modules.bsafe/blob/main/inst/www/statsplan.pdf)
[Test plan](https://github.com/insightsengineering/teal.modules.bsafe/blob/main/inst/www/test_validate_review.pdf)
