
# APD: Average Proportional Distance for Item Analysis from Scales

The **APD** package provides functions to compute the **Average
Proportional Distance**, a measure of internal consistency based on
pairwise proportional differences between item scores.  
This approach focuses on the *average discrepancy* between item
responses (in agreement with Sturman et al., 2009) and complements
inter-item correlation average indices.

## Features

- Compute **Average Proportional Distance (APD)** for item sets.
- Supports **bootstrap confidence intervals** for APD.
- Provides group comparisons based on confidence intervals of the
  difference (MOVER method).
- Computes average inter–item correlation and group comparison, more
  equivalent indices.

## Installation

You can install the development version from GitHub:

``` r
# From official CRAN:
install.packages("APD")
```

## Example

``` r
###### Example 1 ######
library(APD)

## Toy data: 10 persons x 5 items
set.seed(123)
dat.example1 <- matrix(sample(1:5, 50, replace = TRUE), ncol = 5)

## compute APD
APD(dat.example1, ncat = 5, ci = TRUE, level = 0.95, B = 500)

###### Example 2 ######
library(psych)

## Loading data
data("bfi")

## Choosing variables (Neuroticism factor items, more demographics)
data.bfi <- bfi[, c("N1", "N2", "N3", "N4", "N5", "gender", "age")]

## Clean for missing values
data.bfi <- data.bfi[complete.cases(data.bfi), ]


## APD for total sample
APD(data = data.bfi[, 1:5],
      ncat = 5, 
      ci = T, 
      B = 500, 
      cimethod = "perc",
      conf.level = .95)

## Item-level APD
APDitem(data = data.bfi[, 1:5], group = data.bfi$gender, 
        ncat = 5, 
        ci = T)

## Inter-item average correlation (iia) for total sample
iiacor(data = data.bfi[, 1:5])

## APD and iia for sex groups
data.bfi$gender <- as.factor(data.bfi$gender)

APDmg(data = data.bfi[, 1:5],
      ncat = 5, 
      ci = T, 
      B = 1000, 
      cimethod = "perc",
      group = data.bfi$gender, 
      conf.level = .95)

iiacor(data = data.bfi[, 1:5], group = data.bfi$gender)


## APD and iia for customized age groups
DescTools::Freq(data.bfi$age)

table(cut(data.bfi$age, breaks = c(0, 20, 30, 40, 50, 90)))

data.bfi$age4lev <- cut(data.bfi$age, breaks = c(0, 20, 30, 40, 50, 90))

APDmg(data = data.bfi[, 1:5],
      ncat = 5, 
      ci = T, 
      B = 500, 
      cimethod = "perc",
      group = data.bfi$age4lev, 
      conf.level = .95)


iiacor(data = data.bfi[, 1:5], 
       group = data.bfi$age4lev, 
       nboot = 500)
```

## Package Structure

- `APD()` – Compute Average Proportional Distance, for total group.
- `APDmg()` – Compute Average Proportional Distance for multiple groups.
- `APDitem()` – Compute item-level Average Proportional Distance (APD).
- `APDitemmg()` – Compute item-level Average Proportional Distance (APD).
- `aiicor()` – Inter–item average, total and multigroup, more
  supplementary information.
- `rmsiic()` – Root-Mean-Square Inter-Item Correlation, based on the
  squared correlation matrix.
- `aiicorEigen()` – the average inter-item association and equivalent
  first eigenvalue.

## Citation

If you use this package, please cite:

Merino-Soto, C. (2026). APD: Average Proportional Distance for Item
Analysis from Scales. R package version 1.0.1,
<https://CRAN.R-project.org/package=APD>.

You can also obtain the citation in R:

``` r
citation("APD")
```

## References

Sturman, D., Cribbie, R. A., & Flett, G. L. (2009).  
The average distance between item values: A novel approach for
estimating internal consistency.  
*Educational and Psychological Measurement*, 69(6), 913–932.
<https://doi.org/10.1177/0734282908330937>

Briggs, S.R. and Cheek, J.M. (1986).  
The role of factor analysis in the development and evaluation of
personality scales.  
*Journal of Personality*, 54, 106–148.
<https://doi.org/10.1111/j.1467-6494.1986.tb00391.x>

Clark, L. A., & Watson, D. (1995).  
Constructing validity: Basic issues in objective scale development.  
*Psychological Assessment*, 7(3), 309–319.
<https://doi.org/10.1037/1040-3590.7.3.309>

Piedmont, R.L. (2014). Inter-item correlations.  
In A.C. Michalos (Ed.), *Encyclopedia of Quality of Life and Well-Being
Research*.  
Springer, Dordrecht. <https://doi.org/10.1007/978-94-007-0753-5_1493>

## License

This package is released under the **MIT License**.
