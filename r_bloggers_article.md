# Why Standard Conformal Prediction Fails on Spatial Data — and How spconform Fixes It

*By Ahmed Sattar Jabbar (Maintainer and Author of `spconform`)*

------------------------------------------------------------------------

Uncertainty quantification is undergoing a revolution with **Conformal
Prediction**. Unlike traditional Bayesian or parametric methods that
rely on rigid Gaussian assumptions, conformal prediction gives you
**distribution-free prediction intervals with rigorous finite-sample
coverage guarantees** (e.g., exactly 90% confidence), regardless of your
underlying machine learning model.

If you are predicting tabular data or independent time series, split
conformal prediction works like magic.

**But what happens when your data is spatial?**

If you apply standard conformal prediction to environmental
measurements, satellite imagery, soil sensors, or regional disease maps,
**it fundamentally breaks down.**

In this post, we explain why spatial autocorrelation breaks conformal
prediction, and introduce
[**`spconform`**](https://cran.r-project.org/package=spconform) — a new
CRAN package designed specifically for **model-agnostic conformal
prediction on spatial and spatio-temporal data**.

------------------------------------------------------------------------

## 1. The Core Problem: Why Spatial Exchangeability Fails

Standard split conformal prediction assumes that your calibration data
and your test point are **exchangeable**:
``` math
(R_1, R_2, \dots, R_n, R_{\text{new}}) \sim \text{Exchangeable}
```

In spatial statistics, **Tobler’s First Law of Geography** states: \>
*“Everything is related to everything else, but near things are more
related than distant things.”*

Because of continuous spatial autocorrelation and regional
heteroscedasticity: 1. **Variance is not uniform:** Noise levels in one
region might be 5x higher than in another. 2. **Global pooling fails:**
If you pool all calibration residuals globally across the entire map,
your model will **severely under-cover in volatile regions** and
**produce unnecessarily wide, uninformative intervals in smooth
regions**.

    [Standard Conformal]  --> Pools residuals globally --> Regional Under-coverage & Bias!
    [spconform Solution] --> Weights residuals locally  --> Exact Local Conditional Coverage!

------------------------------------------------------------------------

## 2. Introducing `spconform`

The [`spconform`](https://cran.r-project.org/package=spconform) package
bridges this gap by implementing **locally weighted split conformal
prediction**: \* **Point-Referenced / Geostatistical Data:** Uses
spatial-distance Gaussian kernel weights
([`scp_geostatistical()`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md)).
\* **Areal / Lattice Data:** Uses graph-theoretic topological neighbor
weighting
([`scp_areal()`](https://amjed-droid.github.io/spconform/reference/scp_areal.md)).
\* **Spatio-Temporal Regimes:** Joint spatial and temporal lag
calibration. \* **Model-Agnostic:** Plug in Random Forests (`ranger`),
Kriging, GAMs (`mgcv`), Deep Neural Networks, or any predictor! \*
**Zero Heavy Dependencies:** Built purely on base R core packages for
maximum speed and reproducibility.

------------------------------------------------------------------------

## 3. Hands-On Example: 5 Lines of R Code

Let’s see `spconform` in action with a minimal reproducible example.

### Step 1: Install and Load `spconform` from CRAN

``` r

# Install directly from CRAN
install.packages("spconform")
library(spconform)
```

### Step 2: Simulate Non-Stationary Spatial Data

``` r

set.seed(42)
n <- 300

# 2D spatial coordinates
s <- matrix(runif(2 * n, min = 0, max = 10), ncol = 2)

# True non-linear spatial surface with heteroscedastic noise
f_true <- function(s) sin(s[, 1] / 1.5) + cos(s[, 2] / 1.5)
noise_sd <- 0.2 + 0.3 * (s[, 1] / 10)  # noise increases towards the East
y <- f_true(s) + rnorm(n, sd = noise_sd)

# Train / Test split
idx <- sample(n, 200)
s_train <- s[idx, ]; y_train <- y[idx]
s_test  <- s[-idx, ]; y_test <- y[-idx]
```

### Step 3: Define Any Machine Learning Predictor

`spconform` is model-agnostic! You only need to provide a simple wrapper
function:

``` r

# Define your favorite spatial predictor (e.g., Random Forest, GAM, or Linear Trend)
pred_fun <- function(s_train, y_train, s_new) {
  fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] + 
              I(s_train[, 1]^2) + I(s_train[, 2]^2))
  nd <- cbind(1, s_new[, 1], s_new[, 2], s_new[, 1]^2, s_new[, 2]^2)
  as.numeric(nd %*% coef(fit))
}
```

### Step 4: Run Spatial Conformal Prediction

``` r

# Generate 90% localized conformal prediction intervals
res <- scp_geostatistical(
  s_train = s_train,
  y_train = y_train,
  s0      = s_test,
  pred_fun = pred_fun,
  alpha   = 0.10,   # 90% confidence intervals
  seed    = 123
)

# Print summary
print(res)
```

**Output:**

``` text
spconform object: Locally weighted geostatistical conformal prediction
  Target miscoverage (alpha): 0.1 (nominal coverage 90%)
  Prediction locations (m):  100
  Point predictions range:   [-1.685, 1.842]
  Mean interval width:       1.412
```

------------------------------------------------------------------------

## 4. Auditing Results with the Built-in Diagnostic Suite

In spatial uncertainty quantification, you shouldn’t just trust a number
— you need to inspect spatial strata and boundary effects.

`spconform` provides a dedicated multi-panel diagnostic tool:

``` r

# Evaluate empirical coverage and visual diagnostics
diag_report <- diagnose(res, y_true = y_test, s0 = s_test)
```

This single command generates a comprehensive 4-panel diagnostic
plot: 1. **Marginal Coverage vs. Target:** Confirms overall validity. 2.
**Local Spatial Bins:** Checks that coverage is balanced across all
geographic sub-regions. 3. **Convex Hull Boundary Proximity:** Audits
edge effects and spatial extrapolation risk. 4. **Calibration QQ Plot:**
Inspects nonconformity quantile alignment.

``` r

# Get numeric coverage report
cov_info <- coverage_report(res, y_test)
cat(sprintf("Achieved Empirical Coverage: %.1f%%\n", cov_info$coverage * 100))
cat(sprintf("Average Interval Width: %.3f\n", cov_info$mean_width))
```

------------------------------------------------------------------------

## 5. What About Areal (Regional) Data?

If your data is based on counties, states, or polygon networks (e.g.,
disease counts or election statistics),
[`scp_areal()`](https://amjed-droid.github.io/spconform/reference/scp_areal.md)
utilizes topological graph adjacency matrices and breadth-first search
(BFS) neighbor weighting:

``` r

# Fast regional lattice conformal prediction
out_areal <- scp_areal(y = regional_responses, adjacency = county_adj_matrix, alpha = 0.1)
```

------------------------------------------------------------------------

## 6. Summary & Useful Links

`spconform` gives researchers, GIS analysts, and data scientists a
principled, mathematically grounded, and easy-to-use tool for
quantifying spatial prediction uncertainty without fragile assumptions.

- 📦 **CRAN Package:** <https://cran.r-project.org/package=spconform>
- 💻 **GitHub Repository:** <https://github.com/amjed-droid/spconform>
- 📖 **Documentation Site:** <https://amjed-droid.github.io/spconform/>
- 📄 **JSS Manuscript:** Submitted to the *Journal of Statistical
  Software* (September 2026).

Try it out on your spatial datasets with
`install.packages("spconform")`, and feel free to star the repo or open
an issue on GitHub!

------------------------------------------------------------------------

*Feel free to share your thoughts, spatial use cases, or feedback in the
comments below!*
