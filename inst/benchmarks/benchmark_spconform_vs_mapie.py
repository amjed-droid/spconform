import numpy as np
import pandas as pd
import subprocess
import os
import matplotlib.pyplot as plt
from mapie.regression import SplitConformalRegressor, CrossConformalRegressor
from sklearn.linear_model import LinearRegression

# 1. Generate Synthetic Spatial Dataset with Heteroskedastic Spatial Variance
np.random.seed(42)
n_train = 300
n_test = 150

# Spatial coordinates in [0, 10] x [0, 10]
coords_train = np.random.uniform(0, 10, size=(n_train, 2))
coords_test = np.random.uniform(0, 10, size=(n_test, 2))

# True spatial trend
def true_trend(coords):
    return 2.0 * np.sin(coords[:, 0] * 0.5) + 1.5 * np.cos(coords[:, 1] * 0.5)

# Spatial heteroskedastic noise: std increases from 0.4 in the West (x=0) to 2.8 in the East (x=10)
sigma_train = 0.4 + 2.4 * (coords_train[:, 0] / 10.0)
sigma_test  = 0.4 + 2.4 * (coords_test[:, 0] / 10.0)

y_train = true_trend(coords_train) + np.random.normal(0, sigma_train)
y_test  = true_trend(coords_test)  + np.random.normal(0, sigma_test)

# Scratch directory for clean inter-process communication
scratch_dir = "C:/Users/intel/.gemini/antigravity-ide/brain/0c293786-e704-424d-8835-24a93e31b62e/scratch"
os.makedirs(scratch_dir, exist_ok=True)
train_csv = f"{scratch_dir}/bench_train.csv"
test_csv  = f"{scratch_dir}/bench_test.csv"
sp_out_csv = f"{scratch_dir}/spconform_bench_out.csv"
r_script_path = f"{scratch_dir}/run_spconform_bench.R"

# Save data for R script
df_train = pd.DataFrame({'x': coords_train[:, 0], 'y': coords_train[:, 1], 'resp': y_train})
df_test  = pd.DataFrame({'x': coords_test[:, 0],  'y': coords_test[:, 1],  'resp': y_test})
df_train.to_csv(train_csv, index=False)
df_test.to_csv(test_csv, index=False)

print("Dataset generated: N_train = 300, N_test = 150 (Spatial heteroskedastic noise).")

# 2. Run MAPIE 1.6.0 in Python
n_fit = int(0.5 * n_train)
X_fit, y_fit = coords_train[:n_fit], y_train[:n_fit]
X_cal, y_cal = coords_train[n_fit:], y_train[n_fit:]

# MAPIE Split Conformal
sc = SplitConformalRegressor(estimator=LinearRegression(), confidence_level=0.90, prefit=False)
sc.fit(X_fit, y_fit)
sc.conformalize(X_cal, y_cal)
_, y_pis_mapie = sc.predict_interval(coords_test)
lower_mapie = y_pis_mapie[:, 0, 0]
upper_mapie = y_pis_mapie[:, 1, 0]

# MAPIE Cross Conformal (CV+ with 5 folds)
cc = CrossConformalRegressor(estimator=LinearRegression(), confidence_level=0.90, cv=5)
cc.fit_conformalize(coords_train, y_train)
_, y_pis_cv = cc.predict_interval(coords_test)
lower_cv = y_pis_cv[:, 0, 0]
upper_cv = y_pis_cv[:, 1, 0]

# 3. Run spconform in R
r_script = f'''
library(stats)
pkg_dir <- "D:/مسابقة/جديد/spconform/spconform_v0.1.3"
r_files <- list.files(file.path(pkg_dir, "R"), full.names = TRUE)
for (f in r_files) source(f)

d_tr <- read.csv("{train_csv}")
d_te <- read.csv("{test_csv}")

s_tr <- as.matrix(d_tr[, c("x", "y")])
s_te <- as.matrix(d_te[, c("x", "y")])
y_tr <- d_tr$resp
y_te <- d_te$resp

pred_fun <- function(s_train, y_train, s_new) {{
  df_fit <- data.frame(resp = y_train, x = s_train[, 1], y = s_train[, 2])
  df_nw  <- data.frame(x = s_new[, 1], y = s_new[, 2])
  fit <- lm(resp ~ x + y, data = df_fit)
  as.numeric(predict(fit, newdata = df_nw))
}}

# Run spconform geostatistical localized conformal prediction
fit_sp <- scp_geostatistical(s_tr, y_tr, s_te, pred_fun, alpha = 0.10, split = 0.50)

out_df <- data.frame(
  pred = fit_sp$pred,
  lower = fit_sp$lower,
  upper = fit_sp$upper
)
write.csv(out_df, "{sp_out_csv}", row.names = FALSE)
'''

with open(r_script_path, "w", encoding="utf-8") as f:
    f.write(r_script)

res_r = subprocess.run([r"D:\R-4.6.1\bin\Rscript.exe", r_script_path], capture_output=True, text=True)
if res_r.returncode != 0:
    print("R Error:", res_r.stderr)
else:
    print("spconform ran successfully in R!")

sp_out = pd.read_csv(sp_out_csv)
lower_sp = sp_out['lower'].values
upper_sp = sp_out['upper'].values

# 4. Comprehensive Evaluation Metrics
def evaluate_conformal(lower, upper, y_true, coords, alpha=0.10):
    covered = (y_true >= lower) & (y_true <= upper)
    overall_cov = np.mean(covered)
    width = upper - lower
    mean_width = np.mean(width)
    
    # Winkler score: width + (2/alpha)*underage + (2/alpha)*overage
    underage = np.maximum(0, lower - y_true)
    overage  = np.maximum(0, y_true - upper)
    wis = width + (2.0 / alpha) * (underage + overage)
    mean_wis = np.mean(wis)
    
    # Spatial Stratum Analysis:
    # West Stratum (x < 5.0, Low Noise: sigma in [0.4, 1.6])
    # East Stratum (x >= 5.0, High Noise: sigma in [1.6, 2.8])
    west_mask = coords[:, 0] < 5.0
    east_mask = coords[:, 0] >= 5.0
    
    cov_west = np.mean(covered[west_mask])
    cov_east = np.mean(covered[east_mask])
    
    # Spatial conditional gap: |cov_east - 0.90| + |cov_west - 0.90|
    spatial_gap = abs(cov_east - (1.0 - alpha)) + abs(cov_west - (1.0 - alpha))
    
    return {
        'Overall Coverage (%)': overall_cov * 100,
        'West Coverage (Low Noise %)': cov_west * 100,
        'East Coverage (High Noise %)': cov_east * 100,
        'Spatial Coverage Gap': spatial_gap * 100,
        'Mean Interval Width': mean_width,
        'Mean Winkler Score': mean_wis
    }

res_mapie_split = evaluate_conformal(lower_mapie, upper_mapie, y_test, coords_test)
res_mapie_cv    = evaluate_conformal(lower_cv,    upper_cv,    y_test, coords_test)
res_spconform   = evaluate_conformal(lower_sp,    upper_sp,    y_test, coords_test)

summary_df = pd.DataFrame({
    'MAPIE (Split Conformal)': res_mapie_split,
    'MAPIE (CV+ Method)':      res_mapie_cv,
    'spconform (Spatial Localized)': res_spconform
})

print("\n================================================================================")
print("  HEAD-TO-HEAD BENCHMARK: spconform vs. MAPIE (Model Agnostic Prediction Intervals)")
print("================================================================================")
print(summary_df.round(2).to_string())
print("================================================================================")

# 5. Publication-Quality Plot
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Plot 1: Interval Width across Space (West to East)
ord_x = np.argsort(coords_test[:, 0])
x_sorted = coords_test[ord_x, 0]

axes[0].plot(x_sorted, (upper_mapie - lower_mapie)[ord_x], label='MAPIE Split (Constant)', color='#E41A1C', lw=2, ls='--')
axes[0].plot(x_sorted, (upper_cv - lower_cv)[ord_x], label='MAPIE CV+ (Constant)', color='#FF7F00', lw=2, ls=':')
axes[0].plot(x_sorted, (upper_sp - lower_sp)[ord_x], label='spconform (Spatial Adaptive)', color='#2B83BA', lw=2.5)
axes[0].set_xlabel('Spatial X Coordinate (West -> East)', fontsize=11, fontweight='bold')
axes[0].set_ylabel('Interval Width (Uncertainty)', fontsize=11, fontweight='bold')
axes[0].set_title('Spatial Adaptivity of Interval Widths', fontsize=12, fontweight='bold')
axes[0].legend(frameon=True, fontsize=9)

# Plot 2: Stratum Coverage Comparison
categories = ['West (Low Noise)', 'Overall Domain', 'East (High Noise)']
x_cat = np.arange(len(categories))
width_b = 0.25

covs_mapie = [res_mapie_split['West Coverage (Low Noise %)'], res_mapie_split['Overall Coverage (%)'], res_mapie_split['East Coverage (High Noise %)']]
covs_cv    = [res_mapie_cv['West Coverage (Low Noise %)'],    res_mapie_cv['Overall Coverage (%)'],    res_mapie_cv['East Coverage (High Noise %)']]
covs_sp    = [res_spconform['West Coverage (Low Noise %)'],   res_spconform['Overall Coverage (%)'],   res_spconform['East Coverage (High Noise %)']]

axes[1].bar(x_cat - width_b, covs_mapie, width=width_b, label='MAPIE Split', color='#E41A1C', alpha=0.85)
axes[1].bar(x_cat,           covs_cv,    width=width_b, label='MAPIE CV+',   color='#FF7F00', alpha=0.85)
axes[1].bar(x_cat + width_b, covs_sp,    width=width_b, label='spconform',   color='#2B83BA', alpha=0.85)
axes[1].axhline(90, color='black', ls='--', lw=1.5, label='Nominal Target (90%)')
axes[1].set_xticks(x_cat)
axes[1].set_xticklabels(categories, fontsize=10, fontweight='bold')
axes[1].set_ylabel('Empirical Coverage (%)', fontsize=11, fontweight='bold')
axes[1].set_title('Local Spatial Stratum Coverage', fontsize=12, fontweight='bold')
axes[1].set_ylim(60, 105)
axes[1].legend(frameon=True, fontsize=9)

# Plot 3: Winkler Interval Score (Lower is Better)
models = ['MAPIE Split', 'MAPIE CV+', 'spconform']
wis_scores = [res_mapie_split['Mean Winkler Score'], res_mapie_cv['Mean Winkler Score'], res_spconform['Mean Winkler Score']]
bars = axes[2].bar(models, wis_scores, color=['#E41A1C', '#FF7F00', '#2B83BA'], width=0.5, alpha=0.85)
axes[2].set_ylabel('Winkler Score (Lower is Better)', fontsize=11, fontweight='bold')
axes[2].set_title('Strictly Proper Scoring Rule (Winkler)', fontsize=12, fontweight='bold')
for bar in bars:
    yval = bar.get_height()
    axes[2].text(bar.get_x() + bar.get_width()/2.0, yval + 0.1, f'{yval:.2f}', ha='center', va='bottom', fontweight='bold')

plt.tight_layout()
out_fig = 'd:/مسابقة/جديد/spconform/spconform/mapie_vs_spconform_benchmark.png'
plt.savefig(out_fig, dpi=150)
downloads_fig = 'C:/Users/intel/Downloads/mapie_vs_spconform_benchmark.png'
plt.savefig(downloads_fig, dpi=150)
print(f"Saved benchmark figure to {out_fig} and {downloads_fig}")
