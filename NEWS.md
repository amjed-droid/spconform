# spconform 0.1.0

* Initial release of `spconform` on CRAN.
* Implemented locally weighted split conformal prediction for geostatistical (point-referenced) data (`scp_geostatistical()`).
* Implemented neighbourhood-weighted leave-one-out conformal prediction for areal (lattice) data (`scp_areal()`).
* Added kernel weighting utilities for spatial and spatio-temporal predictions (`spatial_kernel_weights()`, `areal_neighbor_weights()`).
* Implemented comprehensive spatial diagnostic tools and multi-panel visualization (`diagnose()`).
* Included S3 methods for printing, summarizing, and plotting conformal prediction intervals (`print()`, `summary()`, `plot()`, `coverage_report()`).
