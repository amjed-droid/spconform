---
output: pdf_document
fontsize: 11pt
geometry: margin=1in
---

\thispagestyle{empty}
\today

\noindent
To the Editors-in-Chief\\
\textbf{The R Journal}

\bigskip

\noindent
Dear Editors,

\bigskip

\noindent
Please consider our original research and software article titled:

\begin{center}
\textbf{spconform: Spatially Localized Conformal Prediction and Spatial Diagnostic Auditing in R}
\end{center}

\noindent
for publication in \textit{The R Journal} under the \textbf{Packages} and \textbf{Comparisons and Benchmarking} categories.

\paragraph{Context and Motivation:}
Quantifying predictive uncertainty is foundational across environmental monitoring, spatial epidemiology, and geosciences. Traditional spatial interpolation methods (e.g., Gaussian kriging) depend heavily on correctly specified spatial covariance structures; when spatial stationarity or covariance assumptions fail, empirical coverage collapses catastrophically (dropping to 73.4\% in our controlled simulations). Standard conformal prediction guarantees finite-sample validity without parametric distribution assumptions, but strictly requires data exchangeability---an assumption systematically violated in the presence of spatial autocorrelation and heteroskedasticity.

The \CRANpkg{spconform} package resolves this challenge by providing a rigorous, model-agnostic, localized conformal calibration framework for spatially and spatio-temporally dependent data on CRAN. 

\paragraph{Why this manuscript is relevant to readers of The R Journal:}
\begin{enumerate}
    \item \textbf{Broad Scope and General Architecture:} \CRANpkg{spconform} is strictly model-agnostic, seamlessly wrapping any point predictor in R (e.g., \CRANpkg{mgcv}, \CRANpkg{ranger}, \CRANpkg{spatial}, \CRANpkg{glm}) with distribution-free prediction intervals.
    \item \textbf{Comprehensive Multi-Modality:} Extends beyond standard geostatistical continuous regression to include areal lattice networks (\code{scp_areal()}), spatial multiclass classification with Adaptive Prediction Sets (\code{scp_classification()}), discrete epidemiology count processes (\code{scp_count()}), and joint multivariate field ellipsoids (\code{scp_multivariate()}).
    \item \textbf{Rigorous Head-to-Head Cross-System Benchmarking:} The paper provides a direct empirical comparison between \CRANpkg{spconform} and Python's leading uncertainty quantification library \textbf{MAPIE 1.6.0}, demonstrating that \CRANpkg{spconform} achieves nominal 90\% coverage across heteroskedastic spatial domains where non-spatial conformal tools fail.
    \item \textbf{Cutting-Edge Spatial Audits:} Introduces five novel diagnostic tools: Auto-Anisotropic coordinate rotation, spatial Mantel permutation tests, bivariate Cross-$L$ Ripley hotspot analysis, strictly proper scoring rules, and local efficiency-coverage trade-off surfaces (\code{conformal_ce_surface()}).
    \item \textbf{CRAN Reliability and Zero Dependencies:} Built entirely in Base R without heavy external C++ or GIS system libraries, ensuring long-term maintenance and zero installation friction across Windows, macOS, and Linux.
\end{enumerate}

All simulation scripts, empirical data, and benchmark code are fully documented, version-controlled, and guaranteed to reproduce all tables and figures down to the exact decimal.

\bigskip

\noindent
Sincerely,

\vspace{1em}

\noindent
\textbf{Ahmed Sattar Jabbar}\\
Department of Statistics, College of Administration and Economics\\
Mustansiriyah University, Baghdad, Iraq\\
Email: \href{mailto:ahmed.state.me@gmail.com}{ahmed.state.me@gmail.com}\\
CRAN Maintainer of \CRANpkg{spconform}
