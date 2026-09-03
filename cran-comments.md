## Resubmission Summary

Dear Leonore Hochhauser, Konstanze Lauseker, and CRAN Team,

Thank you for the review of spconform 0.1.0. We have addressed all feedback:

1. **User Filespace / Working Directory:**
   - Modified `inst/scripts/code.R` and `inst/scripts/full_reproducible_script.R` so all outputs (figures, RDS diagnostics, CSV results) strictly write to `tempdir()` via `file.path(tempdir(), "figures", ...)`.
   - Removed pre-generated CSV and RDS artifacts from `inst/`.
   - Cleaned up rendered HTML/MD secondary files from `inst/scripts/` and `vignettes/`.
   - Confirmed that no functions, examples, vignettes, tests, or scripts write to the user's home filespace or working directory.

2. **Examples / `\dontrun`:**
   - Kept self-contained, rapidly executing (< 0.1s) minimal reproducible examples in `diagnose.R`. The package contains no `\dontrun{}` tags.

3. **URL Validation:**
   - Updated the GNU license URL in `README.md` to `https://www.gnu.org/licenses/gpl-3.0.html` to avoid timeout/redirection issues during CRAN incoming checks.

## Test environments
* local Windows 11, R 4.6.1 ucrt
* Debian GNU/Linux (CRAN incoming pretest)
* win-builder (devel and release)
* R-hub v2 (Ubuntu Linux, macOS)

## R CMD check results
0 ERRORs | 0 WARNINGs | 1 NOTEs

### Explanation of NOTE:
* **New submission:** This is a first-time submission to CRAN.
* **Possibly misspelled words in DESCRIPTION:** The words ('Spatio', 'spatio', 'exchangeability', 'geostatistical') are standard domain-specific spatial statistics terms and are correctly spelled.
* **URL timeout:** Handled by pointing to the direct `.html` endpoint.