## Test environments

- Local: Windows 11, R 4.6.1
- win-builder: R-devel, Windows
- GitHub Actions: Ubuntu-latest, R-devel
- GitHub Actions: macOS-latest (Sequoia 15.7.7), R-devel
- GitHub Actions: Windows-latest, R-devel
- R-hub v2: Linux, Windows, macOS (with `donttest`)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

## Additional comments

This is the first submission of `spconform` to CRAN.

The package has been tested extensively across all major CRAN platforms 
(Linux, macOS, Windows) via both GitHub Actions and R-hub v2 to ensure 
robust, portable behaviour. All checks pass cleanly with 0 errors, 
0 warnings, and 0 notes.