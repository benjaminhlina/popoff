# popr: Fit Hidden Markov models to electronic pop-off satellite transmitters (PSAT) data
> [!CAUTION]
> As of August 2026, development of `popr` is ongoing. Please direct any new issues or pull requests on this repository when finding bugs. 
<!-- badges: start -->

[![R-CMD-check](https://github.com/benjaminhlina/popr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/benjaminhlina/popr/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/benjaminhlina/popr/graph/badge.svg)](https://app.codecov.io/gh/benjaminhlina/popr)
<!-- badges: end -->

Hidden Markov models for pop-off satellite transmitters (PSATs). 

## Installation

You can install the [development
version](https://github.com/benjaminhlina/popr) of {popr} using the
following:

``` r
pak::pak("benjaminhlina/popr")
```

To view the overall progress of the package please see
[news]()

## Vignettes

To follow soon

## Citations

To cite this package please cite the following publications:

To follow soon.


# Contribution and Issues Policy 
If you would like to contribute please feel free to fork the repository, create a branch, implement your changes, and submit a pull request (PR) for review. 
Important things to know when contributing is that repository uses [air]() for code formatting and styling, `snake_case` as I am dyslexic and struggle to 
write and read code that does not use `snake_case`, and uses DRY coding principle compared to WET princples.    

When submitting a PR all changes to functions need to be able to captured in unit tests, pass unit tests, and need to follow DRY prinicples.  

The package has used LLM as a tool mostly in for creating helpful error messages and unit tests.

This repository will not be accepting LLM submitted PRs nor PRs that contain directly copy and pasted code from LLMs that has not be vetted and edited. 

Issues that are submitted by LLMs and/or copied and pasted from LLMs will be deleted and ignored, as the issuer needs to have the ability to be able to write their own 
issues and understand the code that is producing the issue. 
