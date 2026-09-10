# popr: Fit hidden Markov models (HMM) to electronic pop-off satellite transmitter (PSAT) data
> [!CAUTION]
> As of August 2026, development of `{popr}` is ongoing. Please direct any new issues or pull requests on this repository when finding bugs. Even though it uses unit tests, plaese test it 
> and see what does and doesn't work.  
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

To view the overall progress of the package please see [news](https://github.com/benjaminhlina/popr/blob/main/NEWS.md)

## Vignettes

To follow soon

## Citations

To cite this package please cite the following publications:

- Gatti, P., Robert, D., Fisher, J.A.D., Marshall, R.C., and Le Bris, A. 2020. Stock-scale electronic tracking of Atlantic halibut reveals summer site fidelity and winter mixing on common spawning grounds. ICES J Mar Sci 77(7–8): 2890–2904. [doi:10.1093/icesjms/fsaa162](https://doi.org/10.1093/icesjms/fsaa162).

- Le Bris, A., Fisher, J.A.D., Murphy, H.M., Galbraith, P.S., Castonguay, M., Loher, T., and Robert, D. 2018. Migration patterns and putative spawning habitats of Atlantic halibut (<i>Hippoglossus hippoglossus</i>) in the Gulf of St. Lawrence revealed by geolocation of pop-up satellite archival tags. ICES J Mar Sci 75(1): 135–147. [doi:10.1093/icesjms/fsx098](https://doi.org/10.1093/icesjms/fsx098).

- Le Bris, A., Fréchet, A., and Wroblewski, J.S. 2013. Supplementing electronic tagging with conventional tagging to redesign fishery closed areas. Fisheries Research 148: 106–116. [doi:10.1016/j.fishres.2013.08.013](https://doi.org/10.1016/j.fishres.2013.08.013).

- Pedersen, M.W., Righton, D., Thygesen, U.H., Andersen, K.H., and Madsen, H. 2008. Geolocation of North Sea cod (<i>Gadus morhua</i>) using hidden Markov models and behavioural switching. Can. J. Fish. Aquat. Sci. 65(11): 2367–2377. NRC Research Press. [doi:10.1139/F08-144](https://doi.org/10.1139/F08-144).

- Ransier, K.T., Gatti, P., Le Bris, A., Den Heyer, C.E., Claireaux, G., Wringe, B., and Fisher, J.A.D. 2024. Electronic tags reveal high migratory diversity within the largest Atlantic halibut (<i>Hippoglossus hippoglossus</i>) stock. Can. J. Fish. Aquat. Sci. 81(7): 828–846. [doi:10.1139/cjfas-2023-0282](https://doi.org/10.1139/cjfas-2023-0282).

- Thygesen, U.H., Pedersen, M.W., and Madsen, H. 2009. Geolocating Fish Using Hidden Markov Models and Data Storage Tags. In Tagging and Tracking of Marine Animals with Electronic Devices. Edited by J.L. Nielsen, H. Arrizabalaga, N. Fragoso, A. Hobday, M. Lutcavage, and J. Sibert. Springer Netherlands, Dordrecht. pp. 277–293. [doi:10.1007/978-1-4020-9640-2_17](https://link.springer.com/chapter/10.1007/978-1-4020-9640-2_17).


# Contribution and Issues Policy 

If you would like to contribute please feel free to fork the repository, create a branch, implement your changes, and submit a pull request (PR) for review. 

Important things to know when contributing are the following: 

1. The repository uses [air](https://posit-dev.github.io/air/) for code formatting and styling. Please see air's documentation to setup on your machine if you are contributing. 

2. The repository uses `snake_case` as I am dyslexic and struggle to write and read code that does not use `snake_case`. 

3. The repository uses DRY coding principle compared to WET princples.    

4. When submitting a PR all changes to functions need to be able to captured in unit tests, pass unit tests and/or have additonal unit tests, and need to follow DRY prinicples.  

5. The package has used LLM as a tool in for creating helpful error messages and unit tests.

6. This repository will not be accepting LLM submitted PRs nor PRs that contain directly copy and pasted code from LLMs that has not be vetted and edited. 

7. Issues that are submitted by LLMs and/or copied and pasted from LLMs will be deleted and ignored, as the issuer needs to have the ability to be able to write their own 
issues and understand the code that is producing the issue. Preferably the issues use [{reprex}](https://reprex.tidyverse.org/).  
