# QGplasticity_S_purpuratus
Version controlled and editable source for the data and code supporting the paper by Strader et al. on _Strongylocentrotus purpuratus_ quantitative genetics of plasticity.

## Data
### Data citation
If you use the data, please cite the data package TBD<!-- archived in the Dryad Data Repository: https://doi.org/ -->

### Data metadata
Description of variables and abbreviations used in the data for the manuscript "Genetic variation underlies responses to global change drivers in the purple sea urchin, _Strongylocentrotus purpuratus_." The data and code provided here are sufficient to replicate the analyses presented in the above work.

The main files contain information related to either one of two parental environment treatments: non-upwelling (abbreviated `parN`) or upwelling (abbreviated `parU`) conditions. Filenames include either `parN` or `parU` to indicate what is contained within.

Column headings for the datasets reflect variables defined and discussed in the main manuscript and supplementary materials. In brief, these are:

  - `Length.spi` and `Length.bod`: measurements (in mm) of larval spicule length and body length, respectively.
  - `scLength.spi` and `scLength.bod`: scaled values of the above two measurements. These are obtained by multiplying the original measurement by 100.
  - `treat_dev` and `treat_adult`: Either non-upwelling (N) or upwelling (U) environmental conditions in which each larvae was reared or each parent was conditioned.
  - `measurer.spi` and `measurer.bod`: identity of the individual conducting the measurement, coded as -0.5 or 0.5. Two people total made all of the measurements.
  - `animal`: unique code for each individual that was measured.
  - `cultureFac` (sometimes just `culture`): factor indicating the culture tube in which an individual larvae was reared.
  - `blockFac` (sometimes just `block`): factor indicating the experimental block (day) in which an individual larvae was produced.
  - `Dam` and `Sire`: each individual's maternal or paternal identity.
  - `cross`: unique factor indicating the cross between maternal egg and paternal sperm that created each individual.
  - `bucket`: factor indicating the rearing bucket in which all eggs from a given cross were reared.
  - `PercAb`: proportion (0-1) of eggs assayed that were identified as displaying abnormal growth or development.
  - `Average`: average egg diameter for a dam. 
  - `tube`: factor indicating culture tube associated with dam identity from which eggs were sampled.

## Changes
For ease of reference, an overview of significant changes to be noted below. Tag with commits or issues, where appropriate.

### Major

### Minor
