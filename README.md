<img align="left" width="160" height="160" src="https://github.com/peterson-tim-j/HydroSight/blob/master/GUI/icons/icon_webpage.png">  

# _HydroSight_: _Open-source data-driven hydrogeological insights_
 
[![GitHub release](https://img.shields.io/github/release/peterson-tim-j/HydroSight)](https://github.com/peterson-tim-j/HydroSight/releases/) [![Github All Releases](https://img.shields.io/github/downloads/peterson-tim-j/HydroSight/total.svg?style=flat)]()   [![GitHub license](https://img.shields.io/github/license/peterson-tim-j/HydroSight)](https://github.com/peterson-tim-j/HydroSight/blob/master/LICENSE) [![GitHub forks](https://img.shields.io/github/forks/peterson-tim-j/HydroSight)](https://github.com/peterson-tim-j/HydroSight/network) [![GitHub stars](https://img.shields.io/github/stars/peterson-tim-j/HydroSight)](https://github.com/peterson-tim-j/HydroSight/stargazers)

HydroSight is a highly flexible statistical object-oriented toolbox for deriving quantitative insights value from groundwater monitoring data. It includes a powerful graphical interface for analying 100s of groundwater observation bores, allowing data-driven:

* decomposition of a groundwater hydrograph into individual drivers, such as climate and pumping or climate ([Shapoori et al., 2015a](https://github.com/peterson-tim-j/HydroSight/blob/master/documentation/html/papers/Shapoori_2015A.pdf)) and landuse change ([Peterson and Western, 2014](https://doi.org/10.1029/2017WR021838)).
* estimation of aquifer hydraulic properties from the hydrograph ([Shapoori et al., 2015c](https://github.com/peterson-tim-j/HydroSight/blob/master/documentation/html/papers/Shapoori_2015C.pdf), [Peterson et al., in-press](https://doi.org/10.1111/gwat.12946))
* statistical identification of the major groundwater processes ([Shapoori et al., 2015b](https://github.com/peterson-tim-j/HydroSight/blob/master/documentation/html/papers/Shapoori_2015B.pdf)).
* interpolation or extrapolation of the observed hydrographs to a regular time step ([Peterson & Western, 2018](https://doi.org/10.1029/2017WR021838).
* simulation of groundwater head under different climate or, say, pumping scenarios.
* numerical identification of hydrograph monitoring errors and outliers ([Peterson et al., 2018](https://doi.org/10.1007/s10040-017-1660-7)).

## Installation

_HydroSight_ has multiple installation options:
1. Stand-alone app within Windows. The latest .exe is [available here](https://github.com/peterson-tim-j/HydroSight/releases).
1. Load _Hydrosight_ from within Matlab by downloading the latest source code [here](https://github.com/peterson-tim-j/HydroSight/releases).
1. Compile your own stand-alone app from within Matlab bu (i) downloading the [source code](https://github.com/peterson-tim-j/HydroSight/releases) and (ii) running the command: makeStandaloneHydroSight()

For futher details see the [installation wiki page](https://github.com/peterson-tim-j/HydroSight/wiki).

## Examples
Multiple examples are built into the _HydroSight_ GUI, each highlighting aspects of the above papers. Soon, each example will also be supported by online videos. In the meantime major aspects of the graphical interface and the algorithms are outlined on the [wiki page](https://github.com/peterson-tim-j/HydroSight/wiki).

_HydroSight_ can also be run from the Matlab command window. For an example of this [see here](https://github.com/peterson-tim-j/HydroSight/blob/master/algorithms/models/TransferNoise/Example_model/example_TFN_model.m).

## What does _HydroSight_ look like?

The _HydroSight_ interface includes tabs for each step in the modelling of groundwater hydrographs:
1. Project documentation
2. Hydrograph outlier detection
1. Time-series model construction, specifically defining the data and form of the model.
1. Model calibration and tools to examine the internal dynamics of the calibrated model. The screenshot below shows this tab and an estimate of the annual groundwater recharge.  
1. Model simulations, allowing hydrograph decomposition, exploration of scenarios (e.g. different climate or pumping), hindcasting and interpolation 

![_HydroSight_ Recharge estimation](https://user-images.githubusercontent.com/8623994/190363849-d6e8f457-7891-4213-8ace-71076e69e4f6.png)

## Contributing

_HydroSight_ is an ongoing research project and its funding depends upon evidence of impact. To help us demonstrate impact please ensure you cite the relevant papers (using the "Cite Project" option within the GUI). 

And, if _HydroSight_ doesn't do what you need then please consider (i) writing and sharing your own module by extending existing class definitions, (ii) or getting in touch to discuss collaborations.