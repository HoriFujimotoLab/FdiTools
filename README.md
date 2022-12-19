FdiTools
========

Frequency Domain System Identification MATLAB Toolbox.

Main reference:<br>
R. Pintelon and J. Schoukens, System Identification: A Frequency Domain Approach, 2nd ed. Wiley-IEEE Press, 2012.

## Installation 
addpath `src` to MATLAB
### Requred toolbox
* MATLAB
* Control System Toolbox
* Signal Processing Toolbox
* (optional) System Identification Toolbox

# Overview
## ExcitationDesign
* Multi sine
    * Options
        * Frequency range: min, max frequency
        * Frequency resolution
        * Frequency grid: linear/quasi-logarithmic grid
        * Frequency spectrum
* Other signals
    * Chirp sine
    * Pseudo Random Binary Sequence (PRBS)
    * Swept sine

## NonparametricFRF
* Periodic excitation
    * $\hat{G}_{ML}(j\omega)$
    * Sample (co)variance estimation $\sigma_U^2, \sigma_Y^2, \sigma_{YU}^2$
    * Asymptotic variance of $G_{ML}(j\omega)$

## ParametricEstimation
* Deterministic estimators
    * Least squares, Weighted least squares, Nonlinear least squares
* Stochastic estimators
    * Maximum likelihood estimation, Bootstrapped total least squares, Generalized total least squares

# Example
Two-mass system setup

<img src="Examples/plot/twomass.jpg?raw=true" width="400">


## ExcitationDesign
<img src="Examples/plot/1_Multisine.png?raw=true" width="600">

## NonparametricFRF
<img src="Examples/plot/2_FRFest.png?raw=true" width="600">

## NonlinearDistortions
<img src="Examples/plot/3_NL.png?raw=true" width="600">

## ParametricEstimation
<img src="Examples/plot/4_deterministic.png?raw=true" width="600">


<img src="Examples/plot/4_stochastic.png?raw=true" width="600">

