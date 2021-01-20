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

![TwoMass](https://github.com/HoriFujimotoLab/FdiTools/blob/master/Examples/plot/twomass.jpg?raw=true)


## ExcitationDesign
![Multisine](https://github.com/HoriFujimotoLab/FdiTools/blob/master/Examples/plot/1_Multisine.png?raw=true)

## NonparametricFRF
![NonparametricFRF](https://github.com/HoriFujimotoLab/FdiTools/blob/master/Examples/plot/2_FRFest.png?raw=true)

## NonlinearDistortions
![NonlinearDistortions](https://github.com/HoriFujimotoLab/FdiTools/blob/master/Examples/plot/3_NL.png?raw=true)

## ParametricEstimation
![ParametricEstimation1](https://github.com/HoriFujimotoLab/FdiTools/blob/master/Examples/plot/4_deterministic.png?raw=true)
![ParametricEstimation2](https://github.com/HoriFujimotoLab/FdiTools/blob/master/Examples/plot/4_stochastic.png?raw=true)

