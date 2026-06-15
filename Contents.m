% Frequency Domain Identification Toolbox (FdiTools)
% Version 2.1.1 (legacy) 15-Jun-2026
%
% Legacy release. A major-upgraded successor (v3.0) is available at
% https://github.com/WataruOhnishi/FdiTools
%
% Excitation Design
%   multisine       - Multisine excitation signal generation (MIMO)
%   multisine2hdr   - Export multisine reference vector to C header file
%   prbs            - Pseudo-Random-Binary-Sequence signal generation
%   sweptsine       - Swept-Sine excitation signal generation
%
% Nonparametric Estimation
%   pretreat        - Transients, offsets and trends removal (MIMO)
%   splinefit       - Drift/trend removal by b-spline fit
%   time2frf_h1     - Classic least-squares estimation of FRF (H1)
%                     recommended for synchronized arbitrary measurements
%   time2frf_ml     - Stochastic maximum-likelihood estimation of FRF
%                     recommended for synchronised periodic measurements
%   time2frf_log    - Non-linear logarithmic estimation of FRF (Hlog)
%                     for non-synchronised or missing data measurements
%
% Non-Linear Distortions
%   time2nld        - Detect the even/odd non-linear contributions
%                     with random odd-odd multisine measurements
%   time2bla        - Best linear approximation of FRF
%                     with multiple random odd/full multisines
%
% Parametric Estimation
%   lsfdi           - Linear Least Squares analytical estimator
%                     recommended for starting values calculations
%   wlsfdi          - Weighted Least Squares analytical estimator
%                     recommended for improved starting values calc
%   nlsfdi          - Non-linear Least Squares numerical estimator
%                     requires only measured input-output freq data
%   mlfdi           - Maximum-Likelihood stochastic estimator
%                     requires measured non-parametric noise model
%   gtlsfdi         - Generalized Total Least Squares estimator
%                     adviced for starting values calculations
%   btlsfdi         - Bootstrapped Total Least Squares estimator
%                     numerical robust version of maximum-likelihood
%   ssfdi           - Subspace identification for state-space models
%                     (work in progress)
%
% Selection-Validation
%   chi2test        - Chi-Squares test for estimator/model set
%                     validation using Cramer-Rao lower bound
%   costtest        - Maximum-likelihood cost function test for
%                     estimator/model set validation
%   residtest       - Identification residuals whiteness test
%
% Calculation Auxiliary
%   f2t             - Time domain signals from Fourier coefficients
%   t2f             - Fourier serie coefficients from time signals
%   dbm             - FRF magnitude in decibel (20*log10)
%   phs             - FRF phase in degree augmented with glitch removal
%   hfrf            - Frequency response matrix of a model
%   ba2theta        - Rational polynomial to parameter vector (MIMO)
%   theta2ba        - Parameter vector to rational polynomial (MIMO)
%   ba2hm           - Parameter matrices to transfer function array
%   hm2ba           - Transfer function array to parameter matrices
%   cr_rao          - Cramer-Rao Lower Bound of parameter covariance
%   bode_fdi        - Bode plot for FdiTools data/noise objects
%   fdicohere       - Coherence for FRF object from time2frf_ml
%   fcat_fdi        - FdiTools version of fcat (concatenate freq data)
%   fdel_fdi        - FdiTools version of fdel (delete freq points)
%
% Author: Ir. Thomas Beauduin
% University of Tokyo, Hori-Fujimoto Lab
