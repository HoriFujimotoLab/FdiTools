% Frequency Domain Identification Toolbox
% Version 8.1 (R2013a) 06-Aug-2016
%  
% Excitation Design
%   multisine       - Multisine excitation generation
%   msinl2p         - Lp-norm optimization of multisine phase
%   schroed         - Schoeder multisine phase design
%   randph          - Random multisine phase design
%   lpnorm          - Lp-norm vector calculation
%   lin2qlog        - Linear to quasi-logarithmic frequency grid
%   orthogonal      - Transformation matrix for orthogonal multisines
%   effval          - calculate the effective signal value
%   swept           - Swept-Sine excitation signal generation
%   prbs            - Pseudo-Random-Binary-Sequence signal generation
%     
% Nonparametric Estimation
%   pretreat        - transients, offsets and trends removal
%   splinefit       - intricate drift/trend removal by b-spline fit
%   time2frf_h1     - Classic least-squares estimation of frf (H1)
%                     recommended for synchronized arbitrary measurements
%   time2frf_ml     - Stochastic maximum-likelihood estimation of frf
%                     recommended for synchronised periodic measurements
%   time2frf_log    - Non-linear logaritmic estimation of frf (Hlog)
%                     for non-synchronised or missing data measurements
%
% Non-Linear Distortions
%   time2nld        - Detect the even/odd non-linear contributions
%                     with random odd-odd multisine measurements.
%   time2bla        - Reduce the main non-linear contributions
%                     with multiple random odd/full multisines.
%
% Parametric Estimation
%   lsfdi           - Linear Least Squares analytical estimator
%                     recommended for starting values calculations.
%   wlsfdi          - Weighted Least Squares analytical estimator
%                     recommended for improved starting values calc.
%   nlsfdi          - Non-linear Least Squares numerical estimator
%                     requires only measured input-output freq data.
%   mlfdi           - Maximum-Likelihood stochastic estimator
%                     requires measured non-parametric noise model.
%   gtlsfdi         - Generalized Total Least Squares estimator
%                     adviced for starting values calculations.
%   btlsfdi         - Bootstrapped Total Least Squares estimator
%                     numerical robust version of maximum-likelihood.
%   ssfdi           - Subspace Identification for State-Space models.
%
% Selection-Validation
%   chi2test        - Chi-Squares test for estimator/model set
%                     validation using Cramer-Rao lower bound.
%   costtest        - Maximum-likelihood cost function test for
%                     estimator/model set validation.
%   residtest       - Identification Residuals for validation
%
% Calculation Auxiliary
%   f2t             - Calculate time domain from fourrier coefficients
%   t2f             - Calculate frequency domain from time signals
%   dbm             - Calculate FRF magnitude in decibel (20*log10)
%   phs             - Calculate FRF phase in degree augmented with glith removal
%   theta2ba        - Transform parameter vector to rational polynomial
%   ba2theta        - Transform rational polynomial to parameter vector
%   cr_rao          - Cramer-Rao Lower Bound of parameter covariance matrix
%
% Author: Ir. Thomas Beauduin
% University of Tokyo, Hori-Fujimoto Lab