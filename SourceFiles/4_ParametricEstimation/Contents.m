% Parametric Estimation Module
% Frequency Domain Identification Toolbox
% Version 2.0 (R2013a) Feb-2014 
% 
% Deterministic Estimators
%   lsfdi       - Linear Least Squares analytical estimator
%                 adviced for starting values calculations.
%   wlsfdi      - Weighted Least Squares analytical estimator
%                 adviced for improved starting values calc.
%   nlsfdi      - Non-linear Least Squares numerical estimator
%                 requires only measured input-output freq. data.
%
% Stochastic Estimators
%   mlfdi       - Maximum-Likelihood stochastic estimator
%                 requires measured non-parametric noise model.
%   gtlsfdi     - Generalized Total Least Squares estimator
%                 adviced for starting values calculations.
%   btlsfdi     - Bootstrapped Total Least Squares estimator
%                 numerical robust version of maximum-likelihood.
%
% Sub-Space Estimators
%   ssfdi       - Subspace Identification for State-Space models.
% 
% Auxiliary Functions
%   mlfdi_res   - Maximum Likelihood Estimation Residuals.
%   nllsfdi_res - Non Linear Least Squares FDI Residuals.
%   qsvd        - Quadratic Singular Value Decomposition.
%
% Author: Thomas Beauduin
% Copyright (c) PMA, KULeuven, Belgium