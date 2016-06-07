% Nonparametric FRF Module
% Frequency Domain Identification Toolbox
% Version 2.0 (R2013a) Feb-2014 
% 
% Data Treatment
% pretreat      - transients, offsets and trends removal
%                 adviced for time-domain pretreatment of data
% splinefit     - intricate drift/trend removal by b-spline fit
%
% Data Estimators
% time2frf_h1   - Classic least-squares estimation of frf (H1)
%                 adviced for synchronized arbitrary measurements
% time2frf_ml   - Stochastic maximum-likelihood estimation of frf
%                 adviced for synchronised periodic measurements
% time2frf_log  - Non-linear logaritmic estimation of frf (Hlog)
%                 for non-synchronised or missing data measurements
%
% Author: Thomas Beauduin
% Copyright (c) PMA, KULeuven, Belgium