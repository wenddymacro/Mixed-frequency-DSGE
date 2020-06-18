function [best_dlnL, best_mphi, best_msigma] = dsge_likelihood_my(xparam1,DynareDataset,DatasetInfo,DynareOptions,Model,EstimatedParameters,BayesInfo,BoundsInfo,DynareResults,derivatives_info)
% Evaluates the posterior kernel of a dsge model using the specified
% kalman_algo; the resulting posterior includes the 2*pi constant of the
% likelihood function

%@info:
%! @deftypefn {Function File} {[@var{fval},@var{exit_flag},@var{ys},@var{trend_coeff},@var{info},@var{Model},@var{DynareOptions},@var{BayesInfo},@var{DynareResults},@var{DLIK},@var{AHess}] =} dsge_likelihood (@var{xparam1},@var{DynareDataset},@var{DynareOptions},@var{Model},@var{EstimatedParameters},@var{BayesInfo},@var{DynareResults},@var{derivatives_flag})
%! @anchor{dsge_likelihood}
%! @sp 1
%! Evaluates the posterior kernel of a dsge model.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item xparam1
%! Vector of doubles, current values for the estimated parameters.
%! @item DynareDataset
%! Matlab's structure describing the dataset (initialized by dynare, see @ref{dataset_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item Model
%! Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).
%! @item EstimatedParamemeters
%! Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%! @item BayesInfo
%! Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @item derivates_flag
%! Integer scalar, flag for analytical derivatives of the likelihood.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item fval
%! Double scalar, value of (minus) the likelihood.
%! @item info
%! Double vector, second entry stores penalty, first entry the error code.
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==6
%! The jacobian evaluated at the deterministic steady state is complex.
%! @item info==19
%! The steadystate routine thrown an exception (inconsistent deep parameters).
%! @item info==20
%! Cannot find the steady state, info(4) contains the sum of square residuals (of the static equations).
%! @item info==21
%! The steady state is complex, info(4) contains the sum of square of imaginary parts of the steady state.
%! @item info==22
%! The steady has NaNs.
%! @item info==23
%! M_.params has been updated in the steadystate routine and has complex valued scalars.
%! @item info==24
%! M_.params has been updated in the steadystate routine and has some NaNs.
%! @item info==26
%! M_.params has been updated in the steadystate routine and has negative/0 values in loglinear model.
%! @item info==30
%! Ergodic variance can't be computed.
%! @item info==41
%! At least one parameter is violating a lower bound condition.
%! @item info==42
%! At least one parameter is violating an upper bound condition.
%! @item info==43
%! The covariance matrix of the structural innovations is not positive definite.
%! @item info==44
%! The covariance matrix of the measurement errors is not positive definite.
%! @item info==45
%! Likelihood is not a number (NaN).
%! @item info==46
%! Likelihood is a complex valued number.
%! @item info==47
%! Posterior kernel is not a number (logged prior density is NaN)
%! @item info==48
%! Posterior kernel is a complex valued number (logged prior density is complex).
%! @end table
%! @item exit_flag
%! Integer scalar, equal to zero if the routine return with a penalty (one otherwise).
%! @item DLIK
%! Vector of doubles, score of the likelihood.
%! @item AHess
%! Matrix of doubles, asymptotic hessian matrix.
%! @item SteadyState
%! Vector of doubles, steady state level for the endogenous variables.
%! @item trend_coeff
%! Matrix of doubles, coefficients of the deterministic trend in the measurement equation.
%! @item Model
%! Matlab's structure describing the model (initialized by dynare, see @ref{M_}).
%! @item DynareOptions
%! Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%! @item BayesInfo
%! Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%! @item DynareResults
%! Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{dynare_estimation_1}, @ref{mode_check}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{dynare_resolve}, @ref{lyapunov_symm}, @ref{lyapunov_solver}, @ref{compute_Pinf_Pstar}, @ref{kalman_filter_d}, @ref{missing_observations_kalman_filter_d}, @ref{univariate_kalman_filter_d}, @ref{kalman_steady_state}, @ref{getH}, @ref{kalman_filter}, @ref{score}, @ref{AHessian}, @ref{missing_observations_kalman_filter}, @ref{univariate_kalman_filter}, @ref{priordens}
%! @end deftypefn
%@eod:

% Copyright (C) 2004-2018 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT FR

% Initialization of the returned variables and others...
fval        = [];
SteadyState = [];
trend_coeff = [];
exit_flag   = 1;
info        = zeros(4,1);
DLIK        = [];
Hess        = [];

% Ensure that xparam1 is a column vector.
xparam1 = xparam1(:);

% Set flag related to analytical derivatives.
analytic_derivation = DynareOptions.analytic_derivation;

if analytic_derivation && DynareOptions.loglinear
    error('The analytic_derivation and loglinear options are not compatible')
end

if nargout==1
    analytic_derivation=0;
end

if analytic_derivation
    kron_flag=DynareOptions.analytic_derivation_mode;
end

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if isestimation(DynareOptions) && ~isequal(DynareOptions.mode_compute,1) && any(xparam1<BoundsInfo.lb)
    k = find(xparam1<BoundsInfo.lb);
    fval = Inf;
    exit_flag = 0;
    info(1) = 41;
    info(4)= sum((BoundsInfo.lb(k)-xparam1(k)).^2);
    if analytic_derivation
        DLIK=ones(length(xparam1),1);
    end
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if isestimation(DynareOptions) && ~isequal(DynareOptions.mode_compute,1) && any(xparam1>BoundsInfo.ub)
    k = find(xparam1>BoundsInfo.ub);
    fval = Inf;
    exit_flag = 0;
    info(1) = 42;
    info(4)= sum((xparam1(k)-BoundsInfo.ub(k)).^2);
    if analytic_derivation
        DLIK=ones(length(xparam1),1);
    end
    return
end

% Get the diagonal elements of the covariance matrices for the structural innovations (Q) and the measurement error (H).
Model = set_all_parameters(xparam1,EstimatedParameters,Model);

Q = Model.Sigma_e;
H = Model.H;

% Test if Q is positive definite.
if ~issquare(Q) || EstimatedParameters.ncx || isfield(EstimatedParameters,'calibrated_covariances')
    [Q_is_positive_definite, penalty] = ispd(Q(EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness,EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness));
    if ~Q_is_positive_definite
        fval = Inf;
        exit_flag = 0;
        info(1) = 43;
        info(4) = penalty;
        return
    end
    if isfield(EstimatedParameters,'calibrated_covariances')
        correct_flag=check_consistency_covariances(Q);
        if ~correct_flag
            penalty = sum(Q(EstimatedParameters.calibrated_covariances.position).^2);
            fval = Inf;
            exit_flag = 0;
            info(1) = 71;
            info(4) = penalty;
            return
        end
    end
end

% Test if H is positive definite.
if ~issquare(H) || EstimatedParameters.ncn || isfield(EstimatedParameters,'calibrated_covariances_ME')
    [H_is_positive_definite, penalty] = ispd(H(EstimatedParameters.H_entries_to_check_for_positive_definiteness,EstimatedParameters.H_entries_to_check_for_positive_definiteness));
    if ~H_is_positive_definite
        fval = Inf;
        info(1) = 44;
        info(4) = penalty;
        exit_flag = 0;
        return
    end
    if isfield(EstimatedParameters,'calibrated_covariances_ME')
        correct_flag=check_consistency_covariances(H);
        if ~correct_flag
            penalty = sum(H(EstimatedParameters.calibrated_covariances_ME.position).^2);
            fval = Inf;
            exit_flag = 0;
            info(1) = 72;
            info(4) = penalty;
            return
        end
    end
end


%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% my code from here !!!!!!!!!!!!!!

[T, var_number] = size(DynareDataset);

for ip = 1:1 % cycle through the VAR(ip) order
   
   % initialization of the matrices in transition equation
   mphi = zeros(var_number, ip * var_number);
   msigma = eye(var_number);
   % estimation using EM algorithm
   [mphi, msigma, dlnL] = estimate_mixed_frequency(DynareDataset, ip, ...
                                                    mphi, msigma);
   
   % estimation using Quazi-Newton algorithm (as near the maximum 
   % it is better to use this algorithm)
  % [mphi, msigma, dlnL] = estimate_mixed_frequency_qn(DynareDataset, ...
  %                                                 ip, mphi, msigma, dlnL);
   
                                               
   % calculate statistics to choose the model
   dlnL = T * dlnL; % why '* T'
   ck = ip * var_number * var_number;
   dAIC = (dlnL - ck) / T; % comment: do we need to multiply by 2
   dBIC = (dlnL - ck * log(T) / 2) / T; % comment: do we need to multiply by 2, why '/T'
   
   % print the results
   disp('The estimated vectorized polynomial from the VAR(p), i.e. matrix infront of the lag of the state variable in the transition equation')
   disp(mphi)
   disp('Variance of the shocks in the transition equation')
   disp(msigma)
   fprintf('%d, %10.2f, %10.4f, %10.4f \n', ip, dlnL, dAIC, dBIC);
   
   if (ip == 1)
       best_dlnL = dlnL;
       best_mphi = mphi;
       best_msigma = msigma;
       best_dBIC = dBIC;
   else
       if (dlnL > best_dlnL)
    %  if (dBIC < best_dBIC)
           best_dlnL = dlnL;
           best_mphi = mphi;
           best_msigma = msigma;
       end
   end
   
   % comment: should I return the smoothed values also?
   % in Foroni ip = 1, here we allow for higher order VAR
   % do we need to automatically choose the model based on lnL, AIC, BIC?
    
end

