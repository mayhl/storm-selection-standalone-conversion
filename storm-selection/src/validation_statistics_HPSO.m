function [ES_tot, ES_AEP, tgtThreshs_rec] = validation_statistics_HPSO(x, stmRate, Resp, HC_target, opts, lp_model)
% This function provides various error metrics for evaluating HPSO
% performance.
% Input
%  - x: a binary vector representing storm selection in which 1 represent the corresponding storm is selected.
%  - stmRate: a column vector of rates for selected storms
%  - Resp: a matrix of surge responses for storms in the original suite (Ns x Nn matrix;
%          Ns is number of storms in original database not number of selected storms, 
%          Nn is number of nodes)
%  - HC_target: a cell of size Nn with each element corresponding to 
%               a structure having information of target hazard curve. If
%               empty, target hazards curves are estimated using surge
%               responses 'Resp' and probability masses for storms in the original suite
%               'ProbMass'. 
%               Each structure has following fields,
%               .x:  surge levels
%               .Pf: probabilities of exceedance
%               ** 'HC_target' is required to create the structure 'lp_model'
%  - opts: a structure with options for the storm selection, obtained after
%          running 'StormRatesAdjustment'
%          ** 'opts' is required to create the structure 'lp_model'
%  - lp_model: [OPTIONAL] a structure with model for linear programming, obtained
%              after running 'StormRatesAdjustment'.
%              ** If not exist, structures 'opts' and 'HC_target' are required 
%              to create 'lp_model'.
%
% Output
%  - ES_tot: a structure with global error statistics calculated across all locations 
%            and all target AEPs with following fields
%            .devP: weighted average of absolute deviation of estimated AEP levels 
%                   from target AEP levels 
%            .corP: (weighted) correlation coefficient of estimated surge thresholds
%                   to true thresholds
%            .MAE: average of absolute deviation of estimated
%                  surge thresholds from true thresholds
%            .MAPE: average of percentage of absolute deviation of
%                   estimated surge thresholds from true thresholds
%  - ES_AEP: a structure with error statistics per each target AEP level
%            calculated across all locations with following fields
%            .devP: weighted average of absolute deviation of estimated AEP levels 
%                   from target AEP levels 
%            .corP: (weighted) correlation coefficient of estimated surge thresholds
%                   to true thresholds
%            .MAE: median of absolute deviation of estimated surge thresholds 
%                  from true thresholds
%            .MAPE: median of percentage of absolute deviation of estimated 
%                   surge thresholds from true thresholds
%  - tgtThreshs_rec: a matrix of surge thresholds obtained by selected
%                    storms and their adjusted rates (Np x Nn matrix;
%                    Np is number of target AEP levels
%
% V1.0 10/15/2024  WoongHee Jung (wjung2@nd.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5 % Need at least 5 input arguments
    error('Not enough input arguments!')
end

Nn = size(Resp, 2); % number of locations
Nt = length(opts.tgtRate); % number of target AEPs

if ~iscolumn(stmRate) % make sure that 'stmRate' is a column vector
    stmRate = stmRate';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if 'lp_model' is provided, use as it is, otherwise create structure 'lp_model' with 
% only necessary fields based on given information (brought from 'StormRatesAdjustment.m')
if ~exist('lp_model', 'var') 
    %-- Find surge thresholds corresponding to target surge rates  
    % Find anchor points for nodes whose surge responses do not cover the
    % target rates. Anchor points can be either left or right end of the target 
    % hazard curve or both ends 
    % initialization
    tgtThreshs = zeros(Nt, Nn); tgtRates = tgtThreshs; % matrices with each column being target thresholds or rates for each node of interest
    NoThresh = tgtThreshs; % Number of responses over thresholds; to construct P matrix (sparse matrix) in the next section
    Ntsr = zeros(Nn, 1);   % Number of target surge rates for each node of interest
    % re-define variables for using parfor (to avoid large broadcast variables)
    tgtRate = opts.tgtRate;
    interp  = opts.interp;
    HC_target_given = opts.HC_target_given;
    parfor i = 1:Nn
        [tgtRates(:, i), tgtThreshs(:, i), NoThresh(:, i)] = TargetRatesNThreshs_LargeA_mod(Resp(:, i), HC_target{i}.Pf, HC_target{i}.x, tgtRate, interp, HC_target_given);
        Ntsr(i) = sum(~isnan(tgtRates(:, i)));
    end
    Rates = reshape(tgtRates, 1, []); % concatenate target rates for calculating absolute deviations
    Rates(isnan(Rates)) = []; % Remove NaNs
    lp_model.tgtThreshs = tgtThreshs;
    lp_model.anchorpts  = (abs(tgtRates - repmat(tgtRate', 1, Nn)) > eps);
    
    %-- Construct LARGE Indicator matrix P (of size [sum(Ntsr) x Ns])
    % Collect indices in P matrix where 1 (surge response being greater than surge threshold) is present.
    stormi  = zeros(sum(NoThresh, 'all'), 1);
    targeti = stormi;
    tempi   = 1;
    for i = 1:Nn
        ind_temp = find(~isnan(tgtThreshs(:, i)));
        for j = 1:Ntsr(i)
            temp = find(Resp(:, i) >= tgtThreshs(ind_temp(j), i));
            stormi(tempi:tempi+length(temp)-1)  = temp;
            targeti(tempi:tempi+length(temp)-1) = sum(Ntsr(1:i))-Ntsr(i)+j;
            
            tempi = tempi + length(temp);
        end
    end
    % Construct a sparse matrix for P
    %%%%%%%%%%%%%%%%%%% i index, j index, elements,              ,  i size  ,    j size
    lp_model.P = sparse(targeti,  stormi, ones(length(stormi), 1), sum(Ntsr),      Ns);

    %-- Save other variables in lp_model for linear programming
    lp_model.rhs  = [zeros(1, opts.nselect+length(opts.storms_included)), Rates, -Rates]';  % right hand side of the constraint
    lp_model.Ntsr = Ntsr;  % Number of target surge rates for each node of interest
    lp_model.tgtRates = tgtRates; % matrix with each column being target rates for each node of interest
    
    % Update all weights
    if strcmp(opts.cluster.perform, 'yes')
        weight_nodes = repelem(opts.cluster.weight_nodes, Ntsr);
    else
        weight_nodes = repelem(opts.weight_nodes(opts.SS_nodes_use), Ntsr);
    end

    if strcmp(opts.error, 'rel')
        % weights inversly proportional to target rates
        weight_rates = interp1(1./opts.tgtRate, opts.weight_rates./opts.tgtRate*opts.tgtRate(1), 1./Rates);
        weight_rates_unorm = interp1(1./opts.tgtRate, opts.weight_rates, 1./Rates);
    else
        % weights proportional to target rates
        weight_rates = interp1(1./opts.tgtRate, opts.weight_rates, 1./Rates);
        weight_rates_unorm = weight_rates;
    end

    lp_model.weights_nodes = weight_nodes;
    lp_model.weights_rates = weight_rates;
    lp_model.weights_rates_unorm = weight_rates_unorm;
    
    if strcmp(opts.cluster.perform, 'yes')
        lp_model.weights_Corr_nodes = opts.cluster.weight_nodes;
    else
        lp_model.weights_Corr_nodes = opts.weight_nodes(opts.SS_nodes_use);
    end
    lp_model.weights_Corr_rates = opts.weight_rates;
else
    if ~isfield(lp_model, 'weights_Corr_nodes')
        if strcmp(opts.cluster.perform, 'yes')
            lp_model.weights_Corr_nodes = opts.cluster.weight_nodes;
        else
            lp_model.weights_Corr_nodes = opts.weight_nodes(opts.SS_nodes_use);
        end
    end
    if ~isfield(lp_model, 'weights_Corr_rates')
        lp_model.weights_Corr_rates = opts.weight_rates;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Surge thresholds obtained by selected storms and their adjusted rates
Rates = lp_model.tgtRates(:);  % concatenate target rates for calculating absolute deviations
Rates(isnan(Rates)) = [];      % remove NaNs

% Initialize to calculate hazard curves, AEP levels, surge thresholds using
% selected storms
HC_rec = cell(Nn, 1);
tgtRates_rec   = zeros(Nt, Nn);
tgtThreshs_rec = zeros(Nt, Nn);

% Compute surge thresholds using storm subset 'x' and storm rates 'stmRate'
% To avoid large broadcasting variables
tgtRates = lp_model.tgtRates;
interp  = opts.interp;
HC_target_given = opts.HC_target_given;
parfor ii = 1:Nn
    % Calculate hazard curves using selected storms
    Resp_temp = Resp(:, ii);
    Resp_temp(isnan(Resp_temp)) = -inf;

    Resp_temp = Resp_temp(logical(x));
    [~, sorti] = sort(Resp_temp, 'descend');
    HC_rec{ii}.x  = Resp_temp(sorti, 1)';
    HC_rec{ii}.Pf = cumsum(stmRate(sorti))';

    removei = HC_rec{ii}.x == -inf;
    HC_rec{ii}.x(removei)  = [];
    HC_rec{ii}.Pf(removei) = [];
    
    % interpolation for surge thresholds
    Resp_temp(Resp_temp == -inf) = NaN;
    [tgtRates_rec(:, ii), tgtThreshs_rec(:, ii), ~] = TargetRatesNThreshs_LargeA_mod(Resp_temp, HC_rec{ii}.Pf, HC_rec{ii}.x, tgtRates(:, ii), interp, HC_target_given);
end

%% Error statistics
disp('Calculating error statistics ...')
tic
Error_Ra = abs(tgtRates_rec   - lp_model.tgtRates); % absolute deviation of estimated AEP levels from target AEPs
Error_Th = abs(tgtThreshs_rec - lp_model.tgtThreshs); % absolute deviation of estimated surge thresholds from true (target) surge thresholds
% Error_Th(Error_Th == 0) = eps;
Error_Th_rel = Error_Th./lp_model.tgtThreshs; % absolute percentage deviation

% weighted average across all locations
ES_AEP.devP = sum(Error_Ra.*opts.weight_nodes, 2, 'omitnan')./sum(~isnan(Error_Ra).*opts.weight_nodes, 2); 
ES_AEP.corP = spatial_corr(lp_model.tgtThreshs, tgtThreshs_rec, lp_model.anchorpts, lp_model.weights_Corr_nodes);
% Geometric mean
% ES_AEP.MAE  = prod(Error_Th.^opts.weight_nodes, 2, 'omitnan').^(1./sum(~isnan(Error_Th).*opts.weight_nodes, 2));
% ES_AEP.MAPE = prod(Error_Th_rel.^opts.weight_nodes, 2, 'omitnan').^(1./sum(~isnan(Error_Th_rel).*opts.weight_nodes, 2));
% Harmonic mean
% ES_AEP.MAE  = sum(~isnan(Error_Th).*opts.weight_nodes, 2)./sum(repmat(opts.weight_nodes, Nt, 1)./Error_Th, 2, 'omitnan');
% ES_AEP.MAPE = sum(~isnan(Error_Th_rel).*opts.weight_nodes, 2)./sum(repmat(opts.weight_nodes, Nt, 1)./Error_Th_rel, 2, 'omitnan');
% Median
ES_AEP.MAE  = median(Error_Th, 2, 'omitnan');
ES_AEP.MAPE = median(Error_Th_rel, 2, 'omitnan');

% weighted average across all AEP levels
ES_tot.devP = sum((lp_model.weights_nodes.*lp_model.weights_rates/sum(lp_model.weights_nodes.*lp_model.weights_rates_unorm))*abs(lp_model.P(:, logical(x))*stmRate - Rates));
ES_tot.corP = -lp_model.weights_Corr_rates*ES_AEP.corP/sum(lp_model.weights_Corr_rates);
ES_tot.MAE  = mean(ES_AEP.MAE);
ES_tot.MAPE = mean(ES_AEP.MAPE);

t = toc

end

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tgtRates, tgtThresh, NoThresh] = TargetRatesNThreshs_LargeA_mod(Resp, Pf, Thresh, tgtRate, interp, HC_given)
% function to find target rates and corresponding to target surge
% thresholds. This function additionally identify anchor points for nodes
% whose surge responses do not cover the target rates. Anchor points can be 
% either left or right end of the target hazard curve or both ends 
%
% -- Inputs:
%  Resp:  surge response vector for storms in the original suite (Ns x 1 vector;
%         Ns is number of storms)
%
%  Pf:    a vector with probabilies of exceedance that defines a hazard curve
%
%  Thresh:  a vector with corresponding surge thresholds to 'Pf'
%
%  tgtRate: a vector with target rates of size (1 x Nt; Nt is number of targets) 
%           on the hazard curve to match
%
%  interp:  To find surge thresholds corresponding to 'tgtRate' the hazard
%           curve is interpolated. It is conducted with linear
%           interpolation in the original scale (1) or in the logscale (2).
%
%  HC_given: an Indicator that is 1 if target hazard curves (HC_target) is
%            given, else 2 if target hazard curves is calculated based on
%            Resp and ProbMass.
%            Anchor points are identified when HC_given = 2.
%
% -- Outputs:
%  tgtRates:  a vector of size (1 x Nt) with target rates. If there are anchor points, some
%             of 'tgtRate' will be updated.
%
%  tgtThresh: a vector of size (1 x Nt) with corresponding target surge thresholds to
%             'tgtRates'
%
%  NoThresh:  a vector of size (1 x Nt) with number of surge responses over 
%             target surge thresholds 'tgtThresh'
%
%  V1.0  09/01/2023 WoongHee Jung (wjung2@nd.edu)
%  V2.0  09/28/2023 WoongHee Jung (wjung2@nd.edu)
%        Modified linear interpolation to be able to address issues that 
%        can happen when there are some storms having extremely small rates 
%  V2.1  10/11/2023 WoongHee Jung (wjung2@nd.edu)
%        Added an error message when lengths of input arguments 'Pf' and
%        'Thresh' do not match.
%  V2.2  10/19/2023 WoongHee Jung (wjung2@nd.edu)
%        Modified to saturate rates of anchor points to the minimum or 
%        maximum of the target rates.
%  V2.3  02/05/2024 WoongHee Jung (wjung2@nd.edu)
%        Modified to accommodate 'tgtRate' which has some NaN elements.
%  V2.4  03/22/2024 WoongHee Jung (wjung2@nd.edu)
%        Modified to identify anchor points only when the target hazard
%        curves are calculated.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(Pf) ~= length(Thresh)
    error('The lengths of probabilities of exceedance and surge level thresholds do not match.')
end

if length(Pf) <= 2
    tgtRates = NaN(length(tgtRate), 1); tgtThresh = tgtRates;
    NoThresh = zeros(1, length(tgtRate));
else
    % Make sure that Pf is in ascend order.
    [~, ord] = sort(Pf, 'ascend');
    Pf = Pf(ord); Thresh = Thresh(ord);

    ind_int = find(~isnan(tgtRate));
    Nt      = length(ind_int);
    % initialization
    tempR = tgtRate(ind_int);
    tempT = zeros(1, Nt);
    if interp == 1 % Linear interpolation
        for ii = 1:Nt
            ind_temp  = find(Pf > tgtRate(ind_int(ii)), 1, 'first');  % find indices of two Pfs that include target Pf
            if ind_temp == 1  % extrapolation
                tempT(ii) = interp1(Pf(1:2)+(1:2)*1e-14, Thresh(1:2), tgtRate(ind_int(ii)),'linear','extrap');
            elseif isempty(ind_temp)  % extrapolation
                tempT(ii) = interp1(Pf(end-1:end)+(1:2)*1e-14, Thresh(end-1:end), tgtRate(ind_int(ii)),'linear','extrap');
            else
                tempT(ii) = interp1(Pf(ind_temp-1:ind_temp)+(1:2)*1e-14, Thresh(ind_temp-1:ind_temp), tgtRate(ind_int(ii)),'linear','extrap');
            end
        end

        % Identify anchor points
        if HC_given == 1
            [minR, maxR] = bounds(Thresh);
        else
            [minR, maxR] = bounds(Resp);
        end
        % LEFT
        ind_small = tempT < minR;
        if sum(ind_small) > 0
            tempR(ind_small) = NaN;
            tempT(ind_small) = NaN;
            tempT(sum(ind_small)) = minR;
            ind_temp = find(Thresh >= minR, 1, 'last');
            tempR(sum(ind_small)) = interp1(Thresh(ind_temp-1:ind_temp)+(1:2)*1e-14, Pf(ind_temp-1:ind_temp), minR,'linear','extrap');
            tempR(sum(ind_small)) = min(tempR(sum(ind_small)), max(tgtRate(ind_int)));
        end
        % RIGHT
        ind_large = tempT < maxR;
        if sum(ind_large) > 0
            tempR(ind_large) = NaN;
            tempT(ind_large) = NaN;
            tempT(end-sum(ind_large)+1) = maxR;
            ind_temp = find(Thresh <= maxR, 1);
            tempR(end-sum(ind_large)+1) = interp1(Thresh(ind_temp:ind_temp+1)+(1:2)*1e-14, Pf(ind_temp:ind_temp+1), maxR,'linear','extrap');
            tempR(end-sum(ind_large)+1) = max(tempR(end-sum(ind_large)+1), min(tgtRate(ind_int)));
        end
    elseif interp == 2 % interpolation in logscale
        for ii = 1:Nt
            ind_temp  = find(Pf > tgtRate(ind_int(ii)), 1, 'first');  % find indices of two Pfs that include target Pf
            if ind_temp == 1  % extrapolation
                tempT(ii) = interp1(log(Pf(1:2)+(1:2)*1e-14), Thresh(1:2), log(tgtRate(ind_int(ii))),'linear','extrap');
            elseif isempty(ind_temp)  % extrapolation
                tempT(ii) = interp1(log(Pf(end-1:end)+(1:2)*1e-14), Thresh(end-1:end), log(tgtRate(ind_int(ii))),'linear','extrap');
            else
                tempT(ii) = interp1(log(Pf(ind_temp-1:ind_temp)+(1:2)*1e-14), Thresh(ind_temp-1:ind_temp), log(tgtRate(ind_int(ii))),'linear','extrap');
            end
        end

        % Identify anchor points
        if HC_given == 1
            [minR, maxR] = bounds(Thresh);
        else
            [minR, maxR] = bounds(Resp);
        end
        % LEFT
        ind_small = tempT < minR;
        if sum(ind_small) > 0
            tempR(ind_small) = NaN;
            tempT(ind_small) = NaN;
            tempT(sum(ind_small)) = minR;
            ind_temp = find(Thresh >= minR, 1, 'last');
            tempR(sum(ind_small)) = exp(interp1(Thresh(ind_temp-1:ind_temp)+(1:2)*1e-14, log(Pf(ind_temp-1:ind_temp)), minR,'linear','extrap'));
            tempR(sum(ind_small)) = min(tempR(sum(ind_small)), max(tgtRate(ind_int)));
        end
        % RIGHT
        ind_large = tempT > maxR;
        if sum(ind_large) > 0
            tempR(ind_large) = NaN;
            tempT(ind_large) = NaN;
            tempT(end-sum(ind_large)+1) = maxR;
            ind_temp = find(Thresh <= maxR, 1);
            tempR(end-sum(ind_large)+1) = exp(interp1(Thresh(ind_temp:ind_temp+1)+(1:2)*1e-14, log(Pf(ind_temp:ind_temp+1)), maxR,'linear','extrap'));
            tempR(end-sum(ind_large)+1) = max(tempR(end-sum(ind_large)+1), min(tgtRate(ind_int)));
        end
    end
    tgtRates = NaN(length(tgtRate), 1); tgtThresh = tgtRates;
    tgtRates(ind_int) = tempR;
    tgtThresh(ind_int) = tempT;

    NoThresh = zeros(1, length(tgtRate));
    for ii = 1:Nt
        NoThresh(ind_int(ii)) = sum(Resp>=tgtThresh(ind_int(ii)));
    end
end
end

function R = spatial_corr(A, B, indmat, wgt)
index = logical((~indmat).*(~isnan(A)).*(~isnan(B)));
aux = sum(index, 2);
R = zeros(size(A, 1), 1);
parfor ii = 1:size(A, 1)
    R(ii) = (1/aux(ii))*sum(wgt(index(ii, :)).*(A(ii, index(ii, :)) - mean(A(ii, index(ii, :)), 'omitnan')).*(B(ii, index(ii, :)) - mean(B(ii, index(ii, :)), 'omitnan')), 'omitnan')./ ...
            sqrt(mean(wgt(index(ii, :)).*(A(ii, index(ii, :)) - mean(A(ii, index(ii, :)), 'omitnan')).^2, 'omitnan'))./ ...
            sqrt(mean(wgt(index(ii, :)).*(B(ii, index(ii, :)) - mean(B(ii, index(ii, :)), 'omitnan')).^2, 'omitnan'));
end
end