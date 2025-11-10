function [ind_storm, stmRate, lp_model] = FinalAdjustment(ind_storm, stmRate, lp_model, opts, stmRate_pre)
% Function to adjust final results.
% If there are some selected storms with zero rates, they will be sequentially replaced
% with other storms so that all the selected storms will have non-zero rates.
% -- Inputs:
%  ind_storm:    a binary vector with 1 indicating 'included storm'.
%  stmRate:      adjusted rates for selected storms.
%  lp_model:   a structure with model for the original (for entire domain) linear programming 
%              with following fields.
%            .P:   sparse indicator matrix where 1 means that surge response is 
%                  greater than surge threshold
%            .A:   sparse matrix for the constraint of the linear programming
%            .rhs: right hand side of the constraint
%            .obj: vector of coefficients for the objective function
%            .Ntsr:     number of target surge rates for each node of interest
%            .tgtRates: matrix with each column being target rates for each
%                       node of interest
%            .tgtThreshs: matrix with each column being surge thresholds corresponding to
%                         target rates (lp_model.tgtRates) for each node of interest
%
%            .weights_nodes:  weights for nodes considered in 'Sdev' objective function
%            .weights_rates:  weights for rates considered in 'Sdev' objective function
%
%  opts:  a structure with options for the storm selection (detailed
%         explanation is in functions 'StormSelection' and 
%         'StormRatesAdjustment'). This function only use the following
%         fields
%            .gurobi:   Use gurobi optimizer ('yes') or not ('no').  
%                       If you want use gurobi optimizer, the optimizer 
%                       needs to be installed and a license is required. 
%
%  stmRate_pre:  preliminary adjusted rates for selected storms 
%                (ex. adjusted storm rates considering representative nodes 
%                obtained by clustering only).
%
% -- Outputs:
%  ind_storm:  updated binary vector with 1 indicating 'included storm'.
%  stmRate:    updated rates for selected storms.
%  lp_model:   a structure with model for the original (for entire domain) linear programming 
%
%  V1   10/23/2023 WoongHee Jung (wjung2@nd.edu)
%  V1.1 02/21/2024 WoongHee Jung (wjung2@nd.edu)
%       Modified to accommodate updates of fields in 'lp_model'.
%  V1.2 08/28/2024 WoongHee Jung (wjung2@nd.edu)
%       Modified to interrupt 'Final Adjustment' if it is impossible to
%       replace zero-rate storms with storms that were not selected.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns = length(ind_storm);  % Number of total storms

aux = which('gurobi_setup');
if ~isfield(opts, 'gurobi')
    if isempty(aux); opts.gurobi = 'no';
    else;            opts.gurobi = 'yes'; end
end

if ~exist('stmRate_pre', 'var'); stmRate_pre = stmRate; end
if isempty(stmRate_pre);         stmRate_pre = stmRate; end

if opts.nselect > 0
    if (sum(stmRate == 0) > 0)
        if (opts.runtime_rateadjustment*sum(stmRate == 0)*(Ns-sum(ind_storm)) < opts.runtime_stormselection)
            %% Adjustment
            fprintf('Performing final Adjustment to reduce the number of zero-rate storms... \n')
            wb   = PoolWaitbar(sum(stmRate == 0)*(Ns-sum(ind_storm)), 'Adjusting final results');
            while(sum(stmRate == 0) > 0)
                ind_select      = find(ind_storm);  % indices of storms selected
                ind_select_zero = ind_select(stmRate == 0);   % indices of storms selected with zero rates
                ind_select_not  = setdiff(1:Ns, ind_select);  % indices of storms not selected
                
                ind_storm(ind_select_zero(1)) = 0;
                for ii = 1:length(ind_select_not)
                    ind_storm(ind_select_not(ii)) = 1;
                    % tic
                    if strcmp(opts.gurobi, 'yes') % Using gurobi optimizer
                        [d, stmRate] = HC_match_LargeA_gurobi(ind_storm,lp_model); % if you has gurobi optimization
                    else
                        [d, stmRate] = HC_match_LargeA(ind_storm,lp_model);
                    end
                    % toc

                    if (stmRate(sum(ind_storm(1:ind_select_not(ii)))) > 0)&&(sum(stmRate == 0) < length(ind_select_zero))
                        break
                    else
                        ind_storm(ind_select_not(ii)) = 0;
                    end
                    increment(wb)
                end
                if (ii == length(ind_select_not))&&(ind_storm(ind_select_not(ii)) == 0)
                    fprintf('WARNING: Final Adjustment is interrupted \n')
                    fprintf('WARNING: Number of zero-rate storms cannot be reduced \n')
                    fprintf('WARNING: Final Adjustment is not recommended \n')
                    fprintf('WARNING: Add more target or intermediate AEPs and solve again the original storm selection problem \n')
                    fprintf('WARNING: Check also if the number of clusters is too small. \n')
                    break;
                end
            end
            delete(wb)
            fprintf('Now there are %d storms with zero rates. \n', sum(stmRate == 0));
        elseif sum(stmRate_pre == 0) > sum(stmRate == 0)
            fprintf('WARNING: Final Adjustment is not recommended \n')
            fprintf('WARNING: Final Adjustment would take longer time than original storm selection problem \n')
            fprintf('WARNING: Check first if the number of clusters is too small \n')
            fprintf('WARNING: Then consider adding more target or intermediate AEPs and solve again the original storm selection problem. \n')
        else
            fprintf('WARNING: Final Adjustment is not recommended \n')
            fprintf('WARNING: Final Adjustment would take longer time than original storm selection problem \n')
            fprintf('WARNING: Add more target or intermediate AEPs and solve again the original storm selection problem \n')
            fprintf('WARNING: Check also if the number of clusters is too small. \n')
        end
    end
else
    fprintf('WARNING: Change storms included by user to reduce number of storms with zero rates \n')
    fprintf('WARNING: Solving storm selection problem is recommended. \n')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d, stmRate] = HC_match_LargeA(xs, model)
% function to calculate sum of absolute deviations to the original hazard
% curve (objective function for the optimization) and adjust rates of
% selected storms as well
% -- Inputs:
%   xs:  a binary matrix of size [1 x number of entire storms] with 1
%        (selected) of 0 (not selected)
%   model:  a structure with model for linear programming with following fields.
%            .P:   sparse indicator matrix where 1 means that surge response is 
%                  greater than surge threshold
%            .A:   sparse matrix for the constraint of the linear programming
%            .rhs: right hand side of the constraint
%            .obj: vector of coefficients for the objective function
%            .Ntsr:     number of target surge rates for each node of interest
%            .tgtRates: matrix with each column being target rates for each
%                       node of interest
%
%            .weights:  weights considered in objective function
%
% -- Outputs:
%   d:  sum of absolute deviations to the original hazard curve 
%       (objective function for the optimization)
%   stmRate:  adjusted rates of selected storms
%
%  V1    09/05/2023 WoongHee Jung (wjung2@nd.edu)
%  V1.1  02/20/2024 WoongHee Jung (wjung2@nd.edu)
%        Modified to accommodate weights for nodes and rates inside the
%        function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indices to keep in Alarge
indi = logical([xs, ones(1, 2*sum(model.Ntsr))]);
indj = logical([ones(1, sum(model.Ntsr)), xs]);

% matrix A
model.A = model.A(indi, indj);

% options for Linear programming
opt.linopts = optimoptions(@linprog, 'Display', 'off');

% linear programming (objective func: f'*x such that Ax <= b)
xopt = linprog(model.obj, model.A, model.rhs, [], [], [], [], opt.linopts);

d = sum((model.weights_nodes.*model.weights_rates)*xopt(1:sum(model.Ntsr))); % Objective function
stmRate = xopt(sum(model.Ntsr)+1:end);       % Optimal storm rates

end

function [d, stmRate] = HC_match_LargeA_gurobi(xs, model)
% function to calculate sum of absolute deviations to the original hazard
% curve (objective function for the optimization) and adjust rates of
% selected storms as well
% -- Inputs:
%   xs:  a binary matrix of size [1 x number of entire storms] with 1
%        (selected) of 0 (not selected)
%   model:  a structure with model for linear programming with following fields.
%            .P:   sparse indicator matrix where 1 means that surge response is 
%                  greater than surge threshold
%            .A:   sparse matrix for the constraint of the linear programming
%            .rhs: right hand side of the constraint
%            .obj: vector of coefficients for the objective function
%            .Ntsr:     number of target surge rates for each node of interest
%            .tgtRates: matrix with each column being target rates for each
%                       node of interest
%
%            .weights_nodes:  weights for nodes considered in objective function
%            .weights_rates:  weights for rates considered in objective function
%
% -- Outputs:
%   d:  sum of absolute deviations to the original hazard curve 
%       (objective function for the optimization)
%   stmRate:  adjusted rates of selected storms
%
%  V1    09/05/2023 WoongHee Jung (wjung2@nd.edu)
%  V1.1  02/20/2024 WoongHee Jung (wjung2@nd.edu)
%        Modified to accommodate weights for nodes and rates inside the
%        function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indices to keep in Alarge
indi = logical([xs, ones(1, 2*sum(model.Ntsr))]);
indj = logical([ones(1, sum(model.Ntsr)), xs]);

% matrix A
model.A = model.A(indi, indj);

% linear programming (objective func: f'*x such that Ax <= b)
params.OutputFlag = 0;
result = gurobi(model, params);

d = sum((model.weights_nodes.*model.weights_rates)*result.x(1:sum(model.Ntsr))); % Objective function
stmRate = result.x(sum(model.Ntsr)+1:end);       % Optimal storm rates

end