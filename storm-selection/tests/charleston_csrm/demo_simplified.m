% Selection of subset of storms to match a desired hazard curve. The
% hazard curve is defined through the larger set of storms and then a
% smaller subset is selected to match that one. Specific surge rates of interest
% are chosen along the hazard curve.

% The optimization process is as follows:
% Genetic algorithm finds optimal subset of storms with their optimal rates
% obtained by linear programming. The objective is to minimize the sum of
% the absolute deviation of the estimated surge rates for all the nodes or representative nodes 
% to the target surge rates.

% The linear programming is
%   minimize    f'*x
%   such that   Ax <= b
%
% where x = [absolute deviations(1x(NtxNn) vector), storm rates (1xNs vector)]';
%       f = [kron(weight_nodes, weight_rates), zeros(1, Ns)]';  % 'kron(A,B)' is Kronecker product of matrices A and B.
%       A = [zeros(Ns, Nt*Nn),  eye(Ns);     % inequality constraint that storm rates are greater than or equal to 0
%           -eye(Nt*Nn),        P;           % inequality constraint necessary to formulate absolute deviations 
%           -eye(Nt*Nn),       -P];          % inequality constraint necessary to formulate absolute deviations 
%       P = indicator matrix if response is greater than a threshold, 1,
%           otherwise 0. ((NtxNn)xNs) matrix)
%       b = [zeros(1, Ns), repmat(tgtrate, 1, Nn), -repmat(tgtrate, 1, Nn)]'
%    weight_nodes: Nn-dimensional vector of weights for nodes
%    weight_rates: Nt-dimensional vector of weights for rates
%    tgtrate: vector of target rates
%    Nt: number of target rates
%    Nn: number of nodes considered
%    Ns: number of storms
%    Sizes of vectors and matrices can change if anchor points are
%    considred.
%
%  V1   02/22/2024 WoongHee Jung (wjung2@nd.edu)
%  V1.1 10/15/2024 WoongHee Jung (wjung2@nd.edu)
%                  added validation section that runs
%                  "validation_statistics_HPSO.m" to obtain validation
%                  statistics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;

%% User selections
% Selection of nodes of interest within domain
opts.node_interest = [];  % If you want to focuse on a subset of the domain, put indices of nodes in the subset.
                          % If it is empty ([]), all the nodes are used.
opts.node_wetness  = [];  % criterion for the probability of being wet. 
                          % Nodes with the probabilities below than 'opts.node_wetness' will be ignored.
                          % if empty then not used (equivalent to being zero)

% Selection of representative points to reduce computational complexity of
% the storm selection. This only impacts the way storms are chosen and not how
% the probability masses of these storms are assigned (the latter uses all
% the nodes). If you have more than 1000 nodes it is recommended to perform clustering
% if no additional storms need to be identified these options are not necessary
opts.cluster.perform = 'no';  % perform clustering 'yes' or not 'no' (default is 'yes')
                               % if 'yes' then the following quantities are also needed
opts.cluster.k = 0;  % number of clusters (representative nodes) you want to use. A value of 300-500 is a good choice here. 
                       % The larger the value the higher the compexity of choosing storms, but the higher accuracy as well  
opts.cluster.type = 'comb';  % perform clustering using geospatial information only ('geo'), surge response only ('resp') or both ('comb').
                             % if it is 'geo', the grid of the points is also
                             % required. If not porived the type will be overridden
                             % eventually to be 'resp'

% Storm selection parameters
opts.tgtRate  = 1./[1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000];  % target surge rates on hazard curves you wants to match
% To better match hazard curves, intermediate AEPs (with smaller weights than those of actual target AEPs) can be added.  
opts.int_tgtRate.add  = 'yes';   % Add intermediate AEPs ('yes') between the target AEPs
opts.int_tgtRate.sf   = 1.5;     % Safety factor to find number of intermediate
                                 % AEPs needed. default is 1.5.
                                 % If high correlation of responses is
                                 % concerned, increase the value to 3.0
                                 % (5.0, if cluster analysis is performed)                              
opts.int_tgtRate.wgt  = 0.1;     % Relative weight of intermediate AEPs. 
                                 % Set this to a small number to primarily
                                 % focus on the actual target AEPs defined
                                 % above ('opts.tgtRate') in the storm
                                 % selection. default is 0.1.
opts.nselect  = 80;  % number of storms you want to choose
                     % if user wants to adjust only the rates of the storms that
                     % have already been selected (defined by the following
                     % variable 'opts.storms_included'), set to 0.
opts.storms_included = [];   % indices of storms that needs to be included in
                             % the final set. Additional 'opts.nselect' storms
                             % are chosen beyond this set.
                             % if opts.nselect=0 and opts.storms_included=[], 
                             % by default, this will be adjusted to be 1:Ns 
                             % (rates for all Ns storms will be adjusted).

% options for weights for storm rates and nodes
opts.weight_nodes = [];     % weight for nodes. set a row vector of weights for 
                            % nodes of size [1 x Nn]. If it is empty ([]),
                            % uniform weights, ones(1, Nn) is used.
opts.weight_rates = [];     % weight for target rates. set a row vector of weights for 
                            % target rates of size [1 x length(opts.tgtRate)]. 
                            % If it is empty ([]), ones(1, length(opts.tgtRate) is used.
opts.error        = 'rel';  % calculate objective function considering 
                            % relative magnitude of target rates ('rel') 
                            % or not ('abs'). Default is 'rel'.

% options for optimization (Storm Selection)
opts.objective           = 'Sdev'; % Objective for storm selection
                                   % 'Sdev': sum of absolute devations (default)
                                   % 'Scor': spatial correlation

% option for optimizer
opts.gurobi = 'yes';     % Use gurobi optimizer ('yes') or not ('no').
                          % If you want use gurobi optimizer, the optimizer 
                          % needs to be installed and a license is required.


%% NOTE: addpath or other MATLAB path modification not compatible with stand-alone 
%% addpath(genpath('/mnt/hgfs/rdchlmyl/projects/chart/storm-selection/tests/charleston_csrm'))addpath(genpath('/mnt/hgfs/rdchlmyl/projects/chart/storm-selection/tests/charleston_csrm'))

folder_path = "/mnt/hgfs/rdchlmyl/projects/chart/storm-selection/tests/charleston_csrm/";
input_path = folder_path + 'input/final/';

%% Uncomment line for native testing in MATLAB
addpath(genpath('../../src'))

%% load response data
load(input_path + "Charleston_CSRM_Subset_Output_WithTides.mat")       % Response 'Resp' and Indicator matrix 'Ind_wet' telling whether node is wet or dry

% At this point you should have
% 'Resp' is a matrix (number of storms) x (number nodes) with the response
%
% you might also load 'Ind_wet', a matrix indicating whether a node is wet (1) or dry (0) for each storm.
% if not loaded, the 'Ind_wet' matrix will be created.

% Load the hazard curve to match if desired 
% 'HC_target' is a cell of size [number of nodes x 1] with each element being
% a structure with following fields.
%           .x:  surge level thresholds corresponding to probabilities of
%                exceedance (.Pf)
%           .Pf: probabilities of exceedance
% Note that the lengths of HC_target.x and .Pf have to match.
% You can either upload 'HC_target' from a file or define it based on 'Resp' and
% 'ProbMass'. If 'HC_target' is not provided from a file, then it will be calculated. 
% In this case you also need the storm rates 'ProbMass'
% load('Charleston_CSRM_Subset_ProbMass.mat') % Storm rates 'ProbMass'

% 'ProbMass' a vector with the original probability masses of the storms


% if you will choose a subset of nodes based on k-means SPATIAL clustering 
% then also load the latitute/longitude grid
load(input_path + "Charleston_CSRM_Subset_Grid.mat")
% 'gridSM' is a matrix (number of nodes) x 2 with the latitude and longitude
% of the grid 

% load the target hazard curves for the restricted savepoints
load(input_path + "SACS_NCSFL_HC_tbl_ATC_SWL_SLC0.mat", 'HC_tbl_ATC_SWL_SLC0')
load(input_path + "AEF22.mat", 'tbl_aef')
[~, idx, ~] = intersect(tbl_aef, opts.tgtRate);
idx = sort(idx);
HC = HC_tbl_ATC_SWL_SLC0.BE_22;
HC = HC(gridSM(:, 3), idx);
HC_target = cell(height(HC), 1);
for i = 1 : height(HC_target)
    % it seems to be from smallest probability to largest, so fliplr
    temp.x = fliplr(HC(i, :));
    temp.Pf = fliplr(opts.tgtRate);
    HC_target{i} = temp;
end
clear HC_tbl_ATC_SWL_SLC0
%


if ~exist('Ind_wet', 'var');   Ind_wet = []; end % If 'Ind_wet' does not exist, it will be created.
if ~exist('gridSM', 'var');    gridSM  = []; end % If 'gridSM' does not exist, make it empty.
if ~exist('HC_target', 'var') 
    HC_target = []; % If 'HC_target' does not exist, it will be calculated with 'Resp and 'ProbMass'.
    if ~exist('ProbMass','var'); error('The original probability masses of the storms (ProbMass) are necessary'); end
else 
    ProbMass = []; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% After this point you simply run the remaining code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set problem for storm selection
[Resp, opts, HC_target, Ind_wet, gridSM] = set_problem_StormSelection(Resp, opts, HC_target, ProbMass, Ind_wet, gridSM);
% 'Resp' is surge response matrix with some meaningless nodes being removed. 
% 'opts' is a structure that contains all the information and options necessary
%        for the optimization
% 'HC_target' is a cell with target hazard curves. If it was not loaded, it
%             was calculated.
% 'Ind_wet'   is a matrix indicating whether a node is wet (1) or dry (0) for each storm.
%             If it was not loaded, it was calculated.
% 'gridSM'    is grid information with some meaningless nodes being removed 
%             if it was provided.

%% Storm selection
[ind_storm, lp_model, opts] = StormSelection(Resp, HC_target, opts);
% 'ind_storm' is a vector of size [1 x Ns] with 1 (the storm is selected) or 
%                0 (not selected).
% 'lp_model' is a structure with fields that are used in optimization
%               for selecting storms and adjusting their rates (this might be utilized in next section).

select_storms=find(ind_storm); % indices of chosen storms

%% Adjust rates of selected storms
if strcmp(opts.cluster.perform, 'yes')
    opts.cluster.perform = 'no';
    opts.SS_nodes_use    = 1:size(Resp, 2);
end
[stmRate, lp_model, opts, d, SD] = StormRatesAdjustment(ind_storm, Resp, HC_target, opts, lp_model);
% 'stmRate' is a vector with rates for select_storms indices
% 'runtime_rateadjustment' is runtime to complete storm rates adjustment. 
%                          This will be used to make a suggestion in the
%                          Final Adjustment section.

%% Final adjustment (remove zero rate storms)
% If there are selected storms with zero rates, they can be replaced with
% storms that were not selected to achieve nonzero rates for all the storms
% selected.
% WARNING: storms finally selected after this final adjustment are
% suboptimal. If you want to obtain 'better' storms, add even more AEPs (either
% target or intermediate) and do the storm selection again.
[ind_storm, stmRate, lp_model] = FinalAdjustment(ind_storm, stmRate, lp_model, opts);
select_storms=find(ind_storm); % indices of chosen storms

%% Validation
[ES_tot, ES_AEP] = validation_statistics_HPSO(ind_storm, stmRate, Resp, HC_target, opts, lp_model);
% 'ES_tot' is a structure with global error statistics calculated across all locations 
%          and all target AEPs with following fields
%          .devP: weighted average of absolute deviation of estimated AEP levels 
%                 from target AEP levels 
%          .corP: (weighted) correlation coefficient of estimated surge thresholds
%                 to true thresholds
%          .MAE: average of absolute deviation of estimated
%                surge thresholds from true thresholds
%          .MAPE: average of percentage of absolute deviation of
%                 estimated surge thresholds from true thresholds
% 'ES_AEP' is a structure with error statistics per each target AEP level
%          calculated across all locations with following fields
%          .devP: weighted average of absolute deviation of estimated AEP levels 
%                 from target AEP levels 
%          .corP: (weighted) correlation coefficient of estimated surge thresholds
%                 to true thresholds
%          .MAE: median of absolute deviation of estimated surge thresholds 
%                from true thresholds
%          .MAPE: median of percentage of absolute deviation of estimated 
%                 surge thresholds from true thresholds

%% Save Storm Selction results
save('Charleston_CSRM_Storms_and_Stats.mat','select_storms','stmRate','Ind_wet','keep_storms',...
                    'Resp','ES_tot','ES_AEP','HC_target','opts','-v7.3')

%% END