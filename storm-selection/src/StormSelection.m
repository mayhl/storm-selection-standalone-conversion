function [ind_storm, lp_model, opts, d, stmRate, SD] = StormSelection(Resp, HC_target, opts)
% Selection of subset of storms to match a desired hazard curve. The
% hazard curve is defined through the larger set of storms and then a
% smaller subset is selected to match that one. Specific surge rates of interest
% are chosen along the hazard curve.
% -- Inputs:
%  Resp:  surge responses for storms in the original suite (Ns x Nn matrix;
%         Ns is number of storms, Nn is number of nodes)
%
%  HC_target:  a cell of size Nn with each element corresponding to 
%              a structure having information of target hazard curve. If
%              empty, target hazards curves are estimated using surge
%              responses 'Resp' and probability masses for storms in the original suite
%              'ProbMass'. 
%              Each structure has following fields,
%            .x:  surge levels
%            .Pf: probabilities of exceedance
%
%  opts:  a structure with options for the storm selection as follows
%       .tgtRate:     target surge rates.  
%                     default is 1./[10, 20, 50, 100, 200, 500, 1000].
%       .interp:      scale used to interpolate hazard curves to find surge
%                     levels corresponding to target surge rates.
%                     linear interpolation in the original scale (1) or
%                     logscale (2). default is 2.
%       .HC_target_given: if the target hazard curves are given, 1,
%                            otherwise, 2 (Calculated).
%
%      for optimization in general
%       .nselect:         number of storms to select. default is min(50, Ns/10).
%       .storms_included: indices of storms that needs to be included in
%                             the final set. Additional 'opt.nseclose lect' storms
%                             are chosen beyond this set.
%
%      for genetic algorithm
%       .objective:    objective function - ('Sdev') Sum of absolute
%                      deviations and ('Scor') spatial correlation between
%                      hazard maps across different target AEPs. default is
%                      'Sdev'.
%       .npop:         population size. default is 300.
%       .generations:  maximum generation. default is 1000.
%      <STOPPING CRITERIA>
%       .FunctionTolerance: The algorithm stops if the average relative change in 
%                           the best fitness function value over 'MaxStallGenerations' 
%                           is less than or equal to FunctionTolerance. 
%                           default is 1e-9.
%       .MaxStallGenerations: The algorithm stops if the average relative change in 
%                             the best fitness function value over 'MaxStallGenerations' 
%                             is less than or equal to FunctionTolerance. 
%                             default is 50.
%       .MaxTime: the maximum time in seconds the genetic algorithm runs
%                 before stopping. default is Inf.
%
%      for weights for objective function of linear programming
%       .weight_rates  set a row vector of weight for target rates of size that 
%                      has to be same as the size of opts.tgtRate. For
%                      rates of anchor points, weights are linearly
%                      interpolated based on 1./opts.tgtRate. 
%                      default is 1./opts.tgtRate*opts.tgtRate(1).
%       .weight_nodes  set a row vector of weight for nodes of size [1 x Nn].
%                      default is ones(1, Nn).
%       .error         calculate objective function considering target rates relatively
%                      (inversly proportional to target rates) with 'rel'
%                      or not with 'abs'. Default is 'rel'.
%
%      for optimizer (Linear programming)
%       .gurobi        Use gurobi optimizer ('yes') or not ('no').  
%                      If you want use gurobi optimizer, the optimizer 
%                      needs to be installed and a license is required. 
%                      If it is empty, default will be 'yes' if the gurobi optimizer 
%                      is properly setup, or 'no'.
%
% -- Outputs:
%  ind_storm:  binary vector with 1 indicating 'included storm'.
%
%  lp_model:   a structure with model for linear programming with following fields.
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
%            .weights_Corr_nodes:  weights for nodes considered in 'Scor' objective function
%            .weights_Corr_rates:  weights for rates considered in 'Scor' objective function
%
%  opts:  opts with following updated fields
%            .A, .b, .lb, .ub, .INTCON: genetic algorithm characteristics
%            .gaopts:  options for genetic algorithm
%            .runtime_stormselection: runtime to complete storm selection
%
%  d:         objective for outer-loop optimization
%
%  stmRate:   a vector with adjusted storm rates. If clustering has not been
%             performed, this is actual adjusted storm rates, but if
%             clustering has been performed, this is a preliminary result so
%             adjustment of selected storms should be performed again for the
%             original problem (considering all the nodes of interest).
%
%  SD:        The sum of absolute deviations across different target AEPs
%             and different locations (If opts.objective is 'Sdev', 
%             SD is identical to d). 
%
%  V1   10/23/2023 WoongHee Jung (wjung2@nd.edu)
%  V2   02/21/2024 WoongHee Jung (wjung2@nd.edu)
%       1) Modified 'resp' option for cluster analysis by including hazard
%       information (surge thresholds corresponding to target AEPs)
%       2) Modified to consider the number of nodes belonging to each
%       cluster as weight for the representative node in the optimization
%       3) Modified to accommodate a different objective for the outer-loop 
%       optimization such as using spatial correlation across different AEPs 
%       with additional option ("opts.objective"). 
%  V2.1 03/22/2024 WoongHee Jung (wjung2@nd.edu)
%       Modified identification of anchor points. It is conducted only when
%       the target hazard curves are calculated based on given 'Resp' and
%       'ProbMass'.
%  V2.2 03/28/2024 WoongHee Jung (wjung2@nd.edu)
%       Modified two functions "HC_match_LargeA_gurobi" and "HC_match_LargeA"
%       to accommodate the modification in the function
%       "TargetRatesNThreshs_LargeA_mod".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.nselect > 0 % SELECT STORMS
    % Use nodes of interest only
    Resp = Resp(:, opts.SS_nodes_use);
    HC_target = HC_target(opts.SS_nodes_use);
    [Ns, Nn] = size(Resp);
    
    %% Find surge thresholds corresponding to target surge rates
    Nt = length(opts.tgtRate);
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
    Rates(isnan(Rates)) = [];
    lp_model.tgtThreshs = tgtThreshs;
    lp_model.anchorpts  = (abs(tgtRates - repmat(tgtRate', 1, Nn)) > eps);

    %% Construct LARGE Indicator matrix P (of size [sum(Ntsr) x Ns])
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
    %%%%%%%%%%%%%%%%%%% i index, j index,                elements,    i size,  j size
    lp_model.P = sparse(targeti,  stormi, ones(length(stormi), 1), sum(Ntsr),      Ns);

    %% Construct FULL A matrix for linear programming (constraint)
    % Construct a sparse matrix for A by assembling matrices
    lp_model.A = [sparse(Ns, sum(Ntsr)), -speye(Ns);
                      -speye(sum(Ntsr)),  lp_model.P;
                      -speye(sum(Ntsr)), -lp_model.P];

    %% Save other variables in lp_model for linear programming
    lp_model.rhs  = [zeros(1, opts.nselect+length(opts.storms_included)), Rates, -Rates]';  % right hand side of the constraint
    lp_model.Ntsr = Ntsr;  % Number of target surge rates for each node of interest
    lp_model.tgtRates = tgtRates; % matrix with each column being target rates for each node of interest

    % Update weights optimization
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
    
    if strcmp(opts.objective, 'Scor')
        if strcmp(opts.cluster.perform, 'yes')
            lp_model.weights_Corr_nodes = opts.cluster.weight_nodes;
        else
            lp_model.weights_Corr_nodes = opts.weight_nodes(opts.SS_nodes_use);
        end
        lp_model.weights_Corr_rates = opts.weight_rates;
    end

    lp_model.obj  = [(weight_nodes.*weight_rates)/sum(weight_nodes.*weight_rates_unorm), zeros(1, opts.nselect+length(opts.storms_included))]';  % coefficients for the objective function

    %% Set optimization problem
    opts.totn = Ns; % Total number of storms
    % Get optimization ready
    opts = init_opt_GA(opts);

    %% Perform optimization
    if strcmp(opts.objective, 'Sdev')
        ind_obj = 1;
    elseif strcmp(opts.objective, 'Scor')
        ind_obj = 2;
    end
    
    disp('Running Genetic algorithm to select optimal storms ...')
    if strcmp(opts.gurobi, 'yes') % Using gurobi optimizer
        tic
        [ind_storm, ~] = ga(@(x)HC_match_LargeA_gurobi(x,lp_model,Resp,ind_obj,interp,HC_target_given),opts.totn,opts.A,opts.b,[],[],opts.lb,opts.ub,[],opts.INTCON,opts.gaopts);
        t = toc
        opts.runtime_stormselection = t;

        [d, SD, stmRate] = HC_match_LargeA_gurobi(ind_storm,lp_model,Resp,ind_obj,interp,HC_target_given);
    else
        tic
        [ind_storm, ~] = ga(@(x)HC_match_LargeA(x,lp_model,Resp,ind_obj,interp,HC_target_given),opts.totn,opts.A,opts.b,[],[],opts.lb,opts.ub,[],opts.INTCON,opts.gaopts);
        t = toc
        opts.runtime_stormselection = t;

        [d, SD, stmRate] = HC_match_LargeA(ind_storm,lp_model,Resp,ind_obj,interp,HC_target_given);
    end

    %% Postprocessing
    if strcmp(opts.cluster.perform, 'yes')
        lp_model = [];
    end
    disp('Storm selection is completed.')

else % No storm selection
    % If there is no specific storms that user would like to include,
    % include all storms
    Ns = size(Resp, 1);
    if isempty(opts.storms_included); opts.storms_included = 1:Ns; end
    ind_storm = zeros(1, Ns);
    ind_storm(opts.storms_included) = 1;

    lp_model=[]; % Set 'lp_model' empty
    disp('No additional storms are selected.')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function opt = init_opt_GA(opt)
% update opt structure with
% opt.totn: number of discrete design variables 
% opt.{A,B,lb,ub,INTCON}: genetic algorithm characteristics
% opt.gaopts: optimization structure for genetic algorithms
% opt.{linf,linA,linb}: linear programming characteristics
% opt.linopts: optimization structure for linear programming
%
%  V1   09/01/2023 WoongHee Jung (wjung2@nd.edu)
%  V2   09/28/2023 WoongHee Jung (wjung2@nd.edu)
%                  implemented user defined functions for crossover ('pleaseselectstorms') and
%                  mutation ('mutationstormselection') in Genetic algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control the randomness of initial population                  
rng('default')

Ns_included = length(opt.storms_included);
storms_available = setdiff(1:opt.totn, opt.storms_included);

% set up the number of design variables, corresponding lower/upper bound, and dictate discrete variables
opt.A = [ones(1,opt.totn);-ones(1,opt.totn)];  % A and b guarantees selection up to nselect storms
opt.b = [opt.nselect+Ns_included; -opt.nselect-Ns_included];
opt.lb = zeros(opt.totn,1);      % lb, ub and INTCON guarantess binary variable of 0 and 1
opt.lb(opt.storms_included) = 1;
opt.ub = ones(opt.totn,1);
opt.INTCON = 1:opt.totn;       %discrete variables 

%create initial population for genetic algorithm optimization
initPop = zeros(opt.npop,opt.totn);
for i = 1:opt.npop; initPop(i,[randsample(storms_available,opt.nselect), opt.storms_included]) = 1; end

% set characteristics of genetic algorithm optimization
opt.gaopts = optimoptions(@ga, 'FunctionTolerance', opt.FunctionTolerance, 'MaxGenerations', opt.generations, 'MaxStallGenerations', opt.MaxStallGenerations,'MaxTime', opt.MaxTime, ...
                               'InitialPopulationMatrix', initPop, 'PopulationSize', opt.npop, 'PlotFcns', @gaplotbestfun, 'UseParallel', true, ...
                               'CrossoverFcn', @pleaseselectstorms, 'MutationFcn', @mutationstormselection);
end

function xoverkids = pleaseselectstorms(parents, ~, GenomeLength, ~, ~, thisPopulation)
% Crossover function for selecting storms problem
% (1) Storms that both parents have are passed to their kids.
% (2) The rest of storms that parents have are randomly picked and
%     passed to their kids.

% How many kids to produce?
nKids = length(parents)/2;

% Allocate space for the kids
xoverkids = zeros(nKids,GenomeLength);

% Number of ones (number of storms want to select)
Nones = sum(thisPopulation(parents(1), :), 2);

% For each kid
for i = 1:nKids   
    common = all([thisPopulation(parents(2*i-1), :); thisPopulation(parents(2*i), :)]);
    xoverkids(i, :) = double(common);
    if sum(common) < Nones
        pot_ind = setdiff(unique([find(thisPopulation(parents(2*i-1), :)), find(thisPopulation(parents(2*i), :))]), ...
            find(common));
        sel_ind = datasample(pot_ind, Nones-sum(common), 'Replace', false);
        xoverkids(i, sel_ind) = 1;
    end
end

end

function mutationChildren = mutationstormselection(parents,options,GenomeLength,~,~,~,thisPopulation)
% Mutation function for storm selection
% Keep storms that need to be included ('storms_included').
% Randomly replace one selected storm except those kept with one unselected storm.

% Number of mutation children
numMutation = length(parents);

% Variable bounds and range. The mutation function requires finite bounds.
lb = options.LinearConstr.lb;
ub = options.LinearConstr.ub;

% Number of ones (number of storms want to select)
Nones = sum(thisPopulation(parents(1), :), 2);

free_ind = find(abs(ub - lb) > eps); % find free indicies that are indices for 
                                     % storms that are not pre-determined to be 
                                     % included.
Nfree    = length(free_ind);
Nfixed   = length(thisPopulation(parents(1), :))-Nfree;
% Initialize mutation chldren
mutationChildren = zeros(numMutation, GenomeLength);

for i = 1:numMutation
    ind_selected = find(thisPopulation(parents(i), :));
    mutationChildren(i, setdiff(ind_selected, free_ind)) = 1; % set 1 to indices for storms 
                                                              % that are pre-determined to be 
                                                              % included.
    
    % Randomly replace one selected storm with one unselected storm
    ind_selected_selected = datasample(ind_selected, Nones-Nfixed-1, 'Replace', false);
    mutationChildren(i, ind_selected_selected) = 1;
    ind_unselected_selected = datasample(setdiff(free_ind, ind_selected), 1, 'Replace', false);
    mutationChildren(i, ind_unselected_selected) = 1;
end
end

function [d, SD, stmRate] = HC_match_LargeA(xs, model, Resp, ind_obj, interp, HC_given)
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
%            .weights_nodes:  weights for nodes considered in 'Sdev' objective function
%            .weights_rates:  weights for rates considered in 'Sdev' objective function
%            .weights_rates_unorm:  unormalized weights by rates considered in 'Sdev' objective function
%            .weights_Corr_nodes:  weights for nodes considered in 'Scor' objective function
%            .weights_Corr_rates:  weights for rates considered in 'Scor' objective function
%   resp:   surge response. This is necessary only when the objective in
%           the outer-loop optimization corresponds to spatial correlation between
%           hazard maps.
%
%   ind_obj: If 1, the objective in the outer-loop optimization is 'Sdev',
%            elseif 2, 'Scor'.
%
%   interp:  To find surge thresholds corresponding to 'tgtRate' the hazard
%            curve is interpolated. It is conducted with linear
%            interpolation in the original scale (1) or in the logscale (2).
%
%   HC_given: an Indicator that is 1 if target hazard curves (HC_target) is
%             given, else 2 if target hazard curves is calculated based on
%             Resp and ProbMass.
%             Anchor points are identified when HC_given = 2.
%
% -- Outputs:
%   d:        objective for outer-loop optimization
%   SD:       The sum of absolute deviations across different target AEPs
%             and different locations (If opts.objective is 'Sdev', 
%             SD is identical to d). 
%   stmRate:  adjusted rates of selected storms
%
%  V1   09/05/2023 WoongHee Jung (wjung2@nd.edu)
%  V2   02/21/2024 WoongHee Jung (wjung2@nd.edu)
%       Modified to accommodate a different objective for the outer-loop 
%       optimization such as using spatial correlation across different AEPs 
%       with additional option ("opts.objective").
%  V2.1 03/28/2024 WoongHee Jung (wjung2@nd.edu)
%       Modified to handle additional input arguments for function
%       "TargetRatesNThreshs_LargeA_mod".
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

stmRate = xopt(sum(model.Ntsr)+1:end);       % Optimal storm rates

if ind_obj == 1
    d = sum((model.weights_nodes.*model.weights_rates/sum(model.weights_nodes.*model.weights_rates_unorm))*xopt(1:sum(model.Ntsr))); % Objective function
    SD = d;
elseif ind_obj == 2
    SD = sum((model.weights_nodes.*model.weights_rates/sum(model.weights_nodes.*model.weights_rates_unorm))*xopt(1:sum(model.Ntsr))); % Objective function
    if sum(stmRate < 0) > 0
        stmRate(stmRate<0) = eps;
    end

    % Reconstructed Hazard Curve
    Nn = size(Resp, 2);
    Resp = Resp(logical(xs), :);
    Resp(isnan(Resp)) = -Inf;
    tgtThreshs_rec = zeros(size(model.tgtThreshs));
    tgtRates = model.tgtRates;
    parfor ii = 1:Nn
        Resp_temp = Resp(:, ii);
        % sort the response and define hazard statistics
        [~, sorti] = sort(Resp_temp, 'descend');
        HC_rec_x  = Resp_temp(sorti)';
        HC_rec_Pf = cumsum(stmRate(sorti))';
        % address dry isntances
        removei = HC_rec_x == -inf;
        HC_rec_x(removei)  = [];
        HC_rec_Pf(removei) = [];

        Resp_temp(Resp_temp == -inf) = NaN;

        % Thresholds
        [~, tgtThreshs_rec(:, ii), ~] = TargetRatesNThreshs_LargeA_mod(Resp_temp, HC_rec_Pf, HC_rec_x, tgtRates(:, ii), interp, HC_given);
    end
    d = -model.weights_Corr_rates*spatial_corr(model.tgtThreshs, tgtThreshs_rec, model.anchorpts, model.weights_Corr_nodes);
end
end

function [d, SD, stmRate] = HC_match_LargeA_gurobi(xs, model, Resp, ind_obj, interp, HC_given)
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
%            .weights_nodes:  weights for nodes considered in 'Sdev' objective function
%            .weights_rates:  weights for rates considered in 'Sdev' objective function
%            .weights_rates_unorm:  unormalized weights by rates considered in 'Sdev' objective function
%            .weights_Corr_nodes:  weights for nodes considered in 'Scor' objective function
%            .weights_Corr_rates:  weights for rates considered in 'Scor' objective function
%   resp:   surge response. This is necessary only when the objective in
%           the outer-loop optimization corresponds to spatial correlation between
%           hazard maps.
%
%   ind_obj: If 1, the objective in the outer-loop optimization is 'Sdev',
%            elseif 2, 'Scor'.
%
%   interp:  To find surge thresholds corresponding to 'tgtRate' the hazard
%            curve is interpolated. It is conducted with linear
%            interpolation in the original scale (1) or in the logscale (2).
%
%   HC_given: an Indicator that is 1 if target hazard curves (HC_target) is
%             given, else 2 if target hazard curves is calculated based on
%             Resp and ProbMass.
%             Anchor points are identified when HC_given = 2.
%
% -- Outputs:
%   d:        objective for outer-loop optimization
%   SD:       The sum of absolute deviations across different target AEPs
%             and different locations (If opts.objective is 'Sdev', 
%             SD is identical to d). 
%   stmRate:  adjusted rates of selected storms
%
%  V1   09/05/2023 WoongHee Jung (wjung2@nd.edu)
%  V2   02/21/2024 WoongHee Jung (wjung2@nd.edu)
%       Modified to accommodate a different objective for the outer-loop 
%       optimization such as using spatial correlation across different AEPs 
%       with additional option ("opts.objective").  
%  V2.1 03/28/2024 WoongHee Jung (wjung2@nd.edu)
%       Modified to handle additional input arguments for function
%       "TargetRatesNThreshs_LargeA_mod".
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

stmRate = result.x(sum(model.Ntsr)+1:end);       % Optimal storm rates

if ind_obj == 1
    d = sum((model.weights_nodes.*model.weights_rates/sum(model.weights_nodes.*model.weights_rates_unorm))*result.x(1:sum(model.Ntsr))); % Sum of absolute deviations
    SD = d;
elseif ind_obj == 2
    SD = sum((model.weights_nodes.*model.weights_rates/sum(model.weights_nodes.*model.weights_rates_unorm))*result.x(1:sum(model.Ntsr))); % Sum of absolute deviations
    if sum(stmRate < 0) > 0
        stmRate(stmRate<0) = eps;
    end

    % Reconstructed Hazard Curve
    Nn = size(Resp, 2);
    Resp = Resp(logical(xs), :);
    Resp(isnan(Resp)) = -Inf;
    tgtThreshs_rec = zeros(size(model.tgtThreshs));
    tgtRates = model.tgtRates;
    parfor ii = 1:Nn
        Resp_temp = Resp(:, ii);
        % sort the response and define hazard statistics
        [~, sorti] = sort(Resp_temp, 'descend');
        HC_rec_x  = Resp_temp(sorti)';
        HC_rec_Pf = cumsum(stmRate(sorti))';
        % address dry isntances
        removei = HC_rec_x == -inf;
        HC_rec_x(removei)  = [];
        HC_rec_Pf(removei) = [];

        Resp_temp(Resp_temp == -inf) = NaN;

        % Thresholds
        [~, tgtThreshs_rec(:, ii), ~] = TargetRatesNThreshs_LargeA_mod(Resp_temp, HC_rec_Pf, HC_rec_x, tgtRates(:, ii), interp, HC_given);
    end
    d = -model.weights_Corr_rates*spatial_corr(model.tgtThreshs, tgtThreshs_rec, model.anchorpts, model.weights_Corr_nodes);
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