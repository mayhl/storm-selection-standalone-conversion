function [Resp, opts, HC_target, Ind_wet, gridSM, data] = set_problem_StormSelection(Resp, opts, HC_target, ProbMass, Ind_wet, gridSM)
% set problem for Storm Selection
% - find nodes of interest of which target hazard curves will be used by
%   removing meaningless nodes and clustering if desired
% - make the optimization problem better conditioned by adding intermediate
%   AEPs
%
% -- Inputs:
%  Resp:      surge response matrix of size [number of storms x number of nodes] 
%
%  opts:      a structure with following options. If any field is empty, it
%             will be set to default.
%          <related to data preprocessing>
%            .node_interest: indices of nodes of interest if you want to focuse
%                            on a subset of the domain. default is using all the
%                            nodes (1:number of nodes).
%            .wetness_node:  criterion per node for the probability of being wet. 
%                            Nodes with the probability below than 'opts.wetness_node' 
%                            will be removed. default is 0.
%            .wetness_data_thresh:  criterion for the probability of being wet for 
%                                   the entire data. default is 0.98.
%                                   When it is required to perform PCA/PPCA, 
%                                   if the wetness of the entire data is below than 
%                                   'opts.wetness_data_thresh', PCA is performed on 
%                                   data with dry instances being imputed by 
%                                   the minimum recoreded surge values per each node, 
%                                   otherwise PPCA is performed
%
%          <related to intermediate AEPs>
%            .int_tgtRate.pca:  If 'yes', perform PPCA to find number of
%                               principal components that can explain over 99%
%                               original variance. If 'no' the number will be
%                               set to round(5+45*Nn/100000).
%            .int_tgtRate.sf:   Safety factor to find number of intermediate
%                               AEPs needed. default is 1.5. 
%            .int_tgtRate.wgt:  Relative weight of intermediate AEPs to 
%                               weight of actual target AEPs. default is 0.1.
%            .int_tgtRate.ppca: a structure with options for ppca if
%                               performing ppca is requried
%                             .Nlat: dimension of latent space. default is 
%                                    round(10+90*Nn/100000).
%                             .threshold: if relative change in objective
%                                         function is less than this threshold, 
%                                         it is stopped. default is 1e-5.
%                             .MaxCount: Maximum number of iteration.
%                                        default is 2e3.
%                             .dia: if 1, printf objective each step.
%                                         default is 0 (do not print objective).
%                             .randseed:  if user wants to control random seed, 
%                                         set a number. default is not
%                                         specifying random seed.
%
%          <related to clustering>
%            .cluster.perform:  perform cluster ('yes') or not ('no')
%            .cluster.k:        number of clusters. default is min(300, round(Nn*0.05)).
%            .cluster.type:     perform clustering using geospatial information only ('geo'), 
%                               surge response only ('resp') or both ('comb'). 
%                               default is 'comb' but if gridSM is empty 'resp'.
%            .cluster.norm:     normalize values used for clustering 'yes' or not 'no'. 
%                               default is 'yes'.                
%            .cluster.combwgt:  weight between geospatial information and surge
%                               response if opts.cluster_type is 'comb'. 
%                               default is 1.
%            .cluster.iter:     number of max iterations of clustering algorithm.
%                               default is 1000.
%            .cluster.rep:      number of replicates (if you want to repeat with
%                               different starting points). default is 1.
%            .cluster.avoid:    nodes to preferably avoid as centroids. default is [] (empty).
%            .cluster.avoid_t:  nodes to absolutely avoid as centroids. default is [] (empty).
%
%          <related to storm selection>
%            .tgtRate:     target surge rates.  
%                          default is 1./[10, 20, 50, 100, 200, 500, 1000].
%            .interp:      scale used to interpolate hazard curves to find surge
%                          levels corresponding to target surge rates.
%                          linear interpolation in the original scale (1) or
%                          logscale (2). default is 2.
%          for optimization in general
%            .nselect:         number of storms to select. default is min(50, Ns/10).
%            .storms_included: indices of storms that needs to be included in
%                              the final set. Additional 'opt.nselect' storms
%                              are chosen beyond this set. 
%          for genetic algorithm
%            .objective:    objective function - ('Sdev') Sum of absolute
%                           deviations and ('Scor') spatial correlation between
%                           hazard maps across different target AEPs. default is
%                           'Sdev'.
%            .npop:         population size. default is 300.
%            .generations:  maximum generation. default is 1000.
%          <STOPPING CRITERIA>
%            .FunctionTolerance: The algorithm stops if the average relative change in 
%                                the best fitness function value over 'MaxStallGenerations' 
%                                is less than or equal to FunctionTolerance. 
%                                default is 1e-9.
%            .MaxStallGenerations: The algorithm stops if the average relative change in 
%                                  the best fitness function value over 'MaxStallGenerations' 
%                                  is less than or equal to FunctionTolerance. 
%                                  default is 50.
%            .MaxTime: the maximum time in seconds the genetic algorithm runs
%                      before stopping. default is Inf.
%
%          for weights for objective function of linear programming
%            .weight_rates  set a row vector of weight for target rates of size that 
%                           has to be same as the size of opts.tgtRate. For
%                           rates of anchor points, weights are linearly
%                           interpolated based on 1./opts.tgtRate. 
%                           default is 1./opts.tgtRate*opts.tgtRate(1).
%            .weight_nodes  set a row vector of weight for nodes of size [1 x Nn].
%                           default is ones(1, Nn).
%            .error         calculate objective function considering target rates relatively
%                           (inversly proportional to target rates) with 'rel'
%                           or not with 'abs'. Default is 'rel'.
%           for optimizer (Linear programming)
%            .gurobi        Use gurobi optimizer ('yes') or not ('no').  
%                           If you want use gurobi optimizer, the optimizer 
%                           needs to be installed and a license is required. 
%                           If it is empty, default will be 'yes' if the gurobi optimizer 
%                           is properly setup, or 'no'.
%
%  HC_target: a cell of size [number of nodes x 1] with target hazard curve for 
%             each node. Each cell is a structure with following fields.
%           .x:  surge level thresholds corresponding to probabilities of
%            exceedance (.Pf)
%           .Pf: probabilities of exceedance
%             Note that the lengths of HC_target.x and .Pf have to match.
%             If empty, it will be calculated based on 'Resp' and
%             'ProbMass'.
%
%  ProbMass:  a vector of size [number of storms x 1] with the original 
%             probability masses of the storms
%
%  Ind_wet:   a matrix of size [number of storms x number of nodes] 
%             indicating whether a node is wet (1) or dry (0) for each storm
%
%  gridSM:    a matrix of size [number of nodes x 2] with the latitude and
%             longitude of each node. If empty, grid information will not
%             be used anywhere.
%
% -- Outputs:
%  Resp:      preprocessed surge response matrix.
%
%  opts:      a structure with following updated or new fields.
%            .tgtRate:     [updated] intermediate AEPs are added.  
%            .weight_rates:[updated] weights for intermediate AEPs are
%                                    added.
%            .SS_nodes_use:[new] nodes of interest of which target hazard curves 
%                                will be used to select storms.
%            .HC_target_given:[new] if the target hazard curves are given, 1
%                                   otherwise 2 (Calculated).
%
%  HC_target: preprocessed target hazard curve for each node. 
%
%  Ind_wet:   preprocessed indicator matrix with 1 indicating the node is
%             wet or 0 indicating the node is dry.
%
%  gridSM:    preprocessed grid information (latitude and longitude)
%
%  data:      results of pca/ppca with following fields
%            .pca.C: a matrix with each column being an eigenvector
%            .pca.S: a matrix with latent response
%
% V1    10/24/2023 WoongHee Jung (wjung2@nd.edu)
% V2    02/19/2024 WoongHee Jung (wjung2@nd.edu)
%       1) Modified 'resp' option for cluster analysis by including hazard
%       information (surge thresholds corresponding to target AEPs)
%       2) Modified to consider the number of nodes belonging to each
%       cluster as weight for the representative node in the optimization
%       3) Modified to accommodate a different objective for the outer-loop 
%       optimization such as using spatial correlation across different AEPs 
%       with additional option ("opts.objective"). 
% V2.1  03/22/2024 WoongHee Jung (wjung2@nd.edu)
%       Added an indicator for target hazard curves.
% V2.2  08/28/2024 WoongHee Jung (wjung2@nd.edu)
%       1) Adjusted threshold to find the number of effective locations
%       (principal components), used to determine the number of
%       intermediate AEP levels, from 99% to 95%
%       2) Modified to calculate surge thresholds anyway
%       3) Modified distribution of intermediate target rates to place them
%       in more important sections
% V2.3  09/25/2024 WoongHee Jung (wjung2@nd.edu)
%       1) Fixed to activate the option "opts.int_tgtRate.add"
%       (if 1, intermediate AEPs are added else, not added.)
%       2) Modified to apply dimensionality reduction when Nn >= 10
% V2.4  11/21/2024 WoongHee Jung (wjung2@nd.edu)
%       Fixed a bug when adding intermediate AEPs
% V2.5  12/17/2024 WoongHee Jung (wjung2@nd.edu)
%       Fixed a bug happening after imputation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ns, Nn] = size(Resp);

%% set defaults to empty fields in opts
% Related to preprocessing
if ~isfield(opts, 'node_interest'); opts.node_interest = [];               end
if isempty(opts.node_interest);     opts.node_interest = 1:size(Resp, 2);  end
if ~isfield(opts, 'wetness_node');  opts.wetness_node = 0;                 end
if isempty(opts.wetness_node);      opts.wetness_node = 0;                 end

if ~isfield(opts, 'wetness_data_thresh'); opts.wetness_data_thresh = 0.98; end

if ~isfield(opts, 'int_tgtRate');   opts.int_tgtRate = [];                 end
if ~isfield(opts.int_tgtRate, 'pca')
    if Nn >= 10
        opts.int_tgtRate.pca  = 'yes';      
    else
        opts.int_tgtRate.pca  = 'no';
    end
end
if ~isfield(opts.int_tgtRate, 'sf');   opts.int_tgtRate.sf   = 1.5;        end
if ~isfield(opts.int_tgtRate, 'wgt');  opts.int_tgtRate.wgt  = 0.1;        end

% Related to clustering
if ~isfield(opts, 'cluster');         opts.cluster = [];                end
if ~isfield(opts.cluster, 'perform'); opts.cluster.perform = 'yes';     end
if strcmp(opts.cluster.perform, 'yes')
    if ~isfield(opts.cluster, 'k');    opts.cluster.k = min(300, round(Nn*0.05));  end
    if ~isfield(opts.cluster, 'type')
        if (~exist('gridSM', 'var')||isempty(gridSM)); opts.cluster.type = 'resp';
        else; opts.cluster.type = 'comb';           end
    end
    if ~isfield(opts.cluster, 'norm');    opts.cluster.norm = 'yes';    end
    if ~isfield(opts.cluster, 'combwgt'); opts.cluster.combwgt = 1;     end
    if ~isfield(opts.cluster, 'iter');    opts.cluster.iter    = 1000;  end
    if ~isfield(opts.cluster, 'rep');     opts.cluster.rep     = 1;     end
    if ~isfield(opts.cluster, 'avoid');   opts.cluster.avoid   = [];    end
    if ~isfield(opts.cluster, 'avoid_t'); opts.cluster.avoid_t = [];    end
end

% Related to storm selection
if ~isfield(opts, 'tgtRate');    opts.tgtRate = 1./[10, 20, 50, 100, 200, 500, 1000]; end
if ~isfield(opts, 'interp');     opts.interp  = 2;               end
if ~isfield(opts, 'nselect');    opts.nselect = min(50, Ns/10);  end
if ~isfield(opts, 'storms_included'); opts.storms_included = []; end
if ~isfield(opts, 'objective');    opts.objective = 'Sdev'; end
if ~isfield(opts, 'npop');         opts.npop = 300;         end
if ~isfield(opts, 'generations');  opts.generations = 1000; end
if ~isfield(opts, 'FunctionTolerance');   opts.FunctionTolerance = 1e-9; end
if ~isfield(opts, 'MaxStallGenerations'); opts.MaxStallGenerations = 50; end
if ~isfield(opts, 'MaxTime');             opts.MaxTime = Inf;            end
if ~isfield(opts, 'weight_rates'); opts.weight_rates = ones(1, length(opts.tgtRate)); end
if isempty(opts.weight_rates);     opts.weight_rates = ones(1, length(opts.tgtRate)); end
if iscolumn(opts.weight_rates);    opts.weight_rates = opts.weight_rates';            end
if ~isfield(opts, 'weight_nodes'); opts.weight_nodes = ones(1, Nn);                   end
if isempty(opts.weight_nodes);     opts.weight_nodes = ones(1, Nn);                   end
if iscolumn(opts.weight_nodes);    opts.weight_rates = opts.weight_nodes';            end
if ~isfield(opts, 'error');        opts.error = 'rel';           end
if isempty(opts.error);            opts.error = 'rel';           end
aux = which('gurobi_setup');
if ~isfield(opts, 'gurobi')
    if isempty(aux); opts.gurobi = 'no';
    else;            opts.gurobi = 'yes'; end
end
if isempty(aux); opts.gurobi = 'no'; end

%% (0) Define 'Ind_wet' if not provided
if ~exist('Ind_wet', 'var') || isempty(Ind_wet)
    Ind_wet=ones(size(Resp)); % initialize everything wet
    % make dry all Resp instances that are NaN or -99999
    Ind_wet(Resp==-99999)=0; Ind_wet(Resp==-inf)=0; Ind_wet(isnan(Resp))=0; 
end

%% (1) Preprocessing to update the desired nodes within domain
% Nodes of interest
Resp    = Resp(:,    opts.node_interest);
Ind_wet = Ind_wet(:, opts.node_interest);
opts.weight_nodes = opts.weight_nodes(opts.node_interest);
if ~isempty(gridSM); gridSM  = gridSM(opts.node_interest, :); end

% Remove meaningless (almost dry) nodes depending on wetness criterion
probwet = mean(Ind_wet);                  % probabilities of being wet per node
remove_nodes = find(probwet < opts.wetness_node);  % Nodes not satisfying the criterion

% Remove nodes
Resp(:, remove_nodes)    = [];
Ind_wet(:, remove_nodes) = [];
opts.weight_nodes(remove_nodes) = [];
if ~isempty(gridSM); gridSM(remove_nodes, :)  = []; end

Nn = size(Resp, 2);

% necessary if original hazard curve will be defined based on 'Resp'
Resp(~Ind_wet) = -inf;

% Calculate wetness of data:
% Depending on the wetness of the data (# of wet instances / # of all
% instances), it is determined to run PCA or PPCA whenever they need to be performed. 
% If the wetness is greater than a threshold, PCA is performed, otherwise 
% PPCA is performed. The threshold can be set as
% opts.wetness_data_thresh = 0.98; but the default value is 0.98 (98% wet).
opts.wetness_data = sum(Ind_wet, 'all') / numel(Ind_wet);

%% (2) Define the original hazard curve if not provided 
% 'HC_target' is a cell of size [number of nodes x 1] with each element being
% a structure with following fields.
%           .x:  surge level thresholds corresponding to 'HC_target.Pf'
%           .Pf: probabilities of exceedance
% Note that the lengths of HC_target.x and .Pf have to match.
if isempty(HC_target)
    if ~exist('ProbMass', 'var') || isempty(ProbMass)
        error('The original probability masses of the storms (ProbMass) is necessary.')
    end
    disp('Calculating HC_target ...')
    HC_target = cell(Nn, 1);
    for i=1:Nn
        % sort the response and define hazard statistics
        [~, sorti] = sort(Resp(:, i), 'descend');
        HC_target{i}.x  = Resp(sorti, i)';
        HC_target{i}.Pf = cumsum(ProbMass(sorti))';
        % address dry isntances 
        removei = HC_target{i}.x == -inf;
        HC_target{i}.x(removei)  = [];
        HC_target{i}.Pf(removei) = [];
    end
    clear sorti removei
    Resp(Resp == -inf) = NaN;
    opts.HC_target_given = 2; % Switch for HC_target
else
    % lengths of HC_target.x and .Pf have to match.
    if length(HC_target{1}.x) ~= length(HC_target{1}.Pf)
        error('In target hazard curve data, the lengths of probabilities of exceedance and surge level thresholds do not match.')
    end
    HC_target = HC_target(opts.node_interest);
    HC_target(remove_nodes) = []; % If 'HC_target is provided by input, 
                                  % remove nodes that do not satisfy the criterion of 
                                  % the probability of being wet.

    % Remove NaNs in target thresholds
    for ii = 1:length(HC_target)
        temp = isnan(HC_target{ii}.x);
        HC_target{ii}.x(temp) = [];
        HC_target{ii}.Pf(temp) = [];
    end
    Resp(Resp==-99999) = NaN; Resp(Resp == -inf) = NaN;
    opts.HC_target_given = 1;  % Switch for HC_target
end

%% (2-1) Calculate surge thresholds
% re-define variables for using parfor (to avoid large broadcast variables)
tgtRate = opts.tgtRate;
interp  = opts.interp;
HC_target_given = opts.HC_target_given;
% initialization
Nt = length(tgtRate);
tgtThreshs = zeros(Nt, Nn); tgtRates = tgtThreshs; % matrices with each column being target thresholds or rates for each node of interest
parfor i = 1:Nn
    [tgtRates(:, i), tgtThreshs(:, i), ~] = TargetRatesNThreshs_LargeA_mod(Resp(:, i), HC_target{i}.Pf, HC_target{i}.x, tgtRate, interp, HC_target_given);
end
anchorpts  = abs(tgtRates - repmat(tgtRate', 1, Nn)) > eps;
tgtThreshs(anchorpts) = NaN;

%% (3) Perform PCA / PPCA
% Whenever PCA or PPCA has to be performed, it is determined to perform PCA
% or PPCA depending on the wetness of the data.
data = [];
if strcmp(opts.int_tgtRate.pca, 'no')&&(strcmp(opts.cluster.perform, 'no')||(strcmp(opts.cluster.type, 'geo'))) % Not performing PCA / PPCA
    Nindep = round(5+45*Nn/100000); % Number of independent equations
else % Perform PCA / PPCA
    if opts.wetness_data > opts.wetness_data_thresh
        disp('Performing PCA ...')
        % Imputes dry instances with minimum recorded surge
        dryinsts = isnan(Resp);
        Resp(dryinsts) = repelem(min(Resp, [], 'omitnan'), sum(dryinsts));
        [data.pca.C, data.pca.S] = pca(Resp);
        Resp(dryinsts) = NaN;
    else
        disp('Performing PPCA ...')
        if ~isfield(opts.int_tgtRate, 'ppca'); opts.int_tgtRate.ppca = [];                     end
        if ~isfield(opts.int_tgtRate.ppca, 'Nlat');      opts.int_tgtRate.ppca.Nlat = round(10+90*Nn/100000); end
        if ~isfield(opts.int_tgtRate.ppca, 'threshold'); opts.int_tgtRate.ppca.threshold = 1e-5; end
        if ~isfield(opts.int_tgtRate.ppca, 'MaxCount');  opts.int_tgtRate.ppca.MaxCount  = 2e3;  end
        if ~isfield(opts.int_tgtRate.ppca, 'dia');       opts.int_tgtRate.ppca.dia       = 0;    end
        [data.pca.C, ~, ~, data.pca.S] = ppca_mv(Resp, opts.int_tgtRate.ppca);
    end
    L = var(data.pca.S);
    L = cumsum(L)/sum(L);
    Nindep = find(L > 0.95, 1);  % Find the number of principal components corresponding to
    % proportion of explained original variance greater than 95%
end

%%%% WARNING MESSAGE %%%%
if Nindep <= 2
    fprintf('\nWARNING: Locations are very highly correlated.\n\n')
end

%% (4) Add intermediate AEPs
% If needed, some additional intermediate AEPs with small weight
% After running this, 'opts.tgtRate' and 'opts.weight_rates' will be
% changed (including intermediate AEPs)
if strcmp(opts.int_tgtRate.add, 'yes')
    NREq     = (opts.nselect + length(opts.storms_included))*opts.int_tgtRate.sf;  % number of required independent equations
    NtgtRate = length(opts.tgtRate);  % number of original target rates

    Nint_tot = ceil(NREq/Nindep) - NtgtRate;  % number of additional intermediate AEPs

    if Nint_tot > 0
        % calculate difference between target rates
        if opts.interp == 1 % linear interpolation
            d_Pf  = opts.tgtRate(2:end) - opts.tgtRate(1:end-1);
        elseif opts.interp == 2 % interpolation in logscale
            d_Pf  = log(opts.tgtRate(2:end)) - log(opts.tgtRate(1:end-1));
        end
        % Find more informative sections
        temp = mean(tgtThreshs, 2, 'omitnan')';
        temp = abs(d_Pf./(temp(2:end) - temp(1:end-1)));

        % distribute intermediate target rates based on the calculated
        % difference
        weight = temp/sum(temp);
        Nints  = round(weight*Nint_tot);

        % augment uniformly distributed intermediate rates to original rates
        mulfc  = zeros(1, sum(Nints));
        temp = 1;
        for ii = 1:NtgtRate-1
            mulfc(temp:temp+Nints(ii)-1) = 1:Nints(ii);
            temp = temp + Nints(ii);
        end

        if opts.interp == 1 % linear interpolation
            int_tgtRate  = repelem(opts.tgtRate(1:end-1), Nints) + repelem(d_Pf./(Nints+1), Nints).*mulfc;
            opts.tgtRate = sort([opts.tgtRate, int_tgtRate], 'descend');
        elseif opts.interp == 2 % interpolation in logscale
            int_tgtRate  = repelem(log(opts.tgtRate(1:end-1)), Nints) + repelem(d_Pf./(Nints+1), Nints).*mulfc;
            opts.tgtRate = exp(sort([log(opts.tgtRate), int_tgtRate], 'descend'));
        end

        % adjust weight for rates
        weight_rates = zeros(1, length(opts.tgtRate));
        aux1 = 1:NtgtRate;
        aux2 = [0 cumsum(Nints)];
        index_ori_tgtRate = aux1+aux2;
        weight_rates(index_ori_tgtRate) = opts.weight_rates;
        weight_rates(setdiff(1:NtgtRate+sum(Nints), index_ori_tgtRate)) = opts.int_tgtRate.wgt*repelem((opts.weight_rates(1:end-1)+opts.weight_rates(2:end))/2, Nints);
        opts.weight_rates = weight_rates;
    else
        Nints = 0;
    end
    fprintf('In total %d intermediate AEPs are added. \nCheck added intermediate AEPs and their weights. \n', sum(Nints));
end

%% (5) Clustering
if strcmp(opts.cluster.perform, 'yes')&&(opts.nselect > 0)
    disp('Performing clustering analysis')
    % Obtain data matrix that is used to calculate distance between nodes 
    % to perform clustering
    data = data_cluster_analysis(Resp, tgtThreshs, opts, gridSM, data);

    % Perform clustering
    [opts.cluster.cen_id, opts.cluster.member_id,~] = perform_cluster(opts.cluster, data.clustering');
    % 'clust.cen_id' is a vector of size [clust.k x 1] with indices of
    %                   representative nodes we want to use for the optimization
    % 'clust.member_id' is a vector of size [Nn x 1] with the cluster
    %                      each node belongs to

    opts.SS_nodes_use = opts.cluster.cen_id;
    opts.cluster.weight_nodes = zeros(1, opts.cluster.k);
    for ii = 1:opts.cluster.k
        opts.cluster.weight_nodes(ii) = sum(opts.weight_nodes(opts.cluster.member_id == ii));
    end
    data = rmfield(data, 'clustering');
    fprintf('%d representative nodes are selected by clustering analysis. \n', opts.cluster.k);
else
    opts.SS_nodes_use = 1:size(Resp, 2);
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

function [C, ss, M, X, Ye] = ppca_mv(Ye, opts)
%
% implements probabilistic PCA for data with missing values, 
% using a factorizing distribution over hidden states and hidden observations.
%
%  - The entries in Ye that equal NaN are assumed to be missing. - 
%
% [C, ss, M, X, Ye ] = ppca_mv(Y,d,dia)
%
% Y   (N by D)  N data vectors
% d   (scalar)  dimension of latent space
% dia (binary)  if 1: printf objective each step
%
% ss  (scalar)  isotropic variance outside subspace
% C   (D by d)  C*C' +I*ss is covariance model, C has scaled principal directions as cols.
% M   (D by 1)  data mean
% X   (N by d)  expected states
% Ye  (N by D)  expected complete observations (interesting if some data is missing)
%
% J.J. Verbeek, 2006. http://lear.inrialpes.fr/~verbeek
%
% 09/29/2023 WoongHee Jung (wjung2@nd.edu)
%            Modified to accelerate the computation (change inv() --> '/',
%            remove unnecessary for-loop) and added waitbar to see the
%            progress
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N, D]      = size(Ye);       % N observations in D dimensions
threshold   = opts.threshold; % minimal relative change in objective funciton to continue
MaxCount    = opts.MaxCount;  % Maximum number of iteration
hidden      = isnan(Ye); 
missing     = sum(hidden(:));
d           = opts.Nlat;
dia         = opts.dia;

% compute data mean and center data
if missing 
    M    = mean(Ye, 'omitnan');
else
    M    = mean(Ye);                 
end
Ye = Ye - repmat(M,N,1);

if missing
    Ye(hidden)=0; 
end

% =======     Initialization    ======
if isfield(opts, 'randseed')
    rng(opts.randseed); % specifying random seed if desired.
end
C     = randn(D,d);
CtC   = C'*C;
X     = (Ye * C) / CtC;
recon = X*C'; recon(hidden) = 0;
ss    = sum(sum((recon-Ye).^2)) / (N*D-missing);

count = 1; 
old   = Inf;
wb = PoolWaitbar(MaxCount,'Performing PPCA');
while count <= MaxCount          %  ============ EM iterations  ==========      
   
    Sx = eye(d)/( eye(d) + CtC/ss );    % ====== E-step, (co)variances   =====
    ss_old = ss;
    if missing
        proj = X*C'; 
        Ye(hidden) = proj(hidden); 
    end  
    X = Ye*C*(Sx/ss);          % ==== E step: expected values  ==== 
    
    SumXtX = X'*X;                              % ======= M-step =====
    C      = (Ye'*X)  / (SumXtX + N*Sx);    
    CtC    = C'*C;
    ss     = ( sum(sum( (X*C'-Ye).^2 )) + N*sum(sum(CtC.*Sx)) + missing*ss_old ) /(N*D);
    
    objective = N*D + N*(D*log(ss) +trace(Sx)-log(det(Sx)+eps) ) +trace(SumXtX) -missing*log(ss_old);   
           
    rel_ch    = abs( 1 - objective / old );
    old       = objective;
    
    count = count + 1;
    if ( rel_ch < threshold) && (count > 5); break; end
    if dia; fprintf('Objective:  %.2f    relative change: %.7f \n',objective, rel_ch ); end
    increment(wb)
end             %  ============ EM iterations  ==========
delete(wb)

C = orth(C,eps);
[vecs,vals] = eig(cov(Ye*C));
[~,ord] = sort(diag(vals),'descend');
vecs = vecs(:,ord);

C = C*vecs;
X = Ye*C;
 
% add data mean to expected complete data
Ye = Ye + repmat(M,N,1);
end


function [data] = data_cluster_analysis(Resp, tgtThreshs, opts, gridSM, data)
% function to make data matrix for performing cluster. This function has to
% be run after 'set_param_cluster'
% -- Inputs:
%  Resp:     surge response matrix of size [number of storms x number of nodes]  
%
%  tgtThreshs:  target surge thresholds corresponding to target AEPs 
%               obtained from target hazards
%
%  opts:     a structure with options for clustering with following fields.
%            These fields has to be set by running 'set_param_cluster'.
%       .k:     number of clusters. default is min(300, size(Resp, 2)/10).
%       .type:  perform clustering using geospatial information only ('geo'), 
%               surge response only ('resp') or both ('comb'). 
%               default is 'comb' but if gridSM is empty 'resp'.
%       .norm:  normalize values used for clustering 'yes' or not 'no'. 
%               default is 'yes'.                
%       .combwgt:  weight between geospatial information and surge
%                  response if opts.cluster_type is 'comb'. default
%                  is 1.
%       .iter:  number of max iterations of clustering algorithm.
%               default is 1000.
%       .rep:   number of replicates (if you want to repeat with
%               different starting points). default is 1.
%       .avoid:  nodes to preferably avoid as centroids. default is [] (empty).
%       .avoid_t:  nodes to absolutely avoid as centroids. default is [] (empty).
%
%  gridSM:   a matrix of size [number of nodes x 2] with the latitude and longitude
%            of the grid
%
%  data:     a structure with following fields that might exist
%       .pca.C: matrix of which each column is eigenvector.
%       .pca.S: latent response matrix
%
% -- Outputs:
%  data:     a structure with a following added field
%       .clustering: a data matrix used in the clustering analysis
%
%  V1   10/12/2023 WoongHee Jung (wjung2@nd.edu)                  
%  V2   02/19/2024 WOongHee Jung (wjung2@nd.edu)
%       modified to accommodate target surge thresholds in cluster analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define data matrix
if (strcmp(opts.cluster.type, 'resp')||strcmp(opts.cluster.type, 'comb'))
    if ~isfield(data, 'pca')
        if opts.wetness_data > opts.wetness_data_thresh % Perform PCA
            % Imputes dry instances with minimum recorded surge
            dryinsts = isnan(Resp);
            Resp(dryinsts) = repelem(min(Resp, [], 'omitnan'), sum(dryinsts));
            [data.pca.C, data.S] = pca(Resp);

        else % perform PPCA
            if ~isfield(opts.int_tgtRate, 'Nlat'); opts.int_tgtRate.Nlat = round(10+90*Nn/100000); end
            [data.pca.C, ~, ~, data.S] = ppca_mv(Resp, opts.int_tgtRate.Nlat, 0);
        end
    end
    L = var(data.pca.S)';
    L = cumsum(L)/sum(L);
    Npc = find(L > 0.99, 1);  % Find the number of principal components corresponding to
    % proportion of explained original variance greater than 99%

    if isempty(Npc); Npc = length(L); end
    
    % normalization
    Nt = size(tgtThreshs, 1); % number of target AEPs (surge thresholds)
    if strcmp(opts.cluster.norm, 'yes')
        tgtThreshs = tgtThreshs';
        resp_data=[((tgtThreshs - mean(tgtThreshs, 'omitnan'))./std(tgtThreshs, 'omitnan'))';
                   zscore(data.pca.C(:, 1:Npc)',0,2)];
    else
        resp_data=[tgtThreshs; data.pca.C(:, 1:Npc)'];
    end
    weight_resp = [sqrt(sum(L(1:Npc))/Nt)*ones(Nt, 1); sqrt(L(1:Npc))];
    resp_data = resp_data.*repmat(weight_resp, 1, size(resp_data, 2));
end
if strcmp(opts.cluster.type, 'geo')||strcmp(opts.cluster.type, 'comb')
    % normalization
    if strcmp(opts.cluster.norm, 'yes')
        geo_data=zscore(gridSM(:, 1:2)',0,2);
    else
        geo_data=gridSM(:, 1:2)';
    end
    if strcmp(opts.cluster.type, 'comb')
        weight_geo = sqrt(sum(weight_resp.^2))/opts.cluster.combwgt*ones(2,1); 
    else
        weight_geo = [1;1];
    end
    geo_data=geo_data.*repmat(weight_geo,1,size(geo_data,2));
end

if strcmp(opts.cluster.type,'comb')
    data.clustering=[geo_data; resp_data];
elseif strcmp(opts.cluster.type,'geo')
    data.clustering=geo_data;
else
    data.clustering=resp_data;
end
end

function [cen_id,id,cen_id_full]=perform_cluster(cluster,Z)
% perform clustering analysis for matrix Z using characteristic in cluster
% data strucxtures
% cen_id is index of point closest to centroid for each cluster when
%        points in cluster.avoid are avoided
% member_id is for each point the cluster index for its membership
% cen_id_full is index of point closest to centroid for each cluster

tot=Inf; %dummy variable for min value of objective function 
%repeat for different starting points of k-means
if strcmp(cluster.type, 'resp')||strcmp(cluster.type, 'comb') % k-POD clustering
    missing = isnan(Z);
    mean_Z = mean(Z, 'omitnan');
    Z(missing) = 0;
    [Nz, Nv] = size(Z);
    for h=1:cluster.rep
        Z_temp = (~missing).*Z + missing.*(repmat(mean_Z, Nz, 1));
        [idtrial,C_old,sumdtrial,dtrial] = kmeans(Z_temp,cluster.k,'Display','off','MaxIter',cluster.iter); %perform kmeans
        if sum(sumdtrial)<tot %update best solution if betetr than previous best
            id = idtrial; d = dtrial; %keep in memory best solution
            tot=sum(sumdtrial); %update best objective function value
        end
        % idmis_new = randsample(1:ncluster, Nz, true);
        iter = 1;
        while iter <= 1000
            iter  = iter + 1;
            A = sparse(1:Nz, idtrial, ones(1, Nz));
            B = zeros(cluster.k, Nv);
            parfor ii = 1:cluster.k
                B(ii, :) = mean(Z(idtrial == ii, :), 'omitnan');
            end
            Z_temp = (~missing).*Z + missing.*(A*B);
            [idtrial,C_new,sumdtrial,dtrial] = kmeans(Z_temp,cluster.k,'Display','off','Start',C_old,'MaxIter',1000); %perform kmeans

            if sum(sumdtrial)<tot %update best solution if betetr than previous best
                id = idtrial; d = dtrial; %keep in memory best solution
                tot=sum(sumdtrial); %update best objective function value
            end

            if sum((C_new-C_old).^2, 'all') < 1e-3
                break;
            end
            C_old = C_new;
        end
    end
else
    for h=1:cluster.rep
        [idtrial,~,sumdtrial,dtrial] = kmeans(Z,cluster.k,'Display','off','MaxIter',cluster.iter); %perform kmeans
        if sum(sumdtrial)<tot %update best solution if betetr than previous best
            id= idtrial; d = dtrial; %keep in memory best solution
            tot=sum(sumdtrial); %update best objective function value
        end
        clear idtrial sumdtrial dtrial %clear unecessary quantrities from memory
    end
end


%calculate closest point to centroid 
cen_id=zeros(cluster.k,1); %initialize variable 
cen_id_full=cen_id;
for i=1:cluster.k %repeat for each cluster
    % perfomd identification without avoiding any points
    aux=find(id==i); %points that belong in the cluster
    [~, ind_cl]=min(d(aux,i)); %find minimum distance from centroid of cluster of points that belong in the cluster
    cen_id_full(i)=aux(ind_cl); %point with closest distance from centroid
    % perform identification by avoiding cluster.avoid points
    aux=setdiff(find(id==i),cluster.avoid); %points that belong in the cluster and are not included in the "avoided" points
    if ~isempty(aux)
        [~, ind_cl]=min(d(aux,i)); %find minimum distance from centroid of cluster of points that belong in the cluster
        cen_id(i)=aux(ind_cl); %point with closest distance from centroid
    else
        aux_t=setdiff(find(id==i),cluster.avoid_t); %points that belong in the cluster and are not included in the "absolutely avoided" points
        if ~isempty(aux_t)
            [~, ind_cl]=min(d(aux_t,i)); %find minimum distance from centroid of cluster of points that belong in the cluster
            cen_id(i)=aux_t(ind_cl); %point with closest distance from centroid 
        end
    end
end
%now find closest point to centroid for the clusters that do not have yet
%a representative point
aux_t=setdiff(1:length(id),union(cen_id,cluster.avoid_t)); %points to select from (existing cluster centroids removed)
aux_id=find(cen_id==0); %clusters to find representative point 
for i=1:length(aux_id)
    [~, ind_cl]=min(d(aux_t,aux_id(i))); %find minimum distance from centroid of cluster
    cen_id(aux_id(i))=aux_t(ind_cl); %point with closest distance from centroid
    aux_t(ind_cl)=[];%remove point from the points to select from 
end
end
