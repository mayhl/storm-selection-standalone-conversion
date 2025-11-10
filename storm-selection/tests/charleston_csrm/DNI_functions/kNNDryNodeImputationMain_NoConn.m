function SurgeCorrected=kNNDryNodeImputationMain_NoConn(Resp,grid,index_kept)
% Main file for calibration and implementation of kNN based enhancement of dry nodes.
% Needs following auxiliary files
%
% wetinformation: calculate preliminaries for problem
% kNN_distance: calculate indices and distances of nearest neighbors based
%               only on grid distances
% kNN_calibration: optimizes the kNN interpolation weights
% kNN_correction: performs the kNN correction

% V.1 2/23/2020 Alexandros Taflanidis (a.taflanidis@nd.edu)
% v.1.1 2/8/2021 Alexandros Taflanidis (a.taflanidis@nd.edu) corrected bug
%                related to total_k exceeding the number of available nodes for which 
%                distance has been estimated already, also added swtich to 
%                accomodate cases when we do not care if closest neighbor
%                could be node itself
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Define characteristics for database
    data.DBname=[]; %database name used for convenience in saving files 
    
    % Database information, included in 'data' structure
    data.input_file=[]; % Name of input file inside
                        % Input file needs to include input matrix named Grid.
                        % Rows of this matrix correspond to differens save
                        % points. Columns to lat, long and elevation for each save
                        % location. Note that elevation and not depth is used
                        % here (negative elevation corresponds to positive depth). 
                        % Order of geospatial inputs needs to be Latitude/Longitude
                        % Can leave empty if information is pass directlty to
                        % wetinformation m-file
    
    data.output_file=[]; % name of output file 
                        % Output file needs to include output matrix named Resp
                        % Columns of this matrix correspond to surge for differens storms
                        % For peak respones this is a 2D matrix, with rows
                        % corresponding to different storms.
                        % For time-series this is a 3D matrix, with rows 
                        % corresponding to different times and the third
                        % matrix dimension corresponding to different locations
                        % Can leave empty if information is pass directlty to
                        % wetinformation m-file
    
    data.LatLimits=[]; % Limit for Latitude for domain of interest, leave empty if there is no limit 
    data.LongLimits=[];  % Limit for Longitude for domain of interest, leave empty if there is no limit                    
    
    data.ThreshWet=0;    % wetness-threshold, ignore any location that is not 
                             % wet for at least this ratio of storms. These are
                             % removed from database. Use 1/number of storms to include 
                             % all nodes wet in at least one storm                       
    data.k_max=[50 300];      % largest value of nearest neighbor information to keep in memory. 
                             % first number corresponds to callibration of interpolation 
                             % and second to prediction. This simply faciliates
                             % computational efficiency in the implementation 
    data.zer_dist=1;         % set to 1 if you do not want closest identified 
                             % neigbor to be node itself, else set it equal to zero 
    data.ThreshElev=[-2 -20];% threshold for min elevation of points to be used in calibration of interpolation
                             % avoids using deeper points. If empty all points are
                             % used. If not empty, then first value is threshold for points
                             % that will be predicted and second value is threshold for points
                             % that will be used for prediction (needs to be
                             % smaller than first)
    data.ThreshVar=0.05;     % threshold for variation of max-min within database 
                             % of points to be used in calibration of
                             % interpolation. Avoids using points with negligible variability. 
                             % If empty, all points are used                   
    data.error_thresh=[];  % surge threhold to adjust for erroneous results. Any value below
                             % this threshold is assume to be completely
                             % missing in the original data (neither dry nor
                             % wet). Leave empty if not used
    
    data.connectivity=[];   % leave empty if ADCIRC connectivity will not be used. 
                             % Else this is the number of connected nodes to identify using
                             % ADCIRC grid. Should be larger than the larger data.kmax value. 
                             % Suggested value is to use at least twice
                             % data.k_max(2). 
                             % In this case the elem connectivity
                             % of ADCIRC is needed. This needs to be included
                             % in the input_file as elem variable, 
                             % with rows corresponding to connectivity of each element. 
                             % Can be passed also directly to wetinformation m-file. 
    data.ancillary_surge=0;  % value of ancillary surge to accomodate issues with ADCIRC 
                             % connectivity; when part of grid appears
                             % disconnected and all nodes remain dry, then surge value
                             % will be set to ancillary_surge. Use value that actual surge not 
                             % likely to get so that the surge can be
                             % distinguished. Codes included three different options to 
                             % correct surge for these nodes 
                             
    %% Load up matrices and set up preliminaries for wet nodes
    %load up data and calculate wet info for each save location as well as for each storm
    data=wetinformation(data,grid,Resp,index_kept);
    clear Resp
    
    %plot kept grid and grid that needs to be corrected
    figure
    plot(data.input(:,2),data.input(:,1),'bo',data.input(data.p_stormw<1,2),data.input(data.p_stormw<1,1),'g*')
    legend('grid','points dry in at least one storm')
    %use this plot to adjust the data.LatLimits and data.LongLimits to create
    %bounding box if desired
    
    % New outputs defined for 'data' are
    % data.input is input matrix with geographical lat/long and elevation for 
               % all save locations that satisfy wetness criteria. Input order is lat/long and elevation 
    % data.surge is respective surge output matrix 
    % data.n_nodes is number of kept nodes
    % data.n_storms is number of total storms
    % data.i_wet is index for each "kept" location and each storm with 
                 % information whether it is wet(1) or dry (0)
    % data.i_alwaysw is index with "kept" locations that are always wet
    % data.i_elev1 & 2 is index with "kept" locations that satisfy,
                      % respectively the 1st and 2nd elevation thresholds 
    % data.i_var is index with "kept" locations that satisfy the variation
                % max-min threshold
    % data.p_nodew for each storm percentage of "kept" save locations that are wet
    % data.p_stormw for each "kept" save location percentage of storms for which it remained wet
    % data.p_alwaysw percentage of "kept" save locations that were always wet
    % data.index_kept indices of the points kept from the initial database
    % data.elem conenctivity of elements in ADCIRC grid, updated only if
    %           data.connectivity is not empty
    % data.n_missing number of nodes per storm considered to have missing
    %                       values
    % data.index_missing sparse matrix with the indices of nodes with missing
    %                    values. Rows are the different nodes and columns the
    %                    di
    %% Figure out ADCIRC criterion for storm elevation to be considered as dry or wet 
    water_level=data.surge-data.input(:,3); %calculate water level
    %histogram of points with water elevation between 0 and 0.2
    histogram(water_level(water_level>0&water_level<0.2),20)
    title('points with water elevation between 0 and 0.2')
    % if there is a jump around 0.1 (or somewhere around there) then ADCIRC is
    % using a minimum threshold for water elevation to consider node as wet. 
    % This needs to be accomodated for the decision of misclasification.
    % One approach is to set the node elevation to a +0.1 m (or whatever
    % the threhold is), and understand that for wet points the water level needs to be further
    % adjusted by a +0.1 m (since we modified the node elevation). Another is
    % to adjust the misclassification after it is initially calculated 
    clear water_level;
    elev_adjust=0.1; %adjust this value based on histogram 
    
    %% Calculate distances for calibration
    % select what points to use to predict "to" and "from"
    %interection of always wet points that satisfy the thresholds for variability and elevation 
    data.index_ptWet=intersect(data.i_alwaysw,intersect(data.i_elev1,data.i_var));
    data.index_pfWet=intersect(data.i_alwaysw,intersect(data.i_elev2,data.i_var));
    disp('Calculate Distances for Calibration')
    % calculate distances 
    tic
    % calculate nearest neighbors based on distances
    [index_neighbors,dis_neighbors]=kNN_distance(data.input(data.index_ptWet,1:2),data.input(data.index_pfWet,1:2),data.k_max(1),data.zer_dist);
    index_neighbors=data.index_pfWet(index_neighbors); %convert neighbor indexing to global 
    n_neighbors=data.k_max(1)*ones(size(index_neighbors,1),1); %for all points desired neighbors identified 
    toc
    % augment data matrix and remove info from memory
    data.index_neighborsWet=index_neighbors;
    data.dis_neighborsWet=dis_neighbors; data.n_neighborsWet=n_neighbors;
    clear index_neighbors dis_neighbors n_neighbors;
    
    % New outputs defined for 'data' are 
    % data.index_neighborsWet indexes of data.k_max closest neighbors from closer
    %                      to further away
    % data.dis_neighborsWet distances of data.k_max(1) closest neighbors in ascending order
    % data.n_neighborsWet number of closest neighbors identified if it is
    %                       smaller than data.k_max(1) due to ADCIRC conenctivity problems 
    % data.index_ptWet indexes of wet points to predict to 
    % data.index_pfWet indexes of wet points to predict from 
    
    %% Define characteristics for calibraton and perform calibration and cross-validation  
    % calibration information, included in 'opt' structure
    clear opt
    opt.type='number'; %if optimization is performed based on 'number' of closest neighbors or 'distance' or also optimized 'weights'
    opt.lb=[1 0.2 0.1 0.4];% lower bounds for (1) number of closest neighbors (2) distance of closest neighbors allowed and (3-4) length scale and exponent of weights  
    opt.ub=[3 2.2 2.5 2];% upper bounds for same variables
    opt.method='ga'; %perform optimization with genetic algorithms 'ga' or random search 'random'. For 'number' type of optimization search is default exhaustive between bounds 
    opt.func=2000; %number of trials for random search, or number of maximum function evaluations for ga
    opt.dist_kept=[];   %if empty calibration performed using all points, else 
                         % only points for which the opt.dist_kept(1) closest
                         % neighor distance is at least opt.dist_kept(2) are
                         % used. Avoids performing calibration for points that
                         % close neighbors are far away
    opt.statistics='storms'; % whether to select optimal k using statistics estimated across
                             % 'storms' or 'points'
    opt.criterion='NRMSE'; % what criterion to use for select of optimal k value
                           % can be NRMSE or RMSE (normalized or not across
                           % different points/storms)
    opt.x_cand=[]; %if opt.method is 'random' can set random search vector (else leave empty). Should have random search experiments in rows
    
    aux_storms=1:data.n_storms; % modify this to use smaller number of storms if desired for quicker calibration 
                                % entire set of storms is 1:data.n_storms 
    
    % calculate statistics for the different possible k values and select best k value
    % and calculate detailed statistics for it
    disp('Perform Calibration')
    tic
    [opt,val]=kNN_calibration(opt,data.dis_neighborsWet,data.index_neighborsWet,data.index_ptWet,data.surge(:,aux_storms),data.input(data.index_ptWet,3)+elev_adjust);
    toc
    
    % New outputs defined for 'opt' are
    % opt.x_opt %optimal x for kNN weights
    % opt.f_opt %objective function for opt.x_opt
    % opt.x_cand opt.f % if opt.method is 'random' then these are the candidate experiments and corresponding objective function values 
                        
    % validation results in 'val' structure
    % val.Y, Y_hat      %actual and predicted response for all points and storms 
    % val.AR2 , ACC , ANRMSE , AAE, AMC  average across all storms coefficient of 
                        % determination, correlation coefficient,  
                        % normalized root mean squared error, absolute error,
                        % misclasification percentage
    % val.LR2 , LCC , LNRMSE , LAE, LMC  for each points (so statistics across storms) the coefficient of 
                        % determination, correlation coefficient, 
                        % normalized root mean squared error, absolute error,
                        % misclasification  percentage
    % val.SR2 , SCC , SNRMSE , SAE, SMC  for each storm (so statistics across points) the coefficient of 
                        % determination, correlation coefficient, 
                        % normalized root mean squared error, absolute error, 
                        % misclasification percentage 
    
    %plot locations with "lower" accuracy
    figure
    thresh=0.98; %threshold defining lower accuracy 
    plot(data.input(data.i_alwaysw,2),data.input(data.i_alwaysw,1),'bo',data.input(data.index_ptWet(val.LCC<thresh),2),data.input(data.index_ptWet(val.LCC<thresh),1),'g*')
    legend('always wet grid','points with lower accuracy')
    
    %% Calculate distances for correction
    % select what points to use to predict "to" and "from"
    %interection of always wet points that satisfy the thresholds for variability and elevation 
    data.index_pt=setdiff(1:data.n_nodes,data.i_alwaysw)'; %nodes dry in at least one storm 
    data.index_pf=1:data.n_nodes; %all nodes used  
    
    disp('Calculate Distances for Correction')
    % calculate distances 
    tic
    % calculate nearest neighbors based on distances
    [index_neighbors,dis_neighbors]=kNN_distance(data.input(data.index_pt,1:2),data.input(data.index_pf,1:2),data.k_max(2),data.zer_dist);
    index_neighbors=data.index_pf(index_neighbors); %convert neighbor indexing to global
    n_neighbors=data.k_max(2)*ones(size(index_neighbors,1),1); %for all points desired neighbors identified
    toc
    
    % augment data matrix and remove info from memory
    data.index_neighbors=index_neighbors;
    data.dis_neighbors=dis_neighbors; data.n_neighbors=n_neighbors;
    clear index_neighbors dis_neighbors n_neighbors
    % New outputs defined for 'data' are 
    % data.index_neighbors indexes of data.k_max(2) closest neighbors from closer
    %                      to further away
    % data.dis_neighbors distances of data.k_max(2) closest neighbors in ascending order
    % data.n_neighbors number of closest neighbors identified if it is
    % data.index_pt indexes of points to predict to 
    % data.index_pf indexes of points to precit from 
    
    %% Define characteristics for correction and perform correction   
    % information, included in 'cor' structure
    cor=[]; %initialize correction structure 
    cor.total_k=opt.x_opt(1)+4; % value of total NN that will need to include desired number of wet kNN
                    % for points that are far away from current wet points, this
                    % leads to a progressive correction where points closer to
                    % current wet points are corrected first. Smaller value contributes
                    % to smoother variations
    cor.wet='yes'; % whether to predict ('yes') the wet nodes per storm for validaiton purposes. Set to 'no' if you 
                   % do not want to. Provides some minor computational savings
    cor.mc_thresh=0.05; % if point is missclasified (adjusted surge above elevation even though we know it is wet)
                   % the adjusted surge is set to mc_thresh value below point
                   % elevation 
    cor.x=opt.x_opt; %value of x for kNN weights
    cor.mc='no'; %whether to correct or not missclasified nodes. If 'no' you can correct always later since you know their indices  
    % for points initially mis-clasified perform a second round
                        % round of correction using a higher degree of averaging.
                        % cor.total_kup replaces total_k and x_up replaces x
                        % leave cor.x_up empty if you do not want to perform
                        % second round of correction. Generally second round of
                        % correction is not helpful so can skip. This is retained 
                        % for legacy purposes 
    cor.total_kup=[]; cor.x_up=[]; 
    
    cor_up=cor; %set auxiliary variable to use for updating the correction later
    
    disp('Perform Correction')
    tic
    cor=kNN_correction(data.dis_neighbors,data.index_neighbors,data.surge,data.index_pt,...
        data.i_wet,data.input(:,3)+elev_adjust,opt.type,cor,data.k_max(2),data.n_neighbors);
    toc
    % Updated cor fields are
    % Y_hat response for all the PointsTo nodes
    % iter number of iterations required per storm
    % index_mc the indices of misclassified nodes per storm (this is a
    %          structure)
    % MC_or and MC percentage of misclassified nodes either for original 
    %               correction step or after the second step 
    % SCC, SNRMSE, SAE statistics (correlation coefficient, NRMSE, absolute
    %               error) for the surge prediciton for the wet nodes within each storm
    % ACC, ANRMSE, AAE average statistics across all storms 
    
    
    %% plot misclassified nodes for some storms
    i=1; %select storm 
    figure(1)
    plot(data.input(data.index_pt(cor.index_mc{i}),2),data.input(data.index_pt(cor.index_mc{i}),1),'o',data.input(data.index_pt,2),data.input(data.index_pt,1),'x')
    legend('misclassified nodes','all nodes')
    title('misclassified node location for specific storm')
    
    figure(2)
    plot3(data.input(data.index_pt(cor.index_mc{i}),2),data.input(data.index_pt(cor.index_mc{i}),1),data.input(data.index_pt(cor.index_mc{i}),3),'o',...
        data.input(data.index_pt(cor.index_mc{i}),2),data.input(data.index_pt(cor.index_mc{i}),1),cor.Y_hat(cor.index_mc{i},i),'x')
    legend('elevation','surge')
    title('misclassified node location for specific storm')
    
    figure(3)
    plot3(data.input(data.index_pt(cor.index_mc{i}),2),data.input(data.index_pt(cor.index_mc{i}),1),data.input(data.index_pt(cor.index_mc{i}),3),'o',...
        data.input(data.index_pt,2),data.input(data.index_pt,1),data.input(data.index_pt,3),'x');
    legend('misclassified nodes','all nodes')
    title('misclassified node elevation for specific storm')
    
    %find the points that are misclassified across all storms 
    aux=union(cor.index_mc{1},find(data.i_wet(data.index_pt,1)==1)); %misclassified nodes in first storm (or wet) 
    for i=2:length(cor.index_mc)
        aux=intersect(aux,union(cor.index_mc{i},find(data.i_wet(data.index_pt,i)==1))); %intersection with misclassified storms in all other storms (or wet)
    end
    figure(4)
    plot(data.input(data.index_pt(aux),2),data.input(data.index_pt(aux),1),'o',data.input(data.index_pt,2),data.input(data.index_pt,1),'x')
    legend('always misclassified nodes','all nodes')
    title('misclassified node location across storms')
    always_mis=data.index_pt(aux);
    
    %% recalculate misclasification based on the adjustment of the node elevation if desired 
    elev=data.input(data.index_pt,3);
    for i=1:length(cor.index_mc)
        aux=find(cor.Y_hat(cor.index_mc{i},i)<elev(cor.index_mc{i})+elev_adjust);
        cor.MC(i)=cor.MC(i)*(length(cor.index_mc{i})-length(aux))/length(cor.index_mc{i});
        cor.index_mc{i}(aux)=[];
    end
    cor.AMC=mean(cor.MC);
    
    %% Update the original data set and save. If you need to correct first nodes set to 
    % update the surge values
    SurgeCorrected=data.surge(1:data.n_nodes,:);
    SurgeCorrected(data.index_pt,:)=cor.Y_hat;
    % update the indices for misclassified nodes
    index_mc=zeros(data.n_nodes,size(data.surge,2));
    for i=1:length(cor.index_mc)
        index_mc(data.index_pt(cor.index_mc{i}),i)=1;
    end
    index_mc=sparse(index_mc); %make matrix sparse 
    
    % keep grid and wetness information for nodes 
    grid=data.input;
    p_stormw=data.p_stormw;
    index_kept=data.index_kept;  
    % estimate some misclassification information for kNN
    MC_kNN(1)=cor.AMC; %across all corrected nodes
    MC_kNN(2)=mean(mean(index_mc(p_stormw>=5/data.n_storms,:)))/mean(mean(index_mc))*cor.AMC; %across nodes wet in 5 storms at least
    MC_kNN(3)=mean(mean(index_mc(p_stormw>=10/data.n_storms,:)))/mean(mean(index_mc))*cor.AMC; %across nodes wet in 10 storms at least
    MC_kNN(4)=mean(mean(index_mc(p_stormw>=0.1,:)))/mean(mean(index_mc))*cor.AMC; %across nodes wet in 10% of the original storms

end