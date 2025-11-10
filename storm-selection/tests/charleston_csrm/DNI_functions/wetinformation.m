function data=wetinformation(data,grid,Resp,index_kept,elem_f)
% data structure as defined in kNNDryNodeImputationMain
% Grid: Lat/Long of grid if provided directly (else given in
%         data.input_file)
% Resp: Surge for all grid points and storms if provided directly (else given in 
%       data.output_file)
% index_kept: indices of all kept nodes if provided directly, else
%             calculated within the file 
% elem: ADCIRC connectivity matrix, needed if data.connectivity is not
%      empty
% V 1 2/23/2020 Alexandros Taflanidis (a.taflanidis@nd.edu)
% V 1.1 2/27/2021 Alexandros Taflanidis (a.taflanidis@nd.edu) added option
%                 to allow pre-definition of the indices for the kept points 
%                 within the domain, and also accomodate
%                 ADCIRC connectivity information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% load the input and calculate wet information for save locations and for
% each storm

% load up matrices if data not already provided 
if ~exist('Resp','var'); load(data.output_file); elseif isempty(Resp); load(data.output_file); end
if ~exist('grid','var'); load(data.input_file); elseif isempty(grid); load(data.input_file);end
if ~exist('index_kept','var'); cacl_nodes=1; elseif isempty(index_kept); cacl_nodes=1;else; cacl_nodes=0; end

% at this point the followig matrices are needed
% grid: matrix with grid input
% Resp: matrix with responses
data.input=grid; % input with node spatial geometry (latitude and longitude) and elevation
data.surge=Resp; % output of surge responses

%calculate for each location and storm, index for location being wet (1) or no (0)
data.i_wet=1-isnan(data.surge)-(data.surge<-99900);
p_stormw=mean(data.i_wet,2); % for how many storms each kept save location is wet

if cacl_nodes==1 %check if what nodes to keep needs to be estimated
    % identify save locations that are below desired wetness ratio
    index=find(p_stormw<data.ThreshWet);
    %if Lat/Long limits defined update the index
    if ~isempty(data.LongLimits); index=union(index,find(data.input(:,2)<data.LongLimits(1)|data.input(:,2)>data.LongLimits(2))); end
    if ~isempty(data.LatLimits); index=union(index,find(data.input(:,1)<data.LatLimits(1)|data.input(:,1)>data.LatLimits(2))); end
    %remove these locations
    data.input(index,:)=[];data.surge(index,:)=[]; data.i_wet(index,:)=[];
    %indices of points kept
    data.index_kept=setdiff(1:size(Resp,1),index)';
else %if it has been provided use the nodes provided 
    data.index_kept=index_kept;
    data.input=data.input(index_kept,:);
    data.surge=data.surge(index_kept,:);
    data.i_wet=data.i_wet(index_kept,:);
end

if ~isempty(data.connectivity) %check if ADCIRC connectivity needs to be used
    %if it needs save the ADCIRC elem 
    if exist('elem_f','var'); elem=elem_f; clear elem_f; end %update if elements have been provided directly
    % check if element information exist and either update or provide error
    % message
    if exist('elem','var'); data.elem=elem;
    else
        disp('element informaiton not provided so connectivity cannot be identified')
        data.connectivity=[]; %update connectivity to empty
    end
end

data.n_storms=size(data.surge,2); %number of storms in database
data.n_nodes=size(data.surge,1); %number of save locations in database

%calculate for each kept location and storm, index for location being wet
%(1) or no (0) as well as the index of save locations that are always wet
data.p_stormw=mean(data.i_wet,2); % for how many storms each kept save location is wet
data.p_nodew=mean(data.i_wet,1); % for how many storms each kept save location is wet
data.i_alwaysw=find(data.p_stormw==1); % index for always wet save locations
data.p_alwaysw=length(data.i_alwaysw)/data.n_nodes; %percentage of kept save locations always wet

% points that satisty elevation and variation thresholds
data.i_elev1=find(data.input(:,3)>data.ThreshElev(1));
data.i_elev2=find(data.input(:,3)>data.ThreshElev(2));
if data.n_storms>1 %if more than 1 storms exist
    data.i_var=find(max(data.surge,[],2)-min(data.surge,[],2)>=data.ThreshVar);
else
    data.i_var=find(max(data.surge,[],2)-min(data.surge,[],2)>=0);
end
end