function [index_neighbors,dis_neighbors]=kNN_distance(input_t,input_f,k,zer_dist,par)
%estimate nearest neighor indices and corresponding distances for points with 
% input input_t (Lat/Long) from group with input input_f (Lat/Lon)
% --Inputs:
%   input_f: latlon of all origin (from) point [lat lon]
%   input_t: latlon of multiple destination (to) points [lat lon]
%   k: desired number of nearest neighbors 
%   zer_dist: index defining whether you want (1) or not (0) closest neigbor to be node itself
%   par: examine (1) or not (0) parallel implementation for efficiency
%   V.1 2/23/2020 Alexandros Taflanidis (a.taflanidis@nd.edu)
%   V.2 2/3/2021  Alexandros Taflanidis (a.taflanidis@nd.edu)
%                 accomodate cases when we do not care if closest neighbor
%                 could be node itself 
%  V.2.1 3/6/2021 Alexandros Taflanidis (a.taflanidis@nd.edu)
%                 accomodate identification for small number of nodes by doing
%                 only for loop (no parallel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('par','var'); par=1; end % set value of par=1 (examine parfor implementation if not defined)

%initialize variables
index_neighbors=zeros(size(input_t,1),k);
dis_neighbors=zeros(size(input_t,1),k);

%calculate distance for each origin point using for loop 
%check first whether it is faster to perform analysis without or with parallel
%computations.Perform first 500 points iteration both ways. auxfunc is the function that calculates nearest neighbors indices and corresponding distances 
tic
for i=1:min(1000,size(input_t,1)); [index_neighbors(i,:), dis_neighbors(i,:)]=auxfunc(input_t(i,:),input_f,k,zer_dist); end
t1=toc;

if size(input_t,1)>1000 %if many points exist, examine if parallel implementaiton beneficial
    if par==1 %examine parallel loop
        %start parfor loop
        delete(gcp('nocreate')); gcp;
        tic
        parfor i=1:1000;[index_neighbors(i,:), dis_neighbors(i,:)]=auxfunc(input_t(i,:),input_f,k,zer_dist);end
        t2=toc;
    else 
        t2=2*t1; %set t2 largger than t1  
    end
    
    %now select fastest approach and perform remaining iterations
    if t2>t1 %faster to execute without parallel computations
        disp(['estimated time to complete distance calculation is ' num2str(t1/60/1000*size(input_t,1)) ' minutes']);
        for i=1000:size(input_t,1)
            [index_neighbors(i,:), dis_neighbors(i,:)]=auxfunc(input_t(i,1:2),input_f(:,1:2),k,zer_dist);
        end
    else %faster to execute with parallel computations
        if par==1; disp(['estimated time to complete distance calculation is ' num2str(t2/60/1000*size(input_t,1)) ' minutes']); end
        parfor i=1000:size(input_t,1)
            [index_neighbors(i,:), dis_neighbors(i,:)]=auxfunc(input_t(i,:),input_f,k,zer_dist);
        end
    end
   if par==1; delete(gcp('nocreate')); end %close parfor loop if needed
end
end

function [index_neighbors, dis_neighbors]=auxfunc(input_t,input_f,k,zer_dist)
%function to estimate nearest neighor indices and corresponding distances
%for point with input input_t (Lat/Long) within group with input input_f
%(Lat/Lon)
% --Inputs:
%   input_t: latlon of origin point [lat lon]
%   input_f: latlon of multiple destination points [lat lon]
%   zero_dist: indicator for whether (1) or not (0) exclude node itself
%   from closest neighbors
%--------------------------------------------------------------------------
%reduce search domain within destination database, select 5*k points with smallest sum of lat+long distance from the origin point 
[~,aux]=mink(sum(abs(input_t-input_f),2),4*k); % aux is indices of the smaller search domain
dkm=lldistkm(input_t,input_f(aux,:)); %caluclate distance of origin point to all other points in the smaller search domain
[oi,oj]=mink(dkm,k+1); %identify smaller k +1 distance values
if oi(1)==0 && zer_dist==1 %check if closest neighbor is node itself and you do not want it to be included 
    index_neighbors=aux(oj(2:end)); %indices of k nearest neighbors
    dis_neighbors=oi(2:end);%distances of k nearest neighbors
else
    index_neighbors=aux(oj(1:end-1)); %indices of k nearest neighbors
    dis_neighbors=oi(1:end-1);%distances of k nearest neighbors
end
end

function dkm=lldistkm(latlon1,latlon2)
% dkm: distance in km based on Pythagoras theorem
% --Inputs:
%   latlon1: latlon of origin point [lat lon]
%   latlon2: latlon of multiple destination points [lat lon]
%--------------------------------------------------------------------------
radius=6371*pi/180;
x=(latlon2(:,2)-latlon1(2)).*cos((latlon2(:,1)+latlon1(1))*pi/360);
y=latlon2(:,1)-latlon1(1);
dkm=radius*sqrt(x.*x + y.*y); %Pythagoran distance
end

