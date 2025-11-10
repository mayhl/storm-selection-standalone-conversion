%% Charleston_CSRM_Processing_INITIAL.m
%% Author: Luke Aucoin
%% Editor: Kohl Morris
%% Created 2025-01-29 to pre-process Resp and ProbMass arrays to account 
%%     for data features such as missing storms and dry locations.
clear; clc; close all;
addpath(genpath('.\DNI_functions')) % Add path to functions

%% GET script folder path
folder_path = char(mfilename('fullpath'));
idx=strfind(folder_path,'/');
folder_path = folder_path(1:idx(end));
input_proc_path = folder_path + "input/processed/";
input_raw_path = folder_path + "input/raw/";
input_final_path = folder_path + "input/final/";


%% OUTPUT folder
%% out_folder = '.\Inputs_Processed';
if ~exist(input_proc_path,'dir'); mkdir(out_folder); end

%% LOAD inputs
% Load Resp arrays
dataWL=readmatrix(input_raw_path + "SACS_SA_TP_SYN_SLC_0_Tides_0_WAV_1_MaxEle_Table_meters_NAVD88_20200508.csv");
Surge_S0T0 = dataWL(2 : end, 2 : end);
Surge_S0T0(Surge_S0T0 < -999) = NaN;

dataWL=readmatrix(input_raw_path + "SACS_SA_TP_SYN_SLC_0_Tides_2_WAV_1_MaxEle_Table_meters_NAVD88_20200508.csv");
Surge_S0T2 = dataWL(2 : end, 2 : end);
Surge_S0T2(Surge_S0T2 < -999) = NaN;

% Load Probability Masses
ProbMass = dload(input_raw_path + "SACS_NCSFL_ITC_ProbMass_600km.mat");
Param_MT = dload(input_raw_path + "SACS_NCSFL_TC_Param_MasterTable.mat");
keep_storms = (1:size(Param_MT,1))';

% Load Save Point coordinates
SPs = dload(input_raw_path + "SACS_NCSFL_staID.mat"); % [ID LAT LON DEPTH]
npts = size(SPs,1);


%% REMOVE missing storms
% Identify missing storms in Resp data
rem = or(sum(isnan(Surge_S0T0),1)==npts , sum(isnan(Surge_S0T2),1)==npts);
rem = rem';

%% REMOVE far away storms
% Define Study site
P1 = [ 33.00,-80.20 ];
P2 = [ 32.50,-80.20 ];
P3 = [ 32.50,-79.60 ];
P4 = [ 33.00,-79.60 ];
POLY = [P1;P2;P3;P4];
dist_thresh = 250; % km
% Identify missing storms in Resp data
[dist_deg,~] = distance(mean(POLY(:,1)),mean(POLY(:,2)),Param_MT(:,4),Param_MT(:,5));
dist_km = deg2km(dist_deg);
out_rad = dist_km > dist_thresh;
rem = find(or(rem,out_rad));

% Remove missing storms
Surge_S0T0(:,rem) = [];
Surge_S0T2(:,rem) = [];
Param_MT(rem,:) = [];
ProbMass(rem,:) = [];
keep_storms(rem) = [];

%% REMOVE dry locations
nstorms = size(Param_MT,1);
dry_thresh = 0.95*nstorms; % % Node is omitted if dry for > 95% of storms
% Identify missing storms in Resp data
rem = or(sum(isnan(Surge_S0T0),2)>dry_thresh , sum(isnan(Surge_S0T2),2)>dry_thresh);
rem = find(rem);

% Remove dry locations
Surge_S0T0(rem,:) = [];
Surge_S0T2(rem,:) = [];
SPs(rem,:) = [];

%% SAVE
% Resp no tides
Resp = Surge_S0T0';
Ind_wet = Resp;
Ind_wet(~isnan(Ind_wet)) = 1; Ind_wet(isnan(Ind_wet)) = 0;
save(fullfile(input_proc_path, "Surge_SLC_0_Tides_0.mat"),'Resp','Ind_wet','keep_storms','-v7.3');
% % % % Resp no tides (dnc)
% % % SPs2=SPs(:,2:4); SPs2(:,end)=-1*SPs2(:,end);
% % % Resp=kNNDryNodeImputationMain_NoConn(Resp',SPs2,(1:length(SPs2))');
% % % Resp=Resp';
% % % save(fullfile(out_folder,'Surge_SLC_0_Tides_0_dnc.mat'),'Resp','keep_storms','-v7.3');

% Resp with tides
Resp = Surge_S0T2';
Ind_wet = Resp;
Ind_wet(~isnan(Ind_wet)) = 1; Ind_wet(isnan(Ind_wet)) = 0;
save(fullfile(input_proc_path,"Surge_SLC_0_Tides_2.mat"),'Resp','keep_storms','Ind_wet','-v7.3');
% % % % Resp no tides (dnc)
% % % SPs2=SPs(:,2:4); SPs2(:,end)=-1*SPs2(:,end);
% % % Resp=kNNDryNodeImputationMain_NoConn(Resp',SPs2,(1:length(SPs2))');
% % % Resp=Resp';
% % % save(fullfile(out_folder,'Surge_SLC_0_Tides_2_dnc.mat'),'Resp','keep_storms','-v7.3');

% Locations
save(fullfile(input_proc_path,'Savepoints.mat'),'SPs','-v7.3');

% Storm Parameters
save(fullfile(input_proc_path,'Storm_Parameters.mat'),'Param_MT','-v7.3');

% Discreet Storm Weights
save(fullfile(input_proc_path,'DSWs.mat'),'ProbMass','-v7.3');

rmpath(genpath('.\DNI_functions')) % Remove path to functions
%% END