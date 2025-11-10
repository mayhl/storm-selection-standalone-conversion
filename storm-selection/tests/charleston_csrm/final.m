%% Charleston_CSRM_Processing_FINAL.m
%% Author: Luke Aucoin
%% Created 2025-02-26 to prepare files for Storm Selection after final 
%%    savepoints are chosen.
clear; clc; close all;


%% GET script folder path
folder_path = char(mfilename('fullpath'));
idx=strfind(folder_path,'/');
folder_path = folder_path(1:idx(end));

%% Temp folders while reorginizing 
input_proc_path = folder_path + "input/processed/";
input_raw_path = folder_path + "input/raw/";
input_final_path = folder_path + "input/final/";

% Input folder
in_folder = input_proc_path;
if ~exist(in_folder,'dir'); mkdir(in_folder); end

%% LOAD and save
% Locations
load(fullfile(in_folder,'Savepoints.mat'),'SPs');
SPs_all = SPs; clear SPs
load(input_final_path+ 'Charleston_CSRM_Subset_Grid.mat','gridSM');
id = gridSM(:,3); clear gridSM;
[~,idx] = intersect(SPs_all(:,1),id);

% % % % Resp no tides
% % % load(fullfile(in_folder,'Surge_SLC_0_Tides_0.mat'),'Resp','Ind_wet','keep_storms')
% % % Resp = Resp(:,idx);
% % % Ind_wet = Ind_wet(:,idx);
% % % save('Charleston_CSRM_Subset_Output_NoTides','Resp','Ind_wet','keep_storms','-v7.3')

% Resp with tides
load(fullfile(in_folder,'Surge_SLC_0_Tides_2.mat'),'Resp','Ind_wet','keep_storms')
Resp = Resp(:,idx);
Ind_wet = Ind_wet(:,idx);
save(input_final_path + 'Charleston_CSRM_Subset_Output_WithTides','Resp','Ind_wet','keep_storms','-v7.3')

% Discreet Storm Weights
load(fullfile(in_folder,'DSWs.mat'),'ProbMass');
save(input_final_path + 'Charleston_CSRM_Subset_ProbMass','ProbMass','-v7.3');

%% END