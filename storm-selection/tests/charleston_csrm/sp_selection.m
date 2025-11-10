%% Charleston_CSRM_SP_Selection.m
%% Author: Luke Aucoin
%% Editor: Kohl Morris
%% Created 2025-03-03 for performing storm selection for the Charleston CSRM
clear; clc; close all; rng(17);


%% GET script folder path
folder_path = char(mfilename('fullpath'));
idx=strfind(folder_path,'/');
folder_path = folder_path(1:idx(end));
input_proc_path = folder_path + "input/processed/";
input_raw_path = folder_path + "input/raw/";
input_final_path = folder_path + "input/final/";

%% OUTPUT folder for plots
plot_folder = folder_path + "plots/";
if ~exist(plot_folder,'dir'); mkdir(plot_folder); end

%% LOAD Savepoints
SPs_all = dload(input_proc_path + "Savepoints.mat");
ID = SPs_all(:,1); LAT = SPs_all(:,2); LON = SPs_all(:,3); DEPTH = SPs_all(:,4);
clear SPs_all

%% SELECT savepoints in Raritan Bay, NJ region
P1 = [ 33.00,-80.20 ];
P2 = [ 32.50,-80.20 ];
P3 = [ 32.50,-79.60 ];
P4 = [ 33.00,-79.60 ];
POLY = [P1;P2;P3;P4;P1];

grid = [LAT	LON	ID];
bbidx = inpolygon(grid(:, 2),grid(:, 1),POLY(:,2),POLY(:,1));
grid2 = grid(bbidx, :);

% % % %% CLUSTER execution
% % % K = 30;
% % % cidx = kmeans(Resp2,K);stt
% % % sample = [];
% % % for n = 1:Kccccc
% % %     tmp = datasample(find(cidx == n),1);
% % %     sample = [sample;tmp];
% % % end
%% LUKE AND KAZI POINTS IN BOUNDING BOX - LAA 2025-03-03
M = readmatrix(input_raw_path + "Sponsor_Selected_2025-03-03.csv");
selected_IDs = [M(:,3)' , ...
            13017, 12779, 11237,... % Offshore nodes
            11759, 10873, 11518]; % Onshore nodes
[~ , sample] = intersect(grid2(:,3),selected_IDs);
grid3 = grid2(sample,:);
[~,rem] = intersect(grid(:,end),grid3(:,end));
grid(rem,:) = [];
% writematrix(grid3,'Selected.csv')
% writematrix(grid,'NotSelected.csv')

%% SAVE 
fname = "Charleston_CSRM_Subset_Grid";
% Save output for use in Storm Selection
% gridSM = [LAT	LON	ID]
gridSM = grid3;
save(input_final_path + fname + ".mat",'gridSM','-v7.3')

%% PLOT the CSRM locations
% Create geoaxes
ax=geoaxes('basemap','streets','TickLabelFormat','-dd','OuterPosition',[0 0 1 1]);
hold(ax,'on');
% Plot all savepoints
plt1=geoscatter(ax,grid(:,1),grid(:,2),150,[0 0 0],'.');
% Plot selected savepoints
plt2=geoscatter(ax,grid3(:,1),grid3(:,2),300,[204/255 0 0],'.');
% Plot Bounding Box
plt3=geoplot(ax,POLY(:,1),POLY(:,2),'Color',[255 255 0]./255,'LineWidth',2);
% Add Title and legend
title(ax, 'Charleston CSRM Location Selection' , 'FontSize' , 15)
legend([plt1 plt2 plt3],{'Not Selected' , 'Selected' , 'Bounding Box'})
% Add datatips
dt1=datatip(plt1);
plt1.DataTipTemplate.DataTipRows(1)=dataTipTextRow('ID',grid(:,3));
plt1.DataTipTemplate.DataTipRows(2)=dataTipTextRow('Longitude',grid(:,2));
plt1.DataTipTemplate.DataTipRows(3)=dataTipTextRow('Latitude',grid(:,1));
dt2=datatip(plt2);
plt2.DataTipTemplate.DataTipRows(1)=dataTipTextRow('ID',grid3(:,3));
plt2.DataTipTemplate.DataTipRows(2)=dataTipTextRow('Longitude',grid3(:,2));
plt2.DataTipTemplate.DataTipRows(3)=dataTipTextRow('Latitude',grid3(:,1));
% Geolimits
geolimits(ax,[32.3894 33.0308] , [-80.3484 -79.3813])
%If finished plotting, quit writing to axis
hold(ax,'off'); dt1.Visible='off'; dt2.Visible='off';
% Export Plots
% Removed deprecated 'compact'
savefig(gcf,fullfile(plot_folder,fname + ".fig"))
exportgraphics(gcf,fullfile(plot_folder,fname + ".png"),'Resolution',350)
% close all;

%% END