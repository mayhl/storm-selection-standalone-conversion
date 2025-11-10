%% XC Selection by AEF
% Madison Yawn
% 21 Sept 2023

clear; clc; close all;
% 
% %% Load Data
load('R:\Users\ERS\arc_from74_20240311\Task_FY21_SACS\Phase2_SST_Surge\XC_SLC0\Step2_DNC\SACS_Ph2_XC_SLC_0_Tides_1_surge_SPs_cor.mat')
load('R:\Work\CHS_PCHA\arc_from74\PCHA_v2021.03_copy\Atlantic_20200428\8_CHS_Studies\USACE\SACS_NCSEFL\7_HazardCurves_rta\Surge\XC_SLC0_Tide0\out\StormSim_SST_output.mat')
% load('H:\Desktop\SACS\SA\TC_XC_Comparison_Plots\AEF_Data\SWL\StormSim_SST_output_SLC0.mat')

%% User Inputs
svpts = 11902; % svpts of interest 
num_storms = 20;
screening_storms = 5;

tbl_aef = 1./[0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000];

IDs = nan(length(SST_output),1);
for j = 1:length(SST_output)
    IDs(j) = str2double(SST_output(j).staID);
end

v = nan(length(tbl_aef),length(svpts));
for i = 1:length(svpts)
    ix = find(svpts(i,1)==IDs);
    if sum(ix)>0
        v(:,i) = SST_output(ix).HC_tbl(1,1:16)';
    end
end

[resp, ix] = sort(Resp(:,svpts),'descend'); % WLs at specific svpt
vq = interp1(v(4:end,1),log(tbl_aef(4:end))',resp);
eval = [exp(vq),resp,ix];

%% Select Storms
% aefs = (1/100):0.025:(1/2);
aefs = min(eval(:,1)):(max(eval(:,1)-min(eval(:,1)))/(num_storms-1)):max(eval(:,1));
XC_RSS = nan(length(aefs),4);
% col 1: aef of interest
% col 2: interpolated AEF from HC
% col 3: response magnitude
% col 4: storm index

for i = 1:length(aefs)
    if i==or(1,length(aefs))
        ix = find(aefs(i)==eval(:,1));
        XC_RSS(i,:) = [aefs(i),eval(ix,:)];
    else
        diff = aefs(i)-eval(:,1); % create difference in aef and AEFs from eval
        diff(diff<0)=[];
        ix = find(diff==min(diff)); % find which one has the smallest difference
        XC_RSS(i,:) = [aefs(i),eval(ix,:)]; % select storm and then remove from eval list
    end
    eval(ix,:) = [];
end

%% Screening Storms
rng('default')
XC_samp = randsample(1:1:num_storms,screening_storms);
XC_screening = XC_RSS(XC_samp,:);

% % aefs = (1/100):0.025:(1/2);
% aefs = 1./[1, 10, 25, 50, 100];
% XC_screening = nan(length(aefs),4);
% % col 1: aef of interest
% % col 2: interpolated AEF from HC
% % col 3: response magnitude
% % col 4: storm ID
% XC_RSS_temp = XC_RSS;
% for i = 1:length(aefs)
%     diff = abs(aefs(i)-XC_RSS_temp(:,1)); % create difference in aef and AEFs from eval
%     ix = find(diff==min(diff)); % find which one has the smallest difference
%     XC_screening(i,:) = XC_RSS_temp(ix,:); % select storm and then remove from eval list
%     XC_RSS_temp(ix,:) = [];
% end

cstorm_IDs = [1:1:73]';
cstorm_IDs([25,29,45],:) = [];
save('Charleston_XC_Selected_Storms.mat','XC_RSS','XC_screening','cstorm_IDs')

% for .csv table
storm_names = readtable('R:\Archive\6_CSTORM_Modeling\Peaks\CHS-SA\sacs_sa_et_his_list_of_storms.xlsx');
cstorm_rss = storm_names(cstorm_IDs(XC_RSS(:,4)),1:2);
cstorm_screening = storm_names(cstorm_IDs(XC_screening(:,4)),1:2);

%% Check Results
XTickMarks=[1e-2 1e-1 1 10];
% X tick mark labels
XTickMarkLabels={'10^{-2}','10^{-1}','10^{0}','10^{1}'};
% X-axis limits
XAxisLimits=[1e-2 10];
% Calculate Ymax
Ymax=4;
% Y tick marks
if Ymax<=10
    YTickMarks=1:1:Ymax;
else
    if rem(Ymax,2)==0
        YTickMarks=2:2:Ymax;
    else
        YTickMarks=1:2:Ymax;
    end
end
% Y-axis limits
YAxisLimits=[0 Ymax];
 figure('Color',[1 1 1],'visible','on');
        axes('XScale','log','XGrid','on','XMinorTick','on','YGrid','on','YMinorTick','on','FontSize',16,...
            'YTick',YTickMarks,'XTick',XTickMarks,...
            'XTickLabel',XTickMarkLabels);
        xlim(XAxisLimits); ylim(YAxisLimits);
        set(gca,'XDir','reverse')
        hold on
plot(tbl_aef,v','LineWidth',1)
hold on
scatter(XC_RSS(:,1),XC_RSS(:,3),750,'.');
hold on
scatter(XC_screening(:,1),XC_screening(:,3),500,'.');
legend('SST HC','XC RSS','Screening XCs')
