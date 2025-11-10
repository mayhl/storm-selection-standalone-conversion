function [opt,val]=kNN_calibration(opt,dis_neighborsWet,index_neighborsWet,index_ptWet,surge,elev)
% Calibraiton of kNN nearest neighob interpolation   
% For description of inputs and outputc check kNN_DryNodeAdjustment_Main
% V 1 2/23/2020 Alexandros Taflanidis (a.taflanidis@nd.edu)
% V 1.1 2/27/2021 Alexandros Taflanidis (a.taflanidis@nd.edu) simplified 
%                 implementation by translating nodes to map-from to global indexing  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kmax=opt.ub(1); %max k that could be examined
n_storms=size(surge,2); %number of storms in database to use 

%nx is number of variables to be optimized
if strcmp(opt.type,'number'); nx=1;
elseif strcmp(opt.type,'distance'); nx=2;
else; nx=4;
end

% keep only points that satisfy distance criteria for calibration
if isempty(opt.dist_kept)
    aux_kept=1:size(dis_neighborsWet,1);
else
    aux_kept=find(dis_neighborsWet(:,opt.dist_kept(1))<opt.dist_kept(2));
end

Y=surge(index_ptWet(aux_kept),:); % actual response of points to predict to
dist=dis_neighborsWet(aux_kept,1:kmax);%distance of closest neighbors
if nx<=2
    dist=1./dist;% inverse distance 
end
%aux_index=index_pfWet(index_neighborsWet(aux_kept,1:kmax));% indices of the k closest neigbors for each point
aux_index=index_neighborsWet(aux_kept,1:kmax);% indices of the k closest neigbors for each point

%perform optimization by appropriate case
if strcmp(opt.type,'number') %if simple kNN do exhaustive search 
    x_cand=opt.lb:1:opt.ub; %search over all k values
    %perform search
    f=zeros(length(x_cand),1);
    for i=1:length(x_cand)
        f(i)=ObjfkNN_calibration(x_cand(i),dist,aux_index,Y,surge,n_storms,nx,opt.statistics,opt.criterion);
    end
    %keep in meory results
    opt.f=f; opt.x_cand=x_cand; 
    %perform optimization
    [~,j]=min(f);opt.x_opt=x_cand(j); opt.f_opt=f(j); 
elseif strcmp(opt.method,'random')
    if isempty(opt.x_cand)
        x_cand=repmat((opt.ub(1:nx)-opt.lb(1:nx)),opt.func,1).*lhsdesign(opt.func,length(opt.lb(1:nx)))+repmat(opt.lb(1:nx),opt.func,1);    %create ransom samples
        x_cand(:,1)=round(x_cand(:,1)); %set k to be intefer
    else
        x_cand=opt.x_cand;
    end
    %perform random search  
    f=zeros(size(x_cand,1),1);
    parfor i=1:size(x_cand,1)
        f(i)=ObjfkNN_calibration(x_cand(i,(1:nx)),dist,aux_index,Y,surge,n_storms,nx,opt.statistics,opt.criterion);
    end
    %keep in meory results
    opt.f=f; opt.x_cand=x_cand; 
    %perform optimization
    [~,j]=min(f);opt.x_opt=x_cand(j,:); opt.f_opt=f(j); 
else
    gaopts = gaoptimset('TimeLimit',3600,'Generations',opt.func/25/nx,'InitialPopulation',[],'PopulationSize',25*nx,'PlotFcns',@gaplotbestfun,'UseParallel',true);
    [opt.x_opt,opt.f_opt]=ga(@(x)ObjfkNN_calibration(x,dist,aux_index,Y,surge,n_storms,nx,opt.statistics,opt.criterion),...
    nx,[],[],[],[],opt.lb(1:nx),opt.ub(1:nx),[],1);
end

% now calculate resposne for optimized value
Y_hat=ObjfkNN_calibration(opt.x_opt,dist,aux_index,Y,surge,n_storms,nx,opt.statistics,opt.criterion,1);

%keep responses in memory
val.Y=Y; %actual response
val.Y_hat=Y_hat; %predicted response
% calculate validation statistics per point
dy=Y-Y_hat; %difference of prediction to real response
val.LAE=mean(abs(dy),2); %absolute error
val.LR2=1-mean(dy.^2,2)./(var(Y,1,2)+0.0001); % R2
val.LNRMSE=sqrt(mean(dy.^2,2))./(max(Y,[],2)-min(Y,[],2)); %normalized root mean squared error
val.LCC=mean((Y-repmat(mean(Y,2),1,n_storms)).*(Y_hat-repmat(mean(Y_hat,2),1,n_storms)),2)./(std(Y,1,2).*std(Y_hat,1,2)+0.0001); %correlation coefficient
val.LMC=mean(Y_hat<repmat(elev,1,n_storms),2); %misclasification %

% calculate validation statistics per storm
val.SAE=mean(abs(dy),1); %absolute error
val.SR2=1-mean(dy.^2,1)./(var(Y,1,1)+0.0001); % R2
val.SNRMSE=sqrt(mean(dy.^2,1))./(max(Y,[],1)-min(Y,[],1)); %normalized root mean squared error
val.SCC=mean((Y-repmat(mean(Y,1),size(Y,1),1)).*(Y_hat-repmat(mean(Y_hat,1),size(Y,1),1)),1)./(std(Y,1,1).*std(Y_hat,1,1)+0.0001);%correlation coefficient
val.SMC=mean(Y_hat<repmat(elev,1,n_storms),1); %misclasification %

% calculate average statistics across all point
val.AAE=mean(val.SAE); %absolute error
val.AR2=mean(val.SR2); % R2
val.ANRMSE=mean(val.SNRMSE); %normalized root mean squared error
val.ACC=mean(val.SCC); %correlation coefficient
val.AMC=mean(val.SMC); %misclasification 
end

function f=ObjfkNN_calibration(x,dist,aux_index,Y,surge,n_storms,nx,stat,cal,aux_opt)

k=round(x(1));
if nx==1 %kNN optimization 
    weights=dist(:,1:k)./sum(dist(:,1:k),2); % weights for k closest neighbors bases on inverse distance
elseif nx==2
    weights=dist(:,1:k); weights(dist(:,1:k)<1/x(2))=0;
    weights=(weights+eps)./sum(weights+eps,2);
else
    weights=exp(-(dist(:,1:k)/x(3)).^x(4)); weights(dist(:,1:k)>x(2))=0;
    weights=(weights+eps)./sum(weights+eps,2);
end
    
% sum contributions from each neighbor for all points and storms, use a
% for loop to accomodate the triple indexing (across storms, points,
% neighbors). Loop with respect to the neighbors since that is the
% lower dimension
Y_hat=zeros(size(Y)); %initialize variable
for i=1:k
    Y_hat=Y_hat+surge(aux_index(:,i),:).*repmat(weights(:,i),1,n_storms); %predictions for each point
end
% Y_hat is the predicted responce of each point to predict to
if nargin<10 %check if optimization of estimation of response warranted
    if strcmp(stat,'storms'); aux=1;else;aux=2;end %across what dimension to average
    if strcmp(cal,'NRMSE') %what criterion to use
        f=mean(sqrt(mean((Y-Y_hat).^2,aux))./(max(Y,[],aux)-min(Y,[],aux))); %average (over all points) normalized root mean squared error
    else
        f=mean(sqrt(mean(abs((Y-Y_hat).^2),aux)));
    end
else
    f=Y_hat;
end

end
