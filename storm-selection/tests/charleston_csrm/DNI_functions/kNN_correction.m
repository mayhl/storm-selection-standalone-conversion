function cor=kNN_correction(dis_neighbors,index_neighbors,surge,index_pt,i_wet,input,type,cor,k_totalmax,n_neighbors,parp,anc_surge)
% Enhance dry node predictions using calibrated kNN nearest neighob interpolation   
% For description of most inputs and outputc check kNN_DryNodeAdjustment_Main
%   parp: examine (1) or not (0) parallel implementation for efficiency
% V 1 2/23/2020 Alexandros Taflanidis (a.taflanidis@nd.edu)
% V 1.1 2/8/2021 Alexandros Taflanidis (a.taflanidis@nd.edu) corrected bug
%                related to total_k exceeding the number of available nodes for which 
%                distance has been estimated already 
% V 2 2/28/2021 Alexandros Taflanidis (a.taflanidis@nd.edu) Simplified version to 
%                use global indexing for the points to map-from. This
%                accomodates the updates to include ADCIRC connectivity 
%  V.2.1 3/6/2021 Alexandros Taflanidis (a.taflanidis@nd.edu)
%                 accomodate correction for small number of storms by doing
%                 only for loop (no parallel)
%  V.2.2 6/3/2021 Alexandros Taflanidis (a.taflanidis@nd.edu)
%                 accomodate correction for misclasified nodes that does
%                 notinterfere with nodes set to ancillary values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(cor.wet,'yes');wet_aux=0; else; wet_aux=1; end %dummy varriable to facilitate predictions for the storm wet points
if ~exist('parp','var'); parp=1; end; if isempty(parp);parp=1; end % set value of par=1 (examine parfor implementation if not defined)
if ~exist('anc_surge','var');anc_surge=[]; end % set value of ancilary surge if not defined
if isempty(n_neighbors); n_neighbors=k_totalmax*ones(size(dis_neighbors,1),1); end

%nx is number of variables to be optimized
if strcmp(type,'number'); nx=1;
elseif strcmp(type,'distance'); nx=2;
else; nx=4;
end
dist=dis_neighbors;%distance of closest neighbors
if nx<=2
    dist=1./dist;% inverse distance 
end
n_storms=size(surge,2); %number of storms 

total_k=cor.total_k; mc_thresh=cor.mc_thresh;total_kup=cor.total_kup; x=cor.x; x_up=cor.x_up; 
%iniitalize variables
Y_hat=zeros(size(index_pt,1),size(surge,2));

%loop over storms, check if parallel computations are more efficient by
%performing correction for first storms with or without paralell
%computing
if parp==1 %check if parallel will be examined 
    %start parfor loop
    delete(gcp('nocreate')); p=gcp;
    par=p.NumWorkers;%number of storms for the judginf efficiency of parallel computations, equal to number of workers
else 
    par=1;
end

%do computations with parallel
disp('correction performed with parallel implementation')
for j=1:ceil(n_storms/par) %decompose problem to subsets to par storms to improve memory use
    ai_wet=i_wet(:,(j-1)*par+1:min(j*par,n_storms)); %auxiliary variable to reduce memory use
    parfor i=(j-1)*par+1:min(j*par,n_storms)
        disp(i)
        %correct storm
        [Y,MC(i),MC_up(i),iter(i),index_cor]=storm_cor(surge(:,i),wet_aux*ai_wet((index_pt),i-(j-1)*par),...
            ai_wet(:,i-(j-1)*par),x,x_up,index_pt,total_k,dist,index_neighbors,...
            input(index_pt(ai_wet((index_pt),i-(j-1)*par)==0),1),ai_wet((index_pt),i-(j-1)*par)==0,total_kup,nx,k_totalmax,n_neighbors);
        %set misclasified nodes to elevation minus threshold
        %set misclasified nodes to elevation minus threshold
         if strcmp(cor.mc,'yes')
             if ~isempty(anc_surge) %check if there are ancillary surge values that should be avoided
                 aux_s=Y((index_cor))~=anc_surge; %identify misclasified values that do not correspond to ancillary surge
                 Y((index_cor(aux_s)))=input(index_pt(index_cor(aux_s)),1)-mc_thresh; %correct remaining misclassified values
             else
                 Y((index_cor))=input(index_pt(index_cor),1)-mc_thresh; %correct misclassified values
             end
         end
        %keep predictions and misclasified nodes in memory
        Y_hat(:,i)=Y; index_mc{i}=index_cor;
    end
end

%statistics and corrections per storm 
if n_storms>1 %if wet points exist for storm check predicitons for those
    for i=1:n_storms
        %validation statistics (correlation coefficient $ NRMSE) for wet nodes per storm
        aux_i=find(i_wet((index_pt),i)==1); 
        if length(aux_i)>1
            aux=corrcoef(surge(index_pt(aux_i),i),Y_hat(aux_i,i));SCC(i)=aux(2);
            SNRMSE(i)=sqrt(mean((surge(index_pt(aux_i),i)-Y_hat(aux_i,i)).^2,1))/(max(surge(index_pt(aux_i),i))-min(surge(index_pt(aux_i),i))); %normalized root mean squared error
            SAE(i)=mean(abs(surge(index_pt(aux_i),i)-Y_hat(aux_i,i))); %absolute error
            %update the values of the wet nodes with the correct ones
            Y_hat(aux_i,i)=surge(index_pt(aux_i),i);
        else
            SNRMSE(i)=0; SAE(i)=0;
        end
    end
else %set to NaN 
    SAE=NaN;
    SNRMSE=NaN; 
    SCC=NaN;
end
cor.AMC=mean(MC_up); %average missclasification 
cor.ACC=mean(SCC);% average correlation coefficient 
cor.ANRMSE=mean(SNRMSE); %average NRMSE
cor.AAE=mean(SAE); %absolute error

%pack remainign results in cor structure 
cor.Y_hat=Y_hat; cor.MC_or=MC; cor.MC=MC_up;cor.iter=iter; cor.SCC=SCC; cor.index_mc=index_mc; cor.SNRMSE=SNRMSE; cor.SAE=SAE; 

if parp==1; delete(gcp('nocreate')); end %close parfor loop 
end



function [Y_hat,MC,MC_up,iter,index_cor]=storm_cor(surge,i_wett,i_wetf,x,x_up,i_con,total_k,inv_dist,index_neighbors,input,i_wet,total_kup,nx,k_totalmax,n_neighbors)
%function to adjust surge for each storm 
%--------------------------------------- 
i_wetfo=i_wetf; %original wet origin points 
iter=1; %iteration counter
Y_hat=zeros(size(i_wett,1),1); %initialize variables 

while mean(i_wett)<1 % loop if there are still nodes to be corrected
si_wett=sum(i_wett); %total current number of wet nodes
    index_cor=find(i_wett==0); %indices of points to be corrected at the current iteration
    i_wetf_up=i_wetf; %variable with indices for updated wet points 
    %loop over the points
    for j=1:size(index_cor,1)
        [Y,ind]=cor_point(j,x,total_k,i_wetf,inv_dist,index_neighbors,index_cor,surge,nx,n_neighbors);
        %if point can be corrected, update characteristics
        if ind==1
            Y_hat(index_cor(j),1)=Y; % update predictions
            i_wett(index_cor(j))=1; %updated wet index
            %if point was initially dry update the surge value as well
            if i_wetf(i_con(index_cor(j)))==0
                surge(i_con(index_cor(j)))=Y_hat(index_cor(j),1); %updated surge
            end
            i_wetf_up(i_con(index_cor(j)))=1; %adjust indices for updated wet points
        end
    end
    i_wetf=i_wetf_up; %new set of wet points 
    %if sum(i_wett)-si_wett<=1; k=k-1;end %if no corrections have been performed at current iteration reduce value of k by 1
    %if k==0; k=1; total_k=total_k+1; end %if end up with value of k=0 (no point can be corrected) update it to k=1 and increase total_k to get a well-posed problem
    if sum(i_wett)-si_wett<=1 %if no corrections have been performed at current iteration
        total_k=min(k_totalmax,total_k+1); %increase value of neighbors to search (to find k wet) by 1
        if total_k==k_totalmax %if available neighbors exhausted modify k
            x(1)=x(1)-1; %reduce k by 1
            if x(1)==0; x(1)=1; total_k=k_totalmax+1; n_neighbors=n_neighbors+1; end %if reached k=0 then add the anciliary nodes that are always wet
        end
    end
    iter=iter+1; %update iteration counter
end

%calculate misclassification
aux=Y_hat(i_wet,1)-input; MC=mean(aux>0);
%identify misclasified nodes
aa=find(i_wet); index_cor=aa(aux>0);

%pefrorm another round of correction using only the original data and
%potentially updated values for the k NN
if isempty(x_up) %skip second round if x_up is not defined
    MC_up=MC;
else
    i_wetf=i_wetfo; % original wet points
    %perform another correction for these nodes
    i_wett(index_cor)=0; %set nodes to be adjusted as dry
    while mean(i_wett)<1 % loop if there are still nodes to be corrected
        si_wett=sum(i_wett); %total current number of wet nodes
        index_cor=find(i_wett==0); %indices of points to be corrected at the current iteration
        i_wetf_up=i_wetf; %variable with indices for updated wet points
        for j=1:size(index_cor,1)
            [Y,ind]=cor_point(j,x_up,total_kup,i_wetf,inv_dist,index_neighbors,index_cor,surge,length(x_up));
            %if point can be corrected, update characteristics
            if ind==1
                Y_hat(index_cor(j),1)=Y; %predictions
                %if point was initially dry update the surge value as well
                if i_wetf(i_con(index_cor(j)))==0
                    surge(i_con(index_cor(j)))=Y_hat(index_cor(j),1);
                end
                i_wetf_up(i_con(index_cor(j)))=1;
            end
        end
        i_wetf=i_wetf_up; %new set of wet points
        if sum(i_wett)-si_wett<=1; i_wett(index_cor)=1;end %if no corrections have been performed at current iteration quit
    end
    
    %calculate updated misclassification
    aux=Y_hat(i_wet,1)-input; MC_up=mean(aux>0);
    %identify misclasifies nodes
    aa=find(i_wet); index_cor=aa(aux>0);
end
end



function [Y,ind]=cor_point(j,x,total_k,i_wetf,dist,index_neighbors,index_cor,surge,nx,n_neighbors)
%function to adjust surge for each point 
%---------------------------------------
%how many of the k-closest neighbors are wet
i_wetk=find(i_wetf(index_neighbors(index_cor(j),1:min(total_k,n_neighbors(index_cor(j)))))==1);
    k=round(x(1));
%check if sufficient number of neighbors are wet and if yes perform correction, else skip point
if length(i_wetk)>=k
    if nx==1 %kNN optimization
        weights=dist(index_cor(j),i_wetk(1:k)); weights=weights./sum(weights,2);
    elseif nx==2    
         weights=dist(index_cor(j),i_wetk(1:k)); weights(weights<1/x(2))=0;
         weights=(weights+eps)./sum(weights+eps,2);
    else
        weights=exp(-(dist(index_cor(j),i_wetk(1:k))/x(3)).^x(4)); weights(dist(index_cor(j),i_wetk(1:k))>x(2))=0;
        weights=(weights+eps)./sum(weights+eps,2);
    end
    aux_index=index_neighbors(index_cor(j),i_wetk(1:k)); % indices of the k closest neigbors
    Y=weights*surge(aux_index); %predictions
    ind=1; %indicator function for whether point was corrected
else
    Y=0; ind=0; %if point skipped set ind==0
end
end