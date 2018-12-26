% This function is used to perform FVGWAS
% Input:
%   Y are imaging data with size n*V, where n is subject number and V is voxel
%   number
%   SNP are genetic data with size n*C, where C is loci number
%   X are clinical covariates with size n*p, where p is covariate number
%   index is a V*1 vector, which is the index of voxels in an image
%   sizeimg is image size
%   N0 is number of top selected SNPs
%   Flag_S is a flag that indicates if SNP processing is needed
% Output:
%   pv are corrected p-values of significant voxel-locus pairs
%   rawpvalue are raw p-values of significant voxel-locus pairs
%   Cluster is a cell, including information of significant cluster-locus pairs

function [pv,Cluster,rawpvalue]=FVGWAS(Y,SNP,X,index,sizeimg,N0,Flag_S)

lamda=0.001;%%%%%
PX=X/(X'*X+lamda*eye(size(X,2)))*X';
n=size(X,1);%%%%%number of individual
p=size(X,2);%%%%feature dimention of X

%%%%%%%%%%%%%%step 1: Get points and voxels
BatchSize=5000;
deta=zeros(size(Y,2),1);
if length(index)>BatchSize
    for kk=0:fix(length(index)/BatchSize)
        if kk<fix(length(index)/BatchSize)
            temp=Y(:,BatchSize*kk+1:BatchSize*(kk+1))'*(eye(size(PX,1))-PX)*Y(:,BatchSize*kk+1:BatchSize*(kk+1))/(n-p);
            deta(BatchSize*kk+1:BatchSize*(kk+1))=diag(temp);
        else
            if BatchSize*kk+1>length(index)
                break;
            end
            temp=Y(:,BatchSize*kk+1:end)'*(eye(size(PX,1))-PX)*Y(:,BatchSize*kk+1:end)/(n-p);
            deta(BatchSize*kk+1:end)=diag(temp);
        end
    end
end
deta(deta(:)~=0)=deta(deta(:)~=0).^(-1);


%%%%%%%%%step 2: exclude SNPs with MAF values smaller than 0.05. Remaining missing genotype variables were
%%%%%imputed as the modal value
if Flag_S
    num=zeros(3,size(SNP,2));
    for i=1:3
        bw=zeros(size(SNP));
        bw(SNP(:)==i-1)=1;
        num(i,:)=sum(bw);
    end
    [~,maxnum]=max(num);
    maxnum=maxnum-1;
    [r,c]=find(SNP<0);
    for i=1:length(r)
        SNP(r(i),c(i))=maxnum(c(i));
    end
    %%%%%%%%%%%%%%
    minMAF=0.05;%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MAF=sum(SNP)/(2*n);
    temp=find(MAF>0.5);
    MAF(temp)=1-MAF(temp);
    SNP_index=MAF(:)<=minMAF;
    SNP(:,SNP_index)=[];
end

% % %%%%%step 3: use the chi-squared approximation to the observed W(c) to calculate the p-values
tic
[pp, ~]=globalWald(SNP,Y,PX,deta);
fprintf('GSIS: %f\n',toc)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G=1000;

%%%%%%step 4: generate bootstrap samples
tic
[Wgmax,Ngmax]=Wholewild_MaxWgNg(N0,X,PX,Y,deta,SNP,G,index,sizeimg);
fprintf('Generate bootstrap samples: %f\n',toc)


[~,indx]=sort(pp,'descend');
SNPtemp=SNP(:,indx(1:N0));

%%%%%%step 5: detect significant locus-voxel and locus-cluster pairs
tic
[pv,Cluster,rawpvalue]=Genvoxclusnp_New(N0,PX,X,Y,SNPtemp,Wgmax,Ngmax,deta,index,sizeimg);

fprintf('Generate bootstrap samples: %f\n',toc)