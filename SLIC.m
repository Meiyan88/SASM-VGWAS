function [label,time]=SLIC(img_gray,ind_gray,m,cN,flag,img_seg)

% Simple linear iterative clustering (SLIC)
% Input:
%     img_gray, each row is a voxel, each column is a feature
%     ind_gray, voxel indice of img_gray
%     img_seg, brain mask of 90 ROIs
%     m, the tuning parameter
%     cN, cluster number
%     flag=1 uses the anatomical information
% Output:
%     label, the clustering result
%     time, the elapsed time to do clustering by SLIC

tic;

siz=size(img_seg);
a=siz(1);
b=siz(2);
c=siz(3);
num_trial=size(img_gray,2);
img_rs=zeros(prod(siz),num_trial);
img_rs(ind_gray,:)=img_gray;

% unified image to store the image and the coordinates for each
% voxel in the gray matter mask
[A,B,C]=ind2sub(siz,ind_gray);
img_uni=[A,B,C,img_rs(ind_gray,:)];

% initialize the cluster centers as the sphere centers
num_gray=length(ind_gray);
S0=(num_gray/(cN*sqrt(2)*4))^(1/3); % the radius of the spheres
na=a/(S0*sqrt(3));
nb=b/(S0*2);
nc=c/(S0*sqrt(6)*2/3);
nEst=floor(na)*floor(nb)*floor(nc); % estimated cluster number in 3D space
center=zeros(nEst,num_trial+3);
d=1; % a tolerance to search for wider space
nK=0; % to save the exact cluster number
for i=1-d:na+d
    for j=1-d:nb+d
        for k=1-d:nc+d
            % location of the spatial center
            A=S0*(i-1)*sqrt(3)+S0+S0/sqrt(3)*mod(k,3);
            B=S0*(j-mod(i,2)/2)*2+S0*mod(k,3);
            C=S0*(k-1)*sqrt(6)*2/3+S0;
            
            ix=[A,B,C];
            ix=floor(ix);
            
            % box constraint
            if min(ix)>=1 && all(ix<=siz)
                index=sub2ind(siz,ix(1),ix(2),ix(3)); % single index for the subscript
                
                % convert index in 3D space to that in graymatter mask
                [~,index]=ismember(index,ind_gray);
                
                % if in graymatter mask, store the cluster center in double
                if index~=0
                    nK=nK+1;
                    center(nK,:)=[A,B,C,img_uni(index,4:end)];
                end
            end
        end
    end
end
center=center(1:nK,:); % cut the zeros

% create searchlight mask to search in the 2r*2r*2r neighborhood
S=(num_gray/cN)^(1/3); % sampling interval
[riSL,nvSL]=parc_searchlight(1.5*S,'cube');

% initialization for SLIC
label=-1*ones(a,b,c); % set the labels of all voxels to be -1
distance=Inf*ones(a,b,c);

rsd=1;
iter=0;

while rsd>1e-4
    iter=iter+1;
    center_pre=center;
    
    % assignment
    for iK=1:nK
        % searchlight constraint
        ix=riSL+repmat(center(iK,1:3),[nvSL,1]);
        ix=floor(ix); % coordinates should be integers
        
        % box constraint
        tmp=all(ix>=ones(nvSL,3),2) & all(ix<=repmat(siz,nvSL,1),2);
        ix=ix(tmp,:);
        
        % mask constraint
        % convert index in 3D space to that in graymatter mask
        index=sub2ind(siz,ix(:,1),ix(:,2),ix(:,3)); % single index
        [tmp,index]=ismember(index,ind_gray);
        index=index(tmp);    
        ix=ix(tmp,:);
        
        nV=size(ix,1); % number of remaining voxels in this searchlight
        for iV=1:nV
            % for each voxel in the searchlight
            i=ix(iV,1);
            j=ix(iV,2);
            k=ix(iV,3);
            
            tmp=img_uni(index(iV),:)-center(iK,:);
            ds=norm(tmp(1:3)); % spatial distance
            dc=norm(tmp(4:end)); % color distance
            if flag
                % if voxel v and cluster center k have the same anatomical
                % label, distance is 0, otherwise, distance is positive
                % infinite
                temp1=img_uni(index(iV),1:3);
                temp1=sub2ind(siz,temp1(1),temp1(2),temp1(3));
                temp2=floor(center(iK,1:3));
                temp2=sub2ind(siz,temp2(1),temp2(2),temp2(3));
                if img_seg(temp1)==img_seg(temp2)
                    da=0;
                else
                    da=10^6;
                end
                du=sqrt((dc/m)^2+(ds/S)^2+da); % unified distance
            else
                du=sqrt((dc/m)^2+(ds/S)^2); % unified distance
            end
            
            if du<distance(i,j,k)
                distance(i,j,k)=du;
                label(i,j,k)=iK;
            end
        end
    end
    
    % update
    % compute the new cluster centers
    for iK=1:nK
        index=find(label==iK);
        if ~isempty(index)
            % convert index in 3D space to that in graymatter mask
            [~,index]=ismember(index,ind_gray);

            center(iK,:)=mean(img_uni(index,:));
        else
            center(iK,:)=zeros(1,num_trial+3);
        end
    end
    
    % compute residual error between current cluster center and previous
    % cluster center
    % absolute error since it's not appropriate to use relative error
    tmp=center-center_pre;
    rsd=norm(tmp(:));
    
    fprintf('Iteration number: %d, residual error: %0.4f. \n',iter,rsd);
end
label=label(ind_gray); % vector form
time=toc/3600;
fprintf('Time to do clustering by SLIC: %0.2f hours. \n\n',time);