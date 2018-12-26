% this function performs modified SLIC supervoxel segmentation
% Input:
%   Y, each row is a voxel, each column is a feature
%   index, voxel indice of Y in an image A, A is with the same size of T
%   T is a 3D labeled template
%   CN, supervoxel number
% Output:
%   imag_SLIC is a 3D supervoxel segmentation and is with the same size of T
%   K is actual supervoxel number

function [img_SLIC,K]=ModifiedSLIC(Y,index,T,CN)

m=40;
[label]=SLIC(Y,index,m,CN,1,T);
img_SLIC=zeros(size(T));
img_SLIC(index)=label;
K=length(unique(label)); % actual cluster number