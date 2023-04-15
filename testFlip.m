function [flipInd,k] = testFlip(testMeshX,testMeshT,testDeformMeshX)
%TESTFLIP 此处显示有关此函数的摘要
%   此处显示详细说明
% Z1=testMeshX(testMeshT(:,2),:)-testMeshX(testMeshT(:,1),:);Z2=testMeshX(testMeshT(:,3),:)-testMeshX(testMeshT(:,1),:);
% Z1=[Z1,zeros(size(testMeshT,1),1)];Z2=[Z2,zeros(size(testMeshT,1),1)];
% Z=cross(Z1,Z2);
% sigZ=sign(Z(:,3));
% 
% Z1_deform=testDeformMeshX(testMeshT(:,2),:)-testDeformMeshX(testMeshT(:,1),:);Z2_deform=testDeformMeshX(testMeshT(:,3),:)-testDeformMeshX(testMeshT(:,1),:);
% Z1_deform=[Z1_deform,zeros(size(testMeshT,1),1)];Z2_deform=[Z2_deform,zeros(size(testMeshT,1),1)];
% Z_deform=cross(Z1_deform,Z2_deform);
% sigZ_deform=sign(Z_deform(:,3));
% 
% k=sum((sigZ-sigZ_deform)~=0);
% flipInd=find(sigZ-sigZ_deform);
flipInd=zeros(size(testMeshT,1),1);
AT=zeros(size(testMeshT,1),1);
x1=testMeshX(testMeshT(:,2),:)-testMeshX(testMeshT(:,1),:);
x2=testMeshX(testMeshT(:,3),:)-testMeshX(testMeshT(:,1),:);
for i=1:size(x1,1)
    AT(i)=abs(det([x1(i,:);x2(i,:)]));
end 
Jt=zeros(2*size(testMeshT,1),2);
for i=1:size(x1,1)
    tri=testMeshX(testMeshT(i,:),:)';
    J=1/(2*AT(i))*[circshift(tri(2,:),[-1,0])-circshift(tri(2,:),[1,0]);circshift(tri(1,:),[1,0])-circshift(tri(1,:),[-1,0])]*testDeformMeshX(testMeshT(i,:),:);
    Jt(2*i-1:2*i,:)=J;
    if det(J)<0
        flipInd(i)=1;
    end
end
k=sum(flipInd);
end

