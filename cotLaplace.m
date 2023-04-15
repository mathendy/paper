function L = cotLaplace(x,t)
%COTLAPLACE 此处显示有关此函数的摘要
%   此处显示详细说明
l=meshFaceEdgeLen2s(x,t)
[angles,nbroktris]=meshAnglesFromFaceEdgeLen2(l);
L=cotLaplaceFromFaceAngles(angles,t,size(x,1));
end

