function L = cotLaplace(x,t)
%COTLAPLACE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
l=meshFaceEdgeLen2s(x,t)
[angles,nbroktris]=meshAnglesFromFaceEdgeLen2(l);
L=cotLaplaceFromFaceAngles(angles,t,size(x,1));
end

