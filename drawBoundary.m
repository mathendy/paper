function drawBoundary(BD,X)
%DRAWBOUNDARY �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

for i=1:size(BD,2)
    hold on;
    plot(X(BD{i}(1:end),1),X(BD{i}(1:end),2),'r','Marker','o','MarkerSize',2);
end
end

