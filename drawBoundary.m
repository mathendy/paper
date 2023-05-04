function drawBoundary(BD,X)
%DRAWBOUNDARY 此处显示有关此函数的摘要
%   此处显示详细说明

for i=1:size(BD,2)
    hold on;
    plot(X(BD{i}(1:end),1),X(BD{i}(1:end),2),'r','Marker','o','MarkerSize',2);
end
end

