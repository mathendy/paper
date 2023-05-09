function drawBoundary(BD,X)
%DRAWBOUNDARY 此处显示有关此函数的摘要
%   此处显示详细说明
if(size(X,2)==2)
    for i=1:size(BD,2)
        hold on;
%         plot(X(BD{i}(1:(end/4+2)),1),X(BD{i}(1:(end/4+2)),2),'b','Marker','.','MarkerSize',2);
%         plot(X(BD{i}((end/4+1):(end/2+1)),1),X(BD{i}((end/4+1):(end/2+1)),2),'r','Marker','.','MarkerSize',2);
%         plot(X(BD{i}(floor(end/2+1):end),1),X(BD{i}(floor(end/2+1):end),2),'g','Marker','.','MarkerSize',2);
        plot(X(BD{i}(1:end),1),X(BD{i}(1:end),2),'b','Marker','.','MarkerSize',2);
    end
else
    if(size(X,3)==3)
        for i=1:size(BD,2)
            hold on;
            plot3(X(BD{i}(1:end),1),X(BD{i}(1:end),2),X(BD{i}(1:end),3),'b','Marker','.','MarkerSize',2);
        end
    end
end
end

