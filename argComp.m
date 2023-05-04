function arg = argComp(BD,X)
%ARGCOMP 此处显示有关此函数的摘要
%   此处显示详细说明
 arg=zeros(size(BD,2),1);
 for i=1:size(BD,2)
     bdX=X(BD{1,i},:);
     bdE=bdX(2:end,:)-bdX(1:end-1,:);
     a=bdE;
     b=[bdE(2:end,:);bdE(1,:)];
     argE=zeros(size(bdE,1),1);
     for j=1:size(bdE,1)
        argE(j,1)=atan2d(det([a(j,:);b(j,:)]),dot(a(j,:),b(j,:)));
     end
     arg(i)=sum(argE);
 end
end

