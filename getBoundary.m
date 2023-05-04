function BD = getBoundary(TR)
%GETBOUNDARY 此处显示有关此函数的摘要
%   此处显示详细说明
BdIdx=freeBoundary(TR);
i=BdIdx(1:end-1,2)~=BdIdx(2:end,1);
edIdx=find(i);
edIdx=[0;edIdx;size(BdIdx,1)];
BD=cell(1,size(edIdx,1)-1);
for i=1:size(edIdx,1)-1
    BD{i}=[BdIdx((edIdx(i)+1):edIdx(i+1),1);BdIdx(edIdx(i)+1,1)];
end  
end

