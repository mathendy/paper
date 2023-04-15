function TRSub = simpleSubdivision(TR,triNum,subNum)
%SIMPLESUBDIVISION 此处显示有关此函数的摘要
%   此处显示详细说明
T=TR.ConnectivityList;
P=TR.Points;
numP=size(TR.Points,1);
tri=T(triNum,:);
subP=P(T(triNum,:),:);
subtri=[1,2,3];
for i=1:subNum
    [subP,subtri]=LinearSubdivision(subP,subtri);
end
TRSub=triangulation(subtri,subP);
F = freeBoundary(TRSub);
tmpF=[F(:,1);F(1,1)];
subEdgeVer=find(tmpF==1|tmpF==2|tmpF==3);
subEdge=[tmpF(subEdgeVer(1):subEdgeVer(2))';tmpF(subEdgeVer(2):subEdgeVer(3))';tmpF(subEdgeVer(3):subEdgeVer(4))'];

subtri=subtri+numP-3;
subEdge=subEdge+numP-3;
for i=1:numel(subtri)
    if((subtri(i)==numP-2)||(subtri(i)==numP-1)||(subtri(i)==numP)) 
        subtri(i)=tri(subtri(i)-numP+3);
    end        
end
for i=1:3
    subEdge(i,1)=tri(subEdge(i,1)+3-numP);
    subEdge(i,end)=tri(subEdge(i,end)+3-numP);          
end

neiborsTri=[];
neiborsTriIdx=[];
for i=1:3
    ID = edgeAttachments(TR,subEdge(i,1),subEdge(i,end));
    idx=ID{:}(ID{:}~=triNum);
    if(~isempty(idx))
        neiborsTriIdx=[neiborsTriIdx,idx];
        neibor=T(idx,:);
        neiborsVer=neibor(neibor~=subEdge(i,1)&neibor~=subEdge(i,end));
        neiborsTri=[neiborsTri;[neiborsVer*ones(size(subEdge,2)-1,1),subEdge(i,2:end)',subEdge(i,1:end-1)']];
    end
end
P=[P;subP(4:end,:)]; 
T([neiborsTriIdx,triNum],:)=[];
T=[T;subtri;neiborsTri];
TRSub=triangulation(T,P);
end
