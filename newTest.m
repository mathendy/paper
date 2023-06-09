%% start
working_dataset='test';
fC2R = @(x) [real(x) imag(x)];
fR2C = @(x) complex(x(:,1), x(:,2));
datadir = fullfile(cd, 'data', working_dataset, '\');
imgfilepath  = fullfile(datadir, 'image.png');

%% test for a real floor plan
outerBoundary=[515,152;1209,152;1209,442;1143,442;1175,567;1175,873;1195,873;1195,1333;824,1333;824,1581;739,1581;739,1767;708,1767;708,1870;546,1870;546,1811;476,1811;476,1255;420,1255;...
    420,1061;351,1061;351,581;313,581;313,333;515,333;515,152]*2/1870-[1,1];
holes1=[916,442;1049,442;1049,583;1083,583;1083,683;916,683;916,442]*2/1870-[1,1];
holes2=[512,778;640,778;640,949;512,949;512,778]*2/1870-[1,1];
holes3=[998,1109;1093,1109;1093,1176;998,1176;998,1109]*2/1870-[1,1];
holes5=[-0.35,0.9;-0.25,0.9;-0.25,0.7;-0.35,0.7;-0.35,0.9];
plot(outerBoundary(:,1),outerBoundary(:,2),holes1(:,1),holes1(:,2),holes2(:,1),holes2(:,2),holes3(:,1),holes3(:,2),holes5(:,1),holes5(:,2));
%plot(outerBoundary(:,1),outerBoundary(:,2),holes1(:,1),holes1(:,2),holes2(:,1),holes2(:,2),holes3(:,1),holes3(:,2));
holeCenter1=(982+562i)*2/1870-(1+1i);
holeCenter2=(576+863i)*2/1870-(1+1i);
holeCenter3=(1045+1172i)*2/1870-(1+1i);
holeCenter4=-0.05+0.7i;
holeCenter5=-0.3+0.8i;

axis equal
%[X, T] = cdt([{outerBoundary},{holes1,holes2,holes3,holes4}],[], 50000, false);
[X, T] = cdt([{outerBoundary},{holes1,holes2,holes3,holes5}],[], 50000, false);

%% test for 3 holes ring mesh
arg1=(0:2:360)'/180*pi;
holeCenter1=3i/5;
holeCenter2=(1.0-0.6i)/2;
holeCenter3=(-1.0-0.6i)/2;
% holeCenter1=-1i/2;
% holeCenter2=1i/2;
% holeCenter3=1i/2;


outerBoundary=[cos(arg1),sin(arg1)];
r1=0.15;
r2=0.15;
r3=0.15;
arg2=(0:6:360)'/180*pi;
arg3=(0:6:360)'/180*pi;
arg4=(0:6:360)'/180*pi;
holes1=[r1*cos(arg2),r1*sin(arg2)]/2+fC2R(holeCenter1);
holes2=[r2*cos(arg3),r2*sin(arg3)]/2+fC2R(holeCenter2);
holes3=[r3*cos(arg4),r3*sin(arg4)]/2+fC2R(holeCenter3);

holeCenter4=0;
r5=0.25;
arg5=(0:6:360)'/180*pi;
holes4=[r5*cos(arg5),r5*sin(arg5)]/2+fC2R(holeCenter4);
% plot(outerBoundary(:,1),outerBoundary(:,2),holes1(:,1),holes1(:,2),holes2(:,1),holes2(:,2),holes3(:,1),holes3(:,2));

%[X, T] = cdt([{outerBoundary},{holes1,holes2}],[], 10000, false);
%[X, T] = cdt([{outerBoundary},{holes1,holes2,holes3}],[], 10000, false);
[X, T] = cdt([{outerBoundary},{holes1,holes2,holes3,holes4}],[], 200000, false);


%%
n_outerBoundary=[[(-1:0.1:1)';ones(19,1);(1:-0.1:-1)';-ones(20,1)],[ones(20,1);(1:-0.1:-1)';-ones(19,1);(-1:0.1:1)']];
n_hole=0.3*n_outerBoundary;

[X, T] = cdt([{n_outerBoundary},{n_hole}],[], 100000, false);
holeCenter1=0;
holeCenter2=0;
holeCenter3=0;
holeCenter4=0;


%%
n_outerBoundary=[[(-1:0.1:1)';ones(19,1);(1:-0.1:-1)';-ones(20,1)],[ones(20,1);(1:-0.1:-1)';-ones(19,1);(-1:0.1:1)']];
holeCenter1=-0.6+0.6i;
holeCenter2=0.6+0.6i;
holeCenter3=0;
holeCenter4=0.6-0.6i;
holeCenter5=-0.6-0.6i;

n_hole1=0.15*n_outerBoundary+fC2R(holeCenter1);
n_hole2=0.15*n_outerBoundary+fC2R(holeCenter2);
n_hole3=0.15*n_outerBoundary+fC2R(holeCenter3);
n_hole4=0.15*n_outerBoundary+fC2R(holeCenter4);
n_hole5=0.15*n_outerBoundary+fC2R(holeCenter5)

[X, T] = cdt([{n_outerBoundary},{n_hole1,n_hole2,n_hole3,n_hole4,n_hole5}],[], 50000, false);


%%
figure(1);

CX=fR2C(X);
TR=triangulation(T,X);
h=triplot(TR,'Color','#D3D3D3');
BD=getBoundary(TR);
nBD=size(BD,2);
hold on;
drawBoundary(BD,X);
% s=3-3.25;ed=3-3.35;
% plot([s;s;ed;ed;s],[-1;1;1;-1;-1],'r');
axis equal;
axis on;

argHole=argComp(BD,X);




%% integration
%2 holes
% a1=-(holeCenter1+holeCenter2)/2;
% a2=holeCenter1*holeCenter2;
% fDeformInt=@(x) 1/3*x.^3+a1*x.^2+a2*x;

%3 holes
% a1=-(holeCenter1+holeCenter2+holeCenter3)/3;
% a2=(holeCenter1*holeCenter2+holeCenter1*holeCenter3+holeCenter2*holeCenter3)/2;
% a3=-holeCenter1*holeCenter2*holeCenter3;
% fDeformInt=@(x) 1/4*x.^4+a1*x.^3+a2*x.^2+a3*x;

%4 holes
% a1=-(holeCenter1+holeCenter2+holeCenter3+holeCenter4)/4;
% a2=(holeCenter1*holeCenter2+holeCenter1*holeCenter3+holeCenter2*holeCenter3+holeCenter1*holeCenter4+holeCenter4*holeCenter3+holeCenter2*holeCenter4)/3;
% a3=-(holeCenter1*holeCenter2*holeCenter3+holeCenter1*holeCenter2*holeCenter4+holeCenter1*holeCenter4*holeCenter3+holeCenter4*holeCenter2*holeCenter3)/2;
% a4=holeCenter1*holeCenter2*holeCenter3*holeCenter4;
% fDeformInt=@(x) 1/5*x.^5+a1*x.^4+a2*x.^3+a3*x.^2+a4*x;

%5 holes
a1=-(holeCenter1+holeCenter2+holeCenter3+holeCenter4+holeCenter5)/5;
a2=(holeCenter1*holeCenter2+holeCenter1*holeCenter3+holeCenter2*holeCenter3+holeCenter1*holeCenter4+holeCenter4*holeCenter3+holeCenter2*holeCenter4...
    +(holeCenter1+holeCenter2+holeCenter3+holeCenter4)*holeCenter5)/4;
a3=-(holeCenter1*holeCenter2*holeCenter3+holeCenter1*holeCenter2*holeCenter4+holeCenter1*holeCenter4*holeCenter3+holeCenter4*holeCenter2*holeCenter3...
    +(holeCenter1*holeCenter2+holeCenter1*holeCenter3+holeCenter2*holeCenter3+holeCenter1*holeCenter4+holeCenter4*holeCenter3+holeCenter2*holeCenter4)*holeCenter5)/3;
a4=(holeCenter1*holeCenter2*holeCenter3*holeCenter4+(holeCenter1*holeCenter2*holeCenter3+holeCenter1*holeCenter2*holeCenter4+holeCenter1*holeCenter4*holeCenter3...
    +holeCenter4*holeCenter2*holeCenter3)*holeCenter5)/2;
a5=-holeCenter1*holeCenter2*holeCenter3*holeCenter4*holeCenter5;
fDeformInt=@(x) 1/6*x.^6+a1*x.^5+a2*x.^4+a3*x.^3+a4*x.^2+a5*x;


CtestDeformMeshXInt=fDeformInt(CX);
testDeformMeshXInt=fC2R(CtestDeformMeshXInt);


% hold on;
% plot(real(fDeformInt(holeCenter1)),imag(fDeformInt(holeCenter1)),'r','Marker','o','MarkerSize',10);
% plot(testDeformMeshXInt(BdIdx(:,1),1),testDeformMeshXInt(BdIdx(:,1),2),'r','Marker','o','MarkerSize',10);


% testDeformMeshXInt=[testDeformMeshXInt,X(:,2)];
% TRDeformInt=triangulation(T,testDeformMeshXInt);
% trimesh(TRDeformInt);

% scal1=1;
% X=scal1*X;
% CX=scal1*CX;
% CtestDeformMeshXInt=scal1*CtestDeformMeshXInt;
% testDeformMeshXInt=scal1*testDeformMeshXInt;

figure(2);
TRDeformInt=triangulation(T,testDeformMeshXInt);
p=triplot(TRDeformInt,'Color','#D3D3D3');
%p.FaceVertexCData
hold on;
drawBoundary(BD,testDeformMeshXInt);
argHole0=argComp(BD,testDeformMeshXInt);



%% optimal
P2PVtxIds=[];
%CtestDeformMeshXInt=0.25*CtestDeformMeshXInt;
%[XP2PDeform, statsAll] = meshAQP(CX, T, P2PVtxIds, CtestDeformMeshXInt(P2PVtxIds), CtestDeformMeshXInt, 1000);
[XP2PDeform,triEn1, statsAll]=meshNewton(CX,T,P2PVtxIds, CtestDeformMeshXInt(P2PVtxIds),CtestDeformMeshXInt, 1000,100,'SymmDirichlet');
%[XP2PDeform, statsAll] = meshAQP(CX, T, P2PVtxIds, XP2PDeform(P2PVtxIds), XP2PDeform, 1000);
load dat flipInd2

XP2PDeform=fC2R(XP2PDeform);
XP2PDeform=[XP2PDeform,X(:,1)];
argHole1=argComp(BD,XP2PDeform(:,1:2));
%% 
figure(3);
triplot(T,XP2PDeform(:,1),XP2PDeform(:,2),'Color','#D3D3D3');
hold on;
%plot(XP2PDeform(BdIdx(:,1),1),XP2PDeform(BdIdx(:,1),2),'r','Marker','o','MarkerSize',2);
drawBoundary(BD,[XP2PDeform(:,1),XP2PDeform(:,2)]);



fprintf('\n\n\n\n\n\n\n\n\n\n');

%% comparison
direct_scalCX=0.2*CX;
%direct_scalCX=0.5*CX;
[XP2PDeform2,triEn2, statsAll]=meshNewton(CX,T,P2PVtxIds, direct_scalCX(P2PVtxIds),direct_scalCX, 1000,100,'SymmDirichlet');

XP2PDeform2=fC2R(XP2PDeform2);

%%
figure(4);
triplot(T,XP2PDeform2(:,1),XP2PDeform2(:,2),'Color','#D3D3D3');
hold on;
drawBoundary(BD,XP2PDeform2);
argHole2=argComp(BD,XP2PDeform2);

%% draw energy map
figure(3);
trisurf(T,XP2PDeform(:,1),XP2PDeform(:,2),X(:,1),triEn1,'EdgeColor','none');
caxis([0 max(triEn2)]);
colormap('turbo');
colorbar;
hold on;
drawBoundary(BD,[XP2PDeform,X(:,1)]);

figure(4);
trisurf(T,XP2PDeform2(:,1),XP2PDeform2(:,2),zeros(size(XP2PDeform2,1),1),triEn2,'EdgeColor','none');
caxis([0 max(triEn2)]);
colormap('turbo');
colorbar;
hold on;
drawBoundary(BD,[XP2PDeform,zeros(size(XP2PDeform2,1),1)]);

m1=max(triEn1);
m2=max(triEn2);

t1=sum(triEn1);
t2=sum(triEn2);


%% flip
%[flipInd1,k]=testFlip(X,T,XP2PDeform);
k2=find(flipInd2);
if ~isempty(k2)
    figure(1);
    triplot(T,X(:,1),X(:,2));
    hold on;
    plot(X(T(k2(1),:),1),X(T(k2(1),:),2),'r','LineWidth',10);

    TRDeformInt=triangulation(T,testDeformMeshXInt);
    figure(2);
    triplot(TRDeformInt);
    hold on;
    plot(testDeformMeshXInt(T(k2(1),:),1),testDeformMeshXInt(T(k2(1),:),2),'r','LineWidth',10);
    %trimesh(TRDeformInt);
    figure(3);
    V=vertexAttachments(TRDeformInt,T(k2(1),:)');
    triplot(TRDeformInt);
    hold on  
    triplot(TRDeformInt(V{1},:),testDeformMeshXInt(:,1),testDeformMeshXInt(:,2),'Color','r','LineWidth',10);
end
