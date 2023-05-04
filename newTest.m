%% start
working_dataset='test';
fC2R = @(x) [real(x) imag(x)];
fR2C = @(x) complex(x(:,1), x(:,2));
datadir = fullfile(cd, 'data', working_dataset, '\');
imgfilepath  = fullfile(datadir, 'image.png');


%% test for 3 holes ring mesh
arg1=(0:2:360)'/180*pi;
holeCenter1=3i/5;
holeCenter2=(1-0.6i)/2;
holeCenter3=(-1-0.6i)/2;


outerBoundary=[cos(arg1),sin(arg1)];
r=0.3;
arg2=(0:6:360)'/180*pi;
holes1=[r*cos(arg2),r*sin(arg2)]/2+fC2R(holeCenter1);
holes2=[r*cos(arg2),r*sin(arg2)]/2+fC2R(holeCenter2);
holes3=[r*cos(arg2),r*sin(arg2)]/2+fC2R(holeCenter3);

holeCenter4=0;
holes4=[r*cos(arg2),r*sin(arg2)]/2+fC2R(holeCenter4);
% plot(outerBoundary(:,1),outerBoundary(:,2),holes1(:,1),holes1(:,2),holes2(:,1),holes2(:,2),holes3(:,1),holes3(:,2));

%[X, T] = cdt([{outerBoundary},{holes1,holes2,holes3}],[], 10000, false);
[X, T] = cdt([{outerBoundary},{holes1,holes2,holes3,holes4}],[], 10000, false);

figure(1);
CX=fR2C(X);
TR=triangulation(T,X);
triplot(TR);
BD=getBoundary(TR);
nBD=size(BD,2);
argHole=argComp(BD,X);




%% integration
%3 holes
% a1=-(holeCenter1+holeCenter2+holeCenter3)/3;
% a2=(holeCenter1*holeCenter2+holeCenter1*holeCenter3+holeCenter2*holeCenter3)/2;
% a3=-holeCenter1*holeCenter2*holeCenter3;
% fDeformInt=@(x) 1/4*x.^4+a1*x.^3+a2*x.^2+a3*x;

%4 holes
a1=-(holeCenter1+holeCenter2+holeCenter3+holeCenter4)/4;
a2=(holeCenter1*holeCenter2+holeCenter1*holeCenter3+holeCenter2*holeCenter3+holeCenter1*holeCenter4+holeCenter4*holeCenter3+holeCenter2*holeCenter4)/3;
a3=-(holeCenter1*holeCenter2*holeCenter3+holeCenter1*holeCenter2*holeCenter4+holeCenter1*holeCenter4*holeCenter3+holeCenter4*holeCenter2*holeCenter3)/2;
a4=holeCenter1*holeCenter2*holeCenter3*holeCenter4;
fDeformInt=@(x) 1/5*x.^5+a1*x.^4+a2*x.^3+a3*x.^2+a4*x;

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
triplot(TRDeformInt);
hold on;
drawBoundary(BD,testDeformMeshXInt);
argHole0=argComp(BD,testDeformMeshXInt);



%% optimal
P2PVtxIds=[];
CtestDeformMeshXInt=0.6*CtestDeformMeshXInt;
%[XP2PDeform, statsAll] = meshAQP(CX, T, P2PVtxIds, CtestDeformMeshXInt(P2PVtxIds), CtestDeformMeshXInt, 1000);
[XP2PDeform, statsAll]=meshNewton(CX,T,P2PVtxIds, CtestDeformMeshXInt(P2PVtxIds),CtestDeformMeshXInt, 1000,100,'SymmDirichlet');
%[XP2PDeform, statsAll] = meshAQP(CX, T, P2PVtxIds, XP2PDeform(P2PVtxIds), XP2PDeform, 1000);
load dat flipInd2

XP2PDeform=fC2R(XP2PDeform);
XP2PDeform=[XP2PDeform,X(:,2)];
figure(3);
triplot(T,XP2PDeform(:,1),XP2PDeform(:,2));
hold on;
%plot(XP2PDeform(BdIdx(:,1),1),XP2PDeform(BdIdx(:,1),2),'r','Marker','o','MarkerSize',2);
drawBoundary(BD,XP2PDeform);
argHole1=argComp(BD,XP2PDeform(:,1:2));

% trimesh(T,XP2PDeform(:,1),XP2PDeform(:,2),X(:,2));
fprintf('\n\n\n\n\n\n\n\n\n\n');

%% comparison
direct_scalCX=0.2*CX;
[XP2PDeform2, statsAll]=meshNewton(CX,T,P2PVtxIds, direct_scalCX(P2PVtxIds),direct_scalCX, 1000,100,'SymmDirichlet');

XP2PDeform2=fC2R(XP2PDeform2);

figure(4);
triplot(T,XP2PDeform2(:,1),XP2PDeform2(:,2));
hold on;
drawBoundary(BD,XP2PDeform2);
argHole2=argComp(BD,XP2PDeform2);




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
