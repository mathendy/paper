%% start
working_dataset = 'annulus';
deformer_main;

%% test for existing mesh
Y=fC2R(X);
holeCenter1=[3,0];
holeCenter2=[-3,0];
C1=fR2C(holeCenter1);
C2=fR2C(holeCenter2);
triplot(T,Y(:,1),Y(:,2))
hold on;
plot(holeCenter1(1),holeCenter1(2),'o',holeCenter2(1),holeCenter2(2),'o');
hold off;

X_test=(X-C1).*(X-C2);
Y_test=fC2R(X_test);
triplot(T,Y_test(:,1),Y_test(:,2))

Z1=Y(T(:,2),:)-Y(T(:,1),:);Z2=Y(T(:,3),:)-Y(T(:,1),:);
Z1=[Z1,zeros(size(T,1),1)];Z2=[Z2,zeros(size(T,1),1)];
Z=cross(Z1,Z2);
sigZ=sign(Z(:,3));

Z1_test=Y_test(T(:,2),:)-Y_test(T(:,1),:);Z2_test=Y_test(T(:,3),:)-Y_test(T(:,1),:);
Z1_test=[Z1_test,zeros(size(T,1),1)];Z2_test=[Z2_test,zeros(size(T,1),1)];
Z_test=cross(Z1_test,Z2_test);
sigZ_test=sign(Z_test(:,3));

k=sum((sigZ-sigZ_test)~=0);

%% test for 2 holes ring (only boundary)
arg1=(60:300)/180*pi;arg2=(-120:120)/180*pi;
Ring1=[2*cos(arg1)'-1,2*sin(arg1)'];
Ring2=[2*cos(arg2)'+1,2*sin(arg2)'];
Ring3=[0.5*cos((0:360)/180*pi)'-1,0.5*sin((0:360)/180*pi)'];
Ring4=[0.5*cos((0:360)/180*pi)'+1,0.5*sin((0:360)/180*pi)'];
Ring=[Ring1;Ring2;Ring3;Ring4];
plot(Ring1(:,1),Ring1(:,2),Ring2(:,1),Ring2(:,2),Ring3(:,1),Ring3(:,2),Ring4(:,1),Ring4(:,2))

CRing1=fR2C(Ring1);CRing2=fR2C(Ring2);CRing3=fR2C(Ring3);CRing4=fR2C(Ring4);
CRing_test1=CRing1.*CRing1-1;CRing_test2=CRing2.*CRing2-1;CRing_test3=CRing3.*CRing3-1;CRing_test4=CRing4.*CRing4-1;
Ring_test1=fC2R(CRing_test1);Ring_test2=fC2R(CRing_test2);Ring_test3=fC2R(CRing_test3);Ring_test4=fC2R(CRing_test4);
%plot(Ring_test1(:,1),Ring_test1(:,2));
%plot(Ring_test2(:,1),Ring_test2(:,2));
%plot(Ring_test3(:,1),Ring_test3(:,2));
%plot(Ring_test4(:,1),Ring_test4(:,2));
plot(Ring_test1(:,1),Ring_test1(:,2),Ring_test2(:,1),Ring_test2(:,2),Ring_test3(:,1),Ring_test3(:,2),Ring_test4(:,1),Ring_test4(:,2))

%% test for 2 holes ring mesh
[testMeshX, testMeshT] = cdt([{[Ring1;Ring2]},{Ring3,Ring4}],[], 10000, false);
triplot(testMeshT,testMeshX(:,1),testMeshX(:,2));

CtestMeshX=fR2C(testMeshX);
CtestDeformMeshX=CtestMeshX.*CtestMeshX-1;
testDeformMeshX=fC2R(CtestDeformMeshX);
triplot(testMeshT,testDeformMeshX(:,1),testDeformMeshX(:,2));


%% test for 3 holes ring mesh
arg=(0:360)'/180*pi;
holeCenter1=1i;
holeCenter2=0.8-0.5i;
holeCenter3=-0.8-0.5i;
outerBoundary=[2*cos(arg),2*sin(arg)];
holes1=[0.6*cos(arg),0.6*sin(arg)+1];
holes2=[0.6*cos(arg)+0.8,0.6*sin(arg)-0.5];
holes3=[0.6*cos(arg)-0.8,0.6*sin(arg)-0.5];
%plot(outerBoundary(:,1),outerBoundary(:,2),holes1(:,1),holes1(:,2),holes2(:,1),holes2(:,2),holes3(:,1),holes3(:,2));

[testMeshX, testMeshT] = cdt([{outerBoundary},{holes1,holes2,holes3}],[], 100, false);
triplot(testMeshT,testMeshX(:,1),testMeshX(:,2));

CtestMeshX=fR2C(testMeshX);
CtestDeformMeshX=(CtestMeshX-holeCenter1).*(CtestMeshX-holeCenter2).*(CtestMeshX-holeCenter3);
testDeformMeshX=fC2R(CtestDeformMeshX);
%triplot(testMeshT,testDeformMeshX(:,1),testDeformMeshX(:,2));
TR=triangulation(testMeshT,testMeshX);
TRDeform=triangulation(testMeshT,testDeformMeshX);

%% test orientation
Z1=testMeshX(testMeshT(:,2),:)-testMeshX(testMeshT(:,1),:);Z2=testMeshX(testMeshT(:,3),:)-testMeshX(testMeshT(:,1),:);
Z1=[Z1,zeros(size(testMeshT,1),1)];Z2=[Z2,zeros(size(testMeshT,1),1)];
Z=cross(Z1,Z2);
sigZ=sign(Z(:,3));

Z1_deform=testDeformMeshX(testMeshT(:,2),:)-testDeformMeshX(testMeshT(:,1),:);Z2_deform=testDeformMeshX(testMeshT(:,3),:)-testDeformMeshX(testMeshT(:,1),:);
Z1_deform=[Z1_deform,zeros(size(testMeshT,1),1)];Z2_deform=[Z2_deform,zeros(size(testMeshT,1),1)];
Z_deform=cross(Z1_deform,Z2_deform);
sigZ_deform=sign(Z_deform(:,3));

k=sum((sigZ-sigZ_deform)~=0);
flipInd=find(sigZ-sigZ_deform);
figure(1);
triplot(testMeshT,testMeshX(:,1),testMeshX(:,2));
hold on;
plot(testMeshX(testMeshT(flipInd(1),:),1),testMeshX(testMeshT(flipInd(1),:),2),'r','LineWidth',10);
figure(2);
triplot(testMeshT,testDeformMeshX(:,1),testDeformMeshX(:,2));
hold on;
plot(testDeformMeshX(testMeshT(flipInd(1),:),1),testDeformMeshX(testMeshT(flipInd(1),:),2),'r','LineWidth',10);

figure(3);
V=vertexAttachments(TRDeform,testMeshT(flipInd(1),:)');
triplot(TRDeform);
hold on  
triplot(TRDeform(V{1},:),testDeformMeshX(:,1),testDeformMeshX(:,2),'Color','r','LineWidth',10);

%% integration
a1=-(holeCenter1+holeCenter2+holeCenter1)/3;
a2=(holeCenter1*holeCenter2+holeCenter1*holeCenter3+holeCenter2*holeCenter3)/2;
a3=-holeCenter1*holeCenter2*holeCenter3;

fDeformInt=@(x) 1/4*x.^4+a1*x.^3+a2*x.^2+a3*x;
CtestDeformMeshXInt=fDeformInt(CtestMeshX);
testDeformMeshXInt=fC2R(CtestDeformMeshXInt);
TRDeformInt=triangulation(testMeshT,testDeformMeshXInt);
triplot(TRDeformInt);

%% test integration deform orientation
Z1_deformInt=testDeformMeshXInt(testMeshT(:,2),:)-testDeformMeshXInt(testMeshT(:,1),:);Z2_deformInt=testDeformMeshXInt(testMeshT(:,3),:)-testDeformMeshXInt(testMeshT(:,1),:);
Z1_deformInt=[Z1_deformInt,zeros(size(testMeshT,1),1)];Z2_deformInt=[Z2_deformInt,zeros(size(testMeshT,1),1)];
Z_deformInt=cross(Z1_deformInt,Z2_deformInt);
sigZ_deformInt=sign(Z_deformInt(:,3));

k=sum((sigZ-sigZ_deform)~=0);
flipInd=find(sigZ-sigZ_deformInt);
figure(1);
triplot(testMeshT,testMeshX(:,1),testMeshX(:,2));
hold on;
plot(testMeshX(testMeshT(flipInd(1),:),1),testMeshX(testMeshT(flipInd(1),:),2),'r','LineWidth',10);
figure(2);
triplot(testMeshT,testDeformMeshXInt(:,1),testDeformMeshXInt(:,2));
hold on;
plot(testDeformMeshXInt(testMeshT(flipInd(1),:),1),testDeformMeshXInt(testMeshT(flipInd(1),:),2),'r','LineWidth',10);

figure(3);
V=vertexAttachments(TRDeform,testMeshT(flipInd(1),:)');
triplot(TRDeform);
hold on  
triplot(TRDeform(V{1},:),testDeformMeshX(:,1),testDeformMeshX(:,2),'Color','r','LineWidth',10);
