%% start
working_dataset='test';
fC2R = @(x) [real(x) imag(x)];
fR2C = @(x) complex(x(:,1), x(:,2));
datadir = fullfile(cd, 'data', working_dataset, '\');
imgfilepath  = fullfile(datadir, 'image.png');

%% test for 2 holes ring mesh
arg1=(60:300)/180*pi;
arg2=(-120:120)/180*pi;
arg3=(0:360)/180*pi;
holeCenter1=-1;
holeCenter2=1;
Ring1=[2*cos(arg1)'-1,2*sin(arg1)'];
Ring2=[2*cos(arg2)'+1,2*sin(arg2)'];
outerBoundary=[Ring1;Ring2];
holes1=[0.5*cos(arg3)'-1,0.5*sin(arg3)'];
holes2=[0.5*cos(arg3)'+1,0.5*sin(arg3)'];

[X, T] = cdt([{outerBoundary},{holes1,holes2}],[], 100, false);

 CX=fR2C(X);
 TR=triangulation(T,X);
 
%% integration
