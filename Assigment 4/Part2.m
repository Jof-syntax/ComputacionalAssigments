clc
clear
load('INFO_FE.mat');
load('dataP4.mat');
neig = 25;

nnode = size(COOR,1); 
ndim = size(COOR,2);
DOFl = 1:nnode*ndim ;
DOFl(DOFr) = [] ; 

Mll = M(DOFl,DOFl);
Kll = K(DOFl,DOFl);

[MODES, FREQ] = UndampedFREQ(Mll,Kll,neig);

NameFileMesh = 'MyFirstMesh3D.msh';
GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DATA,DOFl);

