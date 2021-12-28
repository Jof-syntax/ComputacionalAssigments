clc
clear
load('INFO_FE.mat');        %these loads contain some useful information
load('dataP4.mat');         %such as the COOR and DORr arrays or the stiffness and mass matrices
neig = 25;                  %number of modes

nnode = size(COOR,1); 
ndim = size(COOR,2);
DOFl = 1:nnode*ndim ;
DOFl(DOFr) = [] ; 

Mll = M(DOFl,DOFl);   % Not restricted matrices of 'M' and 'K'
Kll = K(DOFl,DOFl);

[MODES, FREQ] = UndampedFREQ(Mll,Kll,neig);         %using this function
                                                    %we obtain the natural
                                                    %frequencies and modes

GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,NameFileMesh,DATA,DOFl); 