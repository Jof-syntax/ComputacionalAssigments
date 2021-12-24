function [d,strainGLO,stressGLO,React,posgp] = SolveELAS(K,Fb,Ftrac,dR,DOFr,COOR,CN,TypeElement,celasglo,typePROBLEM,celasgloINV,DATA) ;
% This function returns   the (nnode*ndim x 1) vector of nodal displacements (d),
% as well as the arrays of stresses and strains
%%% points  (qheatGLO)
% Input data
% K = Global stiffness matrix   (nnode*ndim x nnode*ndim)
% Fb = External force vector due to  body forces  (nnode*ndim x 1)
% Ftrac = External force vector due to  boundary tractions    (nnode*ndim x 1)
% DOFr = Set of restricted DOFs
% dR = Vector of prescribed displacements
% ----------------------
if nargin == 0
    load('tmp.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;     %
% Solution of the system of FE equation
% Right-hand side
F = Fb + Ftrac ;
% Set of nodes at which temperature is unknown
DOFl = 1:nnode*ndim ;
DOFl(DOFr) = [] ;


Klr = K(DOFl,DOFr);
Kll = K(DOFl,DOFl);
Krl = K(DOFr,DOFl);
Krr = K(DOFr,DOFr);
Fl = F(DOFl,1);
Fr = F(DOFr,1);

dL = Kll\(Fl-Klr*dR);
Rr = Krl*dL + Krr*dR - Fr;


d = zeros(nnode*ndim,1) ; % Nodal displacements (initialization)
React = zeros(size(d)) ;  %  REaction forces  (initialization)

d(DOFl)= dL ; 
d(DOFr) = dR ;

React(DOFl)= 0 ; 
React(DOFr) = Rr ;

%codi apartat 3
[Rx,Ry,Rz,Mx,My,Mz] = codi_apartat3(DOFr,Rr,COOR);
  

%%%% COmputation of strain and stress vector at each gauss point
disp('Computation of stress and strains at each Gauss point')
[strainGLO,stressGLO,posgp]= StressStrains(COOR,CN,TypeElement,celasglo,d,typePROBLEM,celasgloINV,DATA) ;

