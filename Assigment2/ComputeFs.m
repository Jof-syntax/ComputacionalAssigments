function Fs = ComputeFs(COOR,CN,TypeElement, fNOD)
% This subroutine   returns the  heat source contribution (Fs)    to the
% global flux vector. Inputs
% --------------
% 1. Finite element mesh
% -------------------
% COOR: Coordinate matrix (nnode x ndim)
% CN: Connectivity matrix (nelem x nnodeE)
% TypeElement: Type of finite element (quadrilateral,...)
% -----------
% 2. Vector containing the values of the heat source function at the nodes
% of the mesh
% -----------
%  fNOD (nnode x 1)  %
%%%%

% Dimensions of the problem
nnode = size(COOR,1); 
ndim = size(COOR,2); 
nelem = size(CN,1); 
nnodeE = size(CN,2);
TypeIntegrand = 'RHS';
[ weig , posgp , shapef , dershapef ] = ComputeElementShapeFun( TypeElement , nnodeE , TypeIntegrand ) ;
Fs = zeros(nnode,1);
for e = 1:1:nelem
    CNe = CN( e , : );
    % Source function evaluated at the nodes of element " e "
    fe = fNOD(CNe) ;
    % Coordinates of the nodes of element " e "
    Xe = COOR(CNe, : )';
    % Computation of elemental source fluxvector
    Fse = ComputeFseVector(fe , weig , shapef , dershapef , Xe);
    for a = 1:nnodeE % Assembly of the Fint matrix for the element 'e'
        A = CN(e,a);
        Fs(A) = Fs(A)+Fi(a);
    end
end
end