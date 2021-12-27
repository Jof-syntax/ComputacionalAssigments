function M = ComputeM(COOR,CN,TypeElement, density)
% This subroutine   returns the global stiffness matrix K (ndim*nnode x ndim*nnode)
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CN: Connectivity matrix (nelem x nnodeE), % TypeElement: Type of finite element (quadrilateral,...),  celasglo (nstrain x nstrain x nelem)  % Array of elasticity matrices
% Dimensions of the problem
if nargin == 0
    load('tmp1.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;  
nstrain = size(density,1) ;
% Shape function routines (for calculating shape functions and derivatives)
TypeIntegrand = 'K';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
% Assembly of matrix K
% ----------------
M = sparse([],[],[],nnode*ndim,nnode*ndim,nnodeE*ndim*nelem) ;
for e = 1:nelem
    den = density(e) ;  % Density matrix of element "e"
    CNloc = CN(e,:) ;   % Coordinates of the nodes of element "e"
    Xe = COOR(CNloc,:)' ;     % Computation of elemental stiffness matrix
    Me = ComputeMeMatrix(den,weig,dershapef,Xe) ;
    
    for anod = 1:1:nnodeE
        a = Nod2DOF(anod,ndim) ;
        Anod = CN(e,anod) ; 
        A = Nod2DOF(Anod,ndim) ; 
        for bnod = 1:1:nnodeE
            b = Nod2DOF(bnod,ndim) ;
            Bnod = CN(e,bnod) ; 
            B = Nod2DOF(Bnod,ndim) ;
            M(A,B) = M(A,B) + Me(a,b);
        end
    end
   
    
  
    
    if mod(e,10)==0  % To display on the screen the number of element being assembled
        disp(['e=',num2str(e)])
    end



end