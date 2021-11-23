function [qheatGLO]= HeatFlux(COOR,CN,TypeElement,ConductMglo,d)

nnodeE = size(CN,2) ; %Number of nodes per element
nelem = size(CN,1);   % Number of elements
nnode = size(COOR,1);  % Number of nodes
ndim = size(COOR,2);   % Spatial Dimension of the problem  (2 or 3)

TypeIntegrand = 'RHS'; 

[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ; 

ngaus = length(weig);
qheatGLO = zeros(ngaus*ndim,nelem);


for e = 1:nelem
    
    CN_elem = CN(e,:);
    Xe = (COOR(CN_elem,:))';
    ConductM = ConductMglo(:,:,e);
    ii = 0;
    for  g = 1:ngaus
      % Matrix of derivatives for Gauss point "g"
      BeXi = dershapef(:,:,g) ; 
      % Jacobian Matrix 
      Je = Xe*BeXi' ;
      % Matrix of derivatives with respect to physical coordinates 
      Be = inv(Je)'*BeXi ;
      
      N_u = Be*d(CN_elem,1);
      multiplicacio = ConductM*N_u;
      ii = ii+1;
      qheatGLO(ii,e) = qheatGLO(ii,e) - multiplicacio(1,1);
      ii = ii+1;
      qheatGLO(ii,e) = qheatGLO(ii,e) - multiplicacio(2,1);
    end
end


end