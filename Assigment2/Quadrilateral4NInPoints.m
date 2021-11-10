function [weig,posgp,shapef,dershapef] = Quadrilateral4NInPoints

weig  = [1 1 1 1] ;
posgp(1) = 1/sqrt(3)*[-1 -1 ] ;
posgp(2) = 1/sqrt(3)*[1 -1 ] ;
posgp(3) = 1/sqrt(3)*[1 1 ] ;
posgp(4) = 1/sqrt(3)*[-1 1 ] ;


ndim = 1;
nnodeE = 2 ;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig) 
    xi = posgp(g) ;
    eta = posgp(g) ;%que es eta?
    Ne = 1/4*[(1-xi)*(1-eta)  (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];
    BeXi = 1/4*[-(1-eta)  (1-eta)  (1+eta)  -(1+eta);
                -(1-xi)  -(1+xi)   (1+xi)    (1-xi)] ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end
end





