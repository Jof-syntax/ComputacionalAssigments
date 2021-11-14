function [weig,posgp,shapef,dershapef] = Quadrilateral4NInPoints()
weig  = [1 1 1 1] ;
posgp = 1/sqrt(3)*[-1 -1;
                    1 -1;
                    1  1;
                   -1  1];
ndim = 2;
nnodeE = 4;
ngaus = length(weig) ;
shapef = zeros(ngaus, nnodeE) ;
dershapef = zeros(ndim, nnodeE, ngaus) ;
for g = 1:length(weig)
    xi = posgp(g, 1);
    eta = posgp(g, 2);
    Ne = 1/4*[(1-xi)*(1-eta)  (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];
    BeXi = 1/4*[-(1-eta)  (1-eta)  (1+eta)  -(1+eta);
        -(1-xi)  -(1+xi)   (1+xi)    (1-xi)] ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end
end





