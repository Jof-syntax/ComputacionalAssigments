function [weig,posgp,shapef,dershapef] = Hexahedra8NInPoints()
weig  = [1 1 1 1 1 1 1 1] ;
posgp = (1/sqrt(3))*[-1 -1 -1;
                    1 -1 -1;
                    1  1 -1;
                   -1  1 -1;
                   -1 -1 1;
                    1 -1 1;
                    1  1 1;
                   -1  1 1];
ndim = 3;
nnodeE = 8;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g = 1:length(weig)
    xi = posgp(g, 1);
    eta = posgp(g, 2);
    c = posgp(g, 3);
    Ne = (1/8)*[(1-xi)*(1-eta)*(1-c)  (1+xi)*(1-eta)*(1-c) (1+xi)*(1+eta)*(1-c) (1-xi)*(1+eta)*(1-c) (1-xi)*(1-eta)*(1+c)  (1+xi)*(1-eta)*(1+c) (1+xi)*(1+eta)*(1+c) (1-xi)*(1+eta)*(1+c)];
    BeXi = (1/8)*[  -(1-eta)*(1-c)      (1-eta)*(1-c)    (1+eta)*(1-c)   -(1+eta)*(1-c)     -(1-eta)*(1+c)      (1-eta)*(1+c)    (1+eta)*(1+c)   -(1+eta)*(1+c);
                    -(1-xi)*(1-c)      -(1+xi)*(1-c)     (1+xi)*(1-c)     (1-xi)*(1-c)      -(1-xi)*(1+c)      -(1+xi)*(1+c)     (1+xi)*(1+c)     (1-xi)*(1+c);
                    -(1-xi)*(1-eta)    -(1-eta)*(1+xi)  -(1+eta)*(1+xi)  -(1+eta)*(1-xi)    (1-xi)*(1-eta)      (1-eta)*(1+xi)   (1+eta)*(1+xi)   (1+eta)*(1-xi);] ; 
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end
end