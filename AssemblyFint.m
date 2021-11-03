function [Residual,STRAIN,STRESS] = AssemblyFint(COOR,CN,d_k,stressFUN,AreaFUN)
nElem = size(CN, 1);
nnodeE = size(CN, 2);
xiG = [sqrt(3/5), -sqrt(3/5), 0];
w = [5/9, 5/9, 8/9];

STRAIN = zeros(size(d_k)-1);
Fint = zeros(size(d_k));
for  e = 1:1:nElem
    %compute the coordinates of the element
    
    NODESe = CN( e , : );
    COOR_e = COOR(NODESe );
    he = COOR_e( 2 ) - COOR_e( 1 );

    % compute
    d1 = d_k(e);
    d2 = d_k(e+1);
    dE = [d1; d2];
    Be = 1/he*[-1 1];
    
    sigE = stressFUN(Be*dE);
    STRAIN(e) = Be*dE;
    
    Fe = 0;
    for i = 1:1:length(w) %Gauss integration (3 points)
        N1 = (1-xiG(i))/2; 
        N2 = (1+xiG(i))/2;
        x = N1*COOR_e( 1 )+N2*COOR_e( 2 );
        areaE = AreaFUN(x);
        
        f = Be*areaE*sigE;
        Fe = Fe + he/2*w(i)*f;
    end
    
    for a = 1:nnodeE
        A = CN(e,a);
        Fint(A) = Fint(A)+Fe(a);
    end
    
end
STRESS = stressFUN(STRAIN);
Residual = Fint;
end


