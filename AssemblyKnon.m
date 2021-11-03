
function K = AssemblyKnon(COOR,CN,d_k,AreaFUN,DerStressFUN)

nElem = size(CN, 1 );
nNode = size(COOR, 1 );
nNodeE = size(CN, 2 );
K = sparse( nNode , nNode);

xiG = [sqrt(3/5), -sqrt(3/5), 0];
w = [5/9, 5/9, 8/9];

for e = 1:1:nElem
    
    NODOSe = CN( e , : );
    COOR_e = COOR(NODOSe );
    he = COOR_e( 2 ) - COOR_e( 1 );
    Be = 1/he*[-1 1];
    
    d1 = d_k(e);
    d2 = d_k(e+1);
    de =  [d1; d2];
    strain = Be*de;%*Be'
    Ee = DerStressFUN(strain);
    
    Ke = 0;
    for i = 1:1:length(w) %Gauss integration (3 points)
        N1 = (1-xiG(i))/2;
        N2 = (1+xiG(i))/2;
        x = N1*COOR_e( 1 )+N2*COOR_e( 2 );
        areaE = AreaFUN(x);
        f = Be'.*(areaE*Ee)*Be;
        Ke = Ke + w(i)*f*he/2;
    end
    
    for a = 1:1:nNodeE
        for b = 1:1:nNodeE
            A = CN(e,a);
            B = CN(e,b);
            K(A,B) = K(A,B) + Ke(a,b);
        end
    end
end
end