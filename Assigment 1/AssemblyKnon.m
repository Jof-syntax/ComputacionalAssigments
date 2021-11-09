function K = AssemblyKnon(COOR,CN,d_k,AreaFUN,DerStressFUN)

nElem = size(CN, 1 );
nNode = size(COOR, 1 );
nNodeE = size(CN, 2 );
K = sparse( nNode , nNode);  %Creates K matrix

% Gauss integration parameters
xiG = [sqrt(3/5), -sqrt(3/5), 0];
w = [5/9, 5/9, 8/9];

for e = 1:1:nElem % For each element:
    % Obtain the coordinates, of the element
    NODOSe = CN( e , : );
    COOR_e = COOR(NODOSe );
    
    he = COOR_e( 2 ) - COOR_e( 1 );  % He defined between two coordinates
    Be = 1/he*[-1 1];     % Be definition
    
    % Displacement of the nodes of the element 
    d1 = d_k(e);
    d2 = d_k(e+1);
    de =  [d1; d2];
    
    strain = Be*de; % By definition 
    Ee = DerStressFUN(strain); % 'E' of the element
    
    Ke = 0;
    for i = 1:1:length(w) % Gauss integration (3 points)
        N1 = (1-xiG(i))/2;
        N2 = (1+xiG(i))/2;
        x = N1*COOR_e( 1 )+N2*COOR_e( 2 );
        areaE = AreaFUN(x); % Area of the element
        f = Be'.*(areaE*Ee)*Be; % See demostration of this equation in the report
        Ke = Ke + w(i)*f*he/2; % K for the element 'e'
    end
    
    for a = 1:1:nNodeE % Assembly of the K matrix for the element 'e'
        for b = 1:1:nNodeE
            A = CN(e,a);
            B = CN(e,b);
            K(A,B) = K(A,B) + Ke(a,b);
        end
    end
end
end