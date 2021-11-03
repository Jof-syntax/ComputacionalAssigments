function [Residual,STRAIN,STRESS] = AssemblyFint(COOR,CN,d_k,stressFUN,AreaFUN)

nElem = size(CN, 1);
nNodeE = size(CN, 2);

% Gauss integration parameters
xiG = [sqrt(3/5), -sqrt(3/5), 0];
w = [5/9, 5/9, 8/9];

STRAIN = zeros(nElem, 1); %Creates strain vector
STRESS = zeros(nElem, 1); %Creates stress vector
Fint = zeros(size(d_k)); %Creates internal forces vector

for  e = 1:1:nElem % For each element:
    % Obtain the coordinates, of the element
    NODESe = CN( e , : );
    COOR_e = COOR(NODESe );
    
    he = COOR_e( 2 ) - COOR_e( 1 ); % He defined between two coordinates
    Be = 1/he*[-1 1]; % Be definition

    % Displacement of the nodes of the element 
    d1 = d_k(e);
    d2 = d_k(e+1);
    dE = [d1; d2];
    
    STRAIN(e) = Be*dE; % Computes the strain for the element 'e'

    STRESS(e) = stressFUN(STRAIN(e));  % Computes the stress for the element 'e'

    Fi = 0;
    for i = 1:1:length(w) %Gauss integration (3 points)
        N1 = (1-xiG(i))/2; 
        N2 = (1+xiG(i))/2;
        x = N1*COOR_e( 1 )+N2*COOR_e( 2 );
        areaE = AreaFUN(x); % Area of the element
        
        f = he/2*Be*areaE*STRESS(e); % See demostration of this equation in the report
        Fi = Fi + w(i)*f; % Internal forces for element 'e'
    end
    
    for a = 1:nNodeE % Assembly of the Fint matrix for the element 'e'
        A = CN(e,a);
        Fint(A) = Fint(A)+Fi(a); 
    end
    
end
Residual = Fint; % Associates the Internal forces with the residual
end


