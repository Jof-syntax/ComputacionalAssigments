function [Rx,Ry,Rz,Mx,My,Mz] = codi_apartat3(DOFr,Rr,COOR)
Rx = 0;
Ry = 0;
Rz = 0;
Mx = 0;
My = 0;
Mz = 0;
CG = [0 0.125 -0.125]; %CG respecte els eixos del nodes (COOR)

j=1;
for i = 1:1:length(COOR)
    if COOR(i,1) == 0
        distancia(j,:) = COOR(i,:) - CG;
        j = j+1;
    end
end

%funcions - cross, reshape, sum

A = reshape(Rr',[length(DOFr)/3,3]);
R = A';
S = sum(R,2);

Rx = S(1,1);
Ry = S(2,1);
Rz = S(3,1);

Mom = cross(distancia,A);

SS = sum(Mom',2);

Mx = SS(1,1);
My = SS(2,1);
Mz = SS(3,1);


end