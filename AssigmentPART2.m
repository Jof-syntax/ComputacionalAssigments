classdef AssigmentPART2 < handle
    
    properties
        displacements
        Uexact
        Nelems
    end
    
    methods (Access = public)
        
        function obj = AssigmentPART2(data)
            obj.Nelems = obj.computeNElems();
            obj.displacements = obj.computeDisplacments(data);
            Error = obj.computeError(1); %compute error for U-Ue
            Error2 = obj.computeError(2); %compute error for U-Ue
            Error3 = obj.computeError(3); %compute error for U-Ue
            Error4 = obj.computeError(5); %compute error for U-Ue
        end
        
    end
    
    methods (Access = private)
        
        function Nelems = computeNElems(obj)
            Nelems = [4 8 16 32 64];
        end
        
        function displacements = computeDisplacments(obj, data)
            Nelems = obj.Nelems;
            a = length(Nelems);
            for i = 1:1:a
                b = Nelems(i);
                displacements(i) = Part6(b, data);
            end
        end
        
        function Error = computeError(obj, currentNelem)
            Nelems = obj.Nelems;
            a = length(Nelems);
            coord_approx = obj.displacements(currentNelem).COOR;
            coord_exact = obj.displacements(a).COOR;
            d_approx = obj.displacements(currentNelem).displacement;
            d_exact = obj.displacements(a).displacement;
            Error = 0;
            for i = 1:1:length(Nelems(currentNelem))
                coordx1 = coord_approx(i);
                coordx2 = coord_approx(i+1);
                found = 0;
                k = 1;
                while found == 0
                    if coordx2 == coord_exact(k)
                        found = 1;
                    else
                        k = k+1;
                    end
                end
                j = 1;
                found = 0;
                while found == 0
                    if coordx1 == coord_exact(j)
                        found = 1;
                    else
                        j = j+1;
                    end
                end
                displacement1 = d_approx(i);
                displacement2 = d_approx(i+1);
                displacement1Exact = d_exact(j);
                displacement2Exact = d_exact(k);
                Error = Error + obj.Gaussquadrature(displacement1, displacement2, displacement1Exact, displacement2Exact,  coordx1, coordx2);
            end
            Error = sqrt(Error);
        end
        
    end
    
    methods (Access = private, Static)
        
        function Error = Gaussquadrature(displacement1, displacement2, displacement1Exact, displacement2Exact,  coordx1, coordx2)
            xiG(1) = sqrt(3/5);
            xiG(2) = -sqrt(3/5);
            xiG(3) = 0;
            w(1) = 5/9;
            w(2) = 5/9;
            w(3) = 8/9;
            qFun1 = ((1-xiG)*(displacement1Exact - displacement1) + (1+xiG)*(displacement2Exact - displacement2)).^2;
            int = w*qFun1';
            he = coordx2 - coordx1;
            Error = 1/8*he*int;
        end
        
    end
end

