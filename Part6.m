classdef Part6
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = Part6(inputArg1,inputArg2)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function K = AssemblyK(COOR,CN)
            nelem = size(CN, 1 );
            nnode = size(COOR, 1 );
            nnodeE = size(CN, 2 );
            K = sparse( nnode , nnode );
            for e =1: nelem
                NODOSe = CN( e , : ) ; COOR_e = COOR(NODOSe ) ;
                he = COOR_e( 2 )-COOR_e( 1 ) ;
                Ke = 1/he *[ 1 -1; -1 1] ; % algu...;
                for a = 1: nnodeE
                    for b = 1:nnodeE
                        A = CN(e,a);
                        B = CN(e,b);
                        K(A,B) = K(A,B) + Ke(a,b);
                    end
                end
            end
        end
    end
end

