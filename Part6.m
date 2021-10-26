classdef Part6 < handle
    
    properties (Access = public)
        displacement
        CN
        COOR
    end
    
    properties (Access = private)
        K
        Ff
        g
        L
    end
    
    methods (Access = public)
        
        function obj = Part6(Nelem, data)  %The constructor needs the values of parameters 'g', 'L' and the number of elements 'Nelem'  to work. This function obtains the displacements by means of 'solver'
            obj.g = data.g;
            obj.L = data.L;
            obj.CN = obj.computeCN(Nelem);
            obj.COOR = obj.computeCOOR(Nelem);
            obj.K = obj.AssemblyK();
            obj.Ff = obj.AssemblyF();
            obj.displacement = obj.solver(Nelem);
        end
        
        function plot(obj, time) % Creates a plot of the displacements obtained by function 'solver' 
            obj.createPlot(time);
        end
        
    end
    
    methods (Access = private)
        
        function CN = computeCN(obj, Nelem) %Computes the realation between nodes for a given number of elements 'Nelem'. 
            CN = sparse(Nelem, 2);
            for i = 1:1:Nelem
                CN(i, 1) = i;
                CN(i, 2) = i+1;
            end
        end
        
        function COOR = computeCOOR(obj, Nelem) % Computes the coordinates of the elements for a given number of elements 'Nelem'
            COOR = sparse(Nelem+1,1);
            increase = obj.L/Nelem;
            for i=1:1:Nelem
                COOR(i+1) = increase*i;
            end
        end
        
        
        function K = AssemblyK(obj) %Assembles the matrix K
            COOR = obj.COOR;
            CN = obj.CN;
            nElem = size(CN, 1 );
            nNode = size(COOR, 1 );
            nNodeE = size(CN, 2 );
            K = sparse( nNode , nNode);
            AdditionalTerm = obj.OneDGaussquadratureTermK();
            for e =1: nElem
                NODOSe = CN( e , : );
                COOR_e = COOR(NODOSe );
                he = COOR_e( 2 ) - COOR_e( 1 );
                Ke = - 1/he*[ 1 -1; -1 1] + he*AdditionalTerm; 
                for a = 1: nNodeE
                    for b = 1:nNodeE
                        A = CN(e,a);
                        B = CN(e,b);
                        K(A,B) = K(A,B) + Ke(a,b);
                    end
                end
            end
        end
        
        function Ff = AssemblyF(obj)  %Assembles the column vector F
            L = obj.L;
            g = obj.g;
            rho = pi^2/L^2;
            b = g*pi^2/L;
            COOR = obj.COOR;
            CN = obj.CN;
            nElem = size(CN, 1 );
            nNode = size(COOR, 1 );
            nNodeE = size(CN, 2 );
            Ff =zeros( nNode , 1);
            for e =1: nElem
                NODOSe = CN( e , : );
                COOR_e = COOR(NODOSe );
                he = COOR_e( 2 ) - COOR_e( 1 );
                Fe = he*obj.COMPUTE_Fe_FORCE(COOR_e( 1 ), COOR_e( 2 ));
                for a = 1:nNodeE
                    A = CN(e, a);
                    Ff(A) = Ff(A) + Fe(a);
                end
            end
            Ff(nNode, 1) = Ff(nNode, 1) - b;
        end
        
        function displacement = solver(obj, Nelem) %Solves the displacements of each node, using the matrices F and K
            Ff = obj.Ff;
            K = obj.K;
            R = 1;
            L = 1:Nelem+1;
            L(R) = [];
            KRR = K(R, R);
            KRL = K(R, L);
            KLR = K(L, R);
            KLL = K(L, L);
            dR = -obj.g; %Displacement of the first node (Known)
            FL = Ff(L,1);
            dL = KLL\(FL-KLR*dR); % Solution for the displacements of the other nodes
            displacement = sparse(Nelem +1);
            displacement(1) = dR;
            for i = 1:1:Nelem
                displacement(i+1) = dL(i);
            end
        end
        
        function AdditionalTermK = OneDGaussquadratureTermK(obj)
            L = obj.L;
            rho = pi^2/L^2;
            xiG(1) = sqrt(3/5);
            xiG(2) = -sqrt(3/5);
            xiG(3) = 0;
            w(1) = 5/9;
            w(2) = 5/9;
            w(3) = 8/9;
            qFun1 = (1-xiG).^2;
            int1 = w*qFun1';
            qFun2 = (1-xiG.^2);
            int2 = w*qFun2';
            qFun3 = (1+xiG).^2;
            int3 = w*qFun3';
            AdditionalTermK = rho*1/8*[ int1 int2 ;
                int2 int3];
        end
        
        function TermF = COMPUTE_Fe_FORCE(obj, coordx1, coordx2)
            L = obj.L;
            g = obj.g;
            rho = pi^2/L^2;
            s = g*rho^2;
            xiG(1) = 0.339981043584856;
            xiG(2) = -0.339981043584856;
            xiG(3) = 0.861136311594053;
            xiG(4) = -0.861136311594053;
            w(1) = 0.652145154862546;
            w(2) = 0.652145154862546;
            w(3) = 0.347854845137454;
            w(4) = 0.347854845137454;
            x = ((1-xiG)*coordx1+(1+xiG)*coordx2);
            qFun1 = (1-xiG).*x.^2;
            int1 = w*qFun1';
            qFun2 = (1+xiG).*x.^2;
            int2 = w*qFun2';
            TermF = s*1/16*[int1;
                int2];
        end
        
        function createPlot(obj, time) % Creates the plot and its configuration
            close all;
            figure;
            hold on;
            COOR = obj.COOR;
            Disp = obj.displacement;
            plot(COOR, Disp);
            ylabel('u(x) [m]');
            xlabel('x [m]') ;
            title(' ');
            hold off;
            pause(time)
        end
        
    end
    
end

