classdef Part2 < handle
    
    properties (Access = public)
        displacement
    end
    
    properties (Access = private)
        g
        L
    end
     
    methods (Access = public)
        
        function obj = Part2(data)
            obj.L = data.L;
            obj.g = data.g;
            obj.displacement = obj.computeUExact();
        end
        
        function representSolution(obj)
             obj.computeRepresentSolution(obj.displacement);
        end
        
        function plot(obj, time)
            obj.createPlot(time);
        end
        
    end
    
    methods (Access = private)
        
        function uExact = computeUExact(obj)
            syms x u(x)
            L = obj.L;;
            g = obj.g;
            rho = pi^2/L^2;
            s = g*rho^2;
            b = g*pi^2/L;
            q = (rho*u(x)-s*x^2);
            Du = diff(u,x);
            uExact = dsolve(diff(u,2) == -q, u(0)== -g, Du(1)== b);
            uExact = simplify(uExact);
        end
                
        function createPlot(obj, time)
            close all;
            figure;
            hold on;
            uExact = obj.displacement;
            ezplot(uExact, 0, 1);
            ylabel('u(x) [m]'); 
            xlabel('x [m]') ; 
            title(' ');
            hold off;
            pause(time)
        end
        
    end
    
    methods (Access = private, Static)
        
        function computeRepresentSolution(show)
            disp('----------------------------------------------------Part 2----------------------------------------------------');
            syms x
            disp('Pretty form :');
            pretty(show);
            disp('Latex form :');
            latex(show)
        end
        
    end
end

