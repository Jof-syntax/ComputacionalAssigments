classdef Part2 < handle
    
    properties (Access = public)
        solution %uExact
    end
     
    methods (Access = public)
        
        function obj = Part2()
            obj.solution = obj.computeUExact();
        end
        
        function plot(obj, time)
            obj.createPlot(time);
        end
        
    end
    
    methods (Access = private)
        
        function uExact = computeUExact(obj)
            syms x u(x)
            L = 1;
            g = 0.01;
            rho = pi^2/L^2;
            s = g*rho^2;
            b = g*pi^2/L;
            q = (rho*u(x)-s*x^2);
            Du = diff(u,x);
            uExact = dsolve(diff(u,2) == -q, u(0)== -g, Du(1)== b);
            uExact = simplify(uExact);
            obj.representSolution(uExact);
        end
                
        function createPlot(obj, time)
            close all;
            figure;
            hold on;
            uExact = obj.solution;
            ezplot(uExact, 0, 1);
            ylabel('u(x) [m]'); 
            xlabel('x [m]') ; 
            title(' ');
            hold off;
            pause(time)
        end
        
    end
    
    methods (Access = private, Static)
        
        function representSolution(show)
            disp('----------------------------------------------------Part 2----------------------------------------------------');
            syms x
            disp('Pretty form :');
            pretty(show);
            disp('Latex form :');
            latex(show)
        end
        
    end
end

