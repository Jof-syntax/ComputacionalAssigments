classdef Part2 < handle
    
    properties (Access = public)
        displacement % Exact solution
    end
    
    properties (Access = private)
        g
        L
    end
     
    methods (Access = public)
         
        function obj = Part2(data) %The constructor needs the values of parameters 'g' and 'L' to work. This function obtains the displacements by means of 'computeUExact'
            obj.L = data.L;
            obj.g = data.g;
            obj.displacement = obj.computeUExact();
        end
        
        function representSolution(obj) % Represents the exact solution in the command window
             obj.computeRepresentSolution(obj.displacement);
        end
        
        function plot(obj, time) % Creates a plot for the exact solution
            obj.createPlot(time);
        end
        
    end
    
    methods (Access = private)
        
        function uExact = computeUExact(obj) % Computes the exact solution using 'symbolic'
            syms x u(x)
            L = obj.L;;
            g = obj.g;
            rho = pi^2/L^2;
            s = g*rho^2;
            b = g*pi^2/L;
            q = (rho*u(x)-s*x^2);
            Du = diff(u,x);
            uExact = dsolve(diff(u,2) == -q, u(0)== -g, Du(1)== b);
            uExact = vpa(uExact, 4);
        end
                
        function createPlot(obj, time) % Creates the plot and its configuration
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
        
        function computeRepresentSolution(show) % Computes the representation of the exact solution
            disp('----------------------------------------------------Part 2----------------------------------------------------');
            syms x
            disp('Pretty form :');
            pretty(show);
            disp('Latex form :');
            latex(show)
        end
        
    end
end

