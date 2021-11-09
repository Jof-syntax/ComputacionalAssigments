classdef Part5 < handle
    
    properties (Access = public)
        solution %Approximated solution
    end
    
    properties (Access = private)
        uExact
        Ni
        parameters
        g
        L
    end
    
    methods (Access = public)
        
        function obj = Part5(uExact, data) %The constructor needs the values of parameters 'g' and 'L' and the displacemnts calculated in 'Part2' to work. This function obtains the displacements by means of 'computeDisplacements'
            obj.uExact = uExact;
            obj.createParameters(data)
            obj.createNi();
            obj.computeDisplacements();
        end
        
        function plot1(obj, time) % Plots the approximated solution 'Polynom x' and the exact solution for a given time 'time'
            Sol1 = obj.solution.one;
            obj.createPlot(Sol1, time);
        end
        
        function plot2(obj, time) % Plots the approximated solution 'Polynom x^2' and the exact solution for a given time 'time'
            Sol2 = obj.solution.two;
            obj.createPlot(Sol2, time);
        end
        
        function plot3(obj, time) % Plots the approximated solution 'Polynom x^3' and the exact solution for a given time 'time'
            Sol3 = obj.solution.three;
            obj.createPlot(Sol3, time);
        end
        
        function plot4(obj, time) % Plots the approximated solution 'Polynom x^4' and the exact solution for a given time 'time'
            Sol4 = obj.solution.four;
            obj.createPlot(Sol4, time);
        end
        
    end
    
    methods (Access = private)
        
        function computeDisplacements(obj) %Computes the approximated solution for the 4 cases
            syms x
            Ni1 = obj.Ni.Ni1;
            Ni2 = obj.Ni.Ni2;
            Ni3 = obj.Ni.Ni3;
            Ni4 = obj.Ni.Ni4;
            obj.solution.one = obj.computeGalerkinMethod(Ni1);
            obj.solution.two = obj.computeGalerkinMethod(Ni2);
            obj.solution.three = obj.computeGalerkinMethod(Ni3);
            obj.solution.four = obj.computeGalerkinMethod(Ni4);
        end
        
        function createNi(obj) %Creats the 4 polynomials used to calulate the approximated solution
            syms x
            Ni.Ni1 = [1, x];
            Ni.Ni2 = [1, x, x^2];
            Ni.Ni3 = [1, x, x^2, x^3];
            Ni.Ni4 = [1, x, x^2, x^3, x^4];
            obj.Ni = Ni;
        end
        
        function createParameters(obj, data) % Creates all the parameters used to calulate the approximated solution
            syms x
            obj.parameters.L = data.L;
            L = obj.parameters.L;
            obj.parameters.g = data.g;
            g = obj.parameters.g;
            rho = pi^2/L^2;
            obj.parameters.rho = rho;
            s = g*rho^2;
            obj.parameters.s = s;
            obj.parameters.b = g*pi^2/L;
            obj.parameters.f = s*x^2;
        end
        
        function [u] = computeGalerkinMethod(obj, N) %Computes the displacments using the Garlekin method (Approximated solution) for a given polynomial
            syms x
            f = obj.parameters.f;
            b = obj.parameters.b;
            g = obj.parameters.g;
            rho = obj.parameters.rho;
            B = diff(N,x);
            BtB = B.'*B;
            NtN = N.'*N;
            K = int(BtB,0,1) - int(rho*NtN, 0, 1);
            N_1 = subs(N,1);
            F = -int(N.'*f,0,1) + N_1.'*b;
            r = 1;
            l = 2:length(N) ;
            dl = K(l,l)\(F(l) + K(l,r)*g);
            u  = - g + N(l)*dl;
        end
        
        function createPlot(obj, approxSol, time) % Creates the plot and its configuration
            close all;
            figure;
            hold on;
            uExact = obj.uExact;
            h1 = ezplot(uExact, 0, 1);
            h2 = ezplot(approxSol, [0,1]);
            legend([h1 h2],{['Exact: ',char(vpa(uExact,3))],['Polynomial basis: ',char(vpa(approxSol,3))]})
            title(' ');
            xlabel('x(x) [m]') ;
            ylabel('u(x) [m]');
            pause(time)
        end
        
    end
    
end


