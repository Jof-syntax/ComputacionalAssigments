classdef AssigmentPART2 < handle
    
    properties  (Access = private)
        displacements
        Uexact
        Nelems
        error
    end
    
    methods (Access = public)
        
        function obj = AssigmentPART2(data)
            obj.Nelems = obj.computeNElems();
            obj.displacements = obj.computeDisplacments(data);
            obj.error = obj.computeError();
        end
                
        function plot(obj, time) % Creates a plot for the exact solution
            obj.createPlot(time, 'ErrorU', 1); % creates the plot for U
            obj.createPlot(time, 'ErrorDU', 2); % creates the plot for dU
        end
        
    end
    
    methods (Access = private)
        
        function Nelems = computeNElems(obj)
            Nelems = [4 8 16 32 64 2048]; %last element is considered as exact
        end
        
        function displacements = computeDisplacments(obj, data)
            Nelems = obj.Nelems;
            a = length(Nelems);
            for i = 1:1:a
                b = Nelems(i);
                displacements(i) = Part6(b, data);
            end
        end
        
        function Error = computeError(obj)
            Nelems = obj.Nelems;
            a = length(Nelems);
            Error = zeros(a-1, 2); % Matrix of errors
            for i = 1:1:a-1
                [Error(i, 1) Error(i, 2)]  = obj.computeErrorUdU(i);
            end
        end
        
        function [ErrorU ErrorDU] = computeErrorUdU(obj, currentNelem)
            Nelems = obj.Nelems;
            a = length(Nelems);
            coord_approx = obj.displacements(currentNelem).COOR;
            coord_exact = obj.displacements(a).COOR;
            d_approx = obj.displacements(currentNelem).displacement;
            d_exact = obj.displacements(a).displacement;
            Error = zeros(1,2);
            Error1 = 0;
            Error2 = 0;
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
                Error1 = Error1 + obj.GaussquadratureU(displacement1, displacement2, displacement1Exact, displacement2Exact,  coordx1, coordx2);
                Error2 = Error2 + obj.GaussquadratureDU(displacement1, displacement2, displacement1Exact, displacement2Exact,  coordx1, coordx2);
            end
            ErrorU = sqrt(Error1);
            ErrorDU = sqrt(Error2);
        end
        
        function sizeElement = computeSizeElement(obj)
            Nelems = obj.Nelems;
            a = length(Nelems);
            Nelems = Nelems(1:a-1);
            sizeElement = zeros(a-1);
            for i = 1:1:a-1
            sizeElement(i) = 1/Nelems(i);
            end
        end
        
        function createPlot(obj, time, xLabel, typeError) % Creates the plot and its configuration
            sizeElement = obj.computeSizeElement();
            Error = obj.error(:,typeError); % 1 ---> error u    2 ---->  error du
            close all;
            figure;
            hold on;            
            loglog(sizeElement, Error);
            xlabel('Element size');
            ylabel(xLabel) ;
            title(' ');
            hold off;
            pause(time)
        end
        
    end
    
    methods (Access = private, Static)
        
        function Error = GaussquadratureU(displacement1, displacement2, displacement1Exact, displacement2Exact,  coordx1, coordx2)
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
        
        function Error = GaussquadratureDU(displacement1, displacement2, displacement1Exact, displacement2Exact,  coordx1, coordx2)
            qFun1 = ((-displacement1Exact + displacement2Exact) - (-displacement1 + displacement2))^2;
            he = coordx2 - coordx1;
            Error = 1/he*qFun1;
        end
        
    end
end

