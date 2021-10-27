classdef AssigmentPART2 < handle
    
    properties  (Access = private)
        displacements
        Uexact
        nElems
        error
    end
    
    methods (Access = public)
        
        function obj = AssigmentPART2(data) % This constructor requires the same data as class Part6.
            obj.nElems = obj.computeNElems();
            obj.displacements = obj.computeDisplacments(data);
            obj.error = obj.computeError();
        end
        
        function plot(obj, time) % Creates a plot for the exact solution
            obj.createPlot(time, 'log( ErrorU )', 1); % creates the plot for U
            obj.createPlot(time, 'log( ErrorDU )', 2); % creates the plot for dU
        end
        
    end
    
    methods (Access = private)
        
        function displacements = computeDisplacments(obj, data) %computes the displacments for each Number of element
            nElems = obj.nElems;
            a = length(nElems);
            for i = 1:1:a
                b = nElems(i);
                displacements(i) = Part6(b, data);
            end
        end
        
        function Error = computeError(obj) %computes the error for 'U' and 'dU' in a matrix of Error
            nElems = obj.nElems;
            a = length(nElems);
            Error = zeros(a-1, 2); % Matrix of errors
            for i = 1:1:a-1
                [Error(i, 1) Error(i, 2)]  = obj.computeErrorUdU(i);
            end
        end
        
        function [ErrorU ErrorDU] = computeErrorUdU(obj, currentNelem) % computes the error for 'U' and 'DU'
            nElems = obj.nElems;
            l = length(nElems);
            coordApprox = obj.displacements(currentNelem).COOR;
            coordExact = obj.displacements(l).COOR;
            dApprox = obj.displacements(currentNelem).displacement;
            dExact = obj.displacements(l).displacement;
            ErrorU = 0;
            ErrorDU = 0;
            for i = 1:1:length(nElems(currentNelem)) % for each element of the beam:
                coordx1 = coordApprox(i);
                coordx2 = coordApprox(i+1);
                found = 0;
                k = 1;
                while found == 0 %finds the second coordinate of the exact solution that is equivalent to the approximated solution
                    if coordx2 == coordExact(k)
                        found = 1;
                    else
                        k = k+1;
                    end
                end
                j = 1;
                found = 0;
                while found == 0 %finds the first coordinate of the exact solution that is equivalent to the approximated solution
                    if coordx1 == coordExact(j)
                        found = 1;
                    else
                        j = j+1;
                    end
                end
                displacementN1 = dApprox(i);
                displacementN2 = dApprox(i+1);
                displacementN1Exact = dExact(j);
                displacementN2Exact = dExact(k);
                ErrorU = ErrorU + obj.GaussQuadratureU(displacementN1, displacementN2, displacementN1Exact, displacementN2Exact,  coordx1, coordx2);
                ErrorDU = ErrorDU + obj.GaussQuadratureDU(displacementN1, displacementN2, displacementN1Exact, displacementN2Exact,  coordx1, coordx2);
            end
            %See demostration of these equations in the report in order to understand the process
            ErrorU = sqrt(ErrorU);
            ErrorDU = sqrt(ErrorDU);
        end
        
        function sizeElement = computeSizeElement(obj) %computes the size of the elements for each number of element
            nElems = obj.nElems;
            a = length(nElems);
            nElems = nElems(1:a-1);
            sizeElement = zeros(a-1, 1);
            for i = 1:1:a-1
                sizeElement(i) = 1/nElems(i); % 1 is the length of the beam
            end
        end
        
        function createPlot(obj, time, xLabel, typeError) % Creates the plot and its configuration
            sizeElement = obj.computeSizeElement();
            Error = obj.error(:,typeError); % typeError:  1 ---> error u  or  2 ---->  error du
            close all;
            figure;
            hold on;
            plot(log(sizeElement), log(Error), 'b--o');
            xlabel('log( Element size )');
            ylabel(xLabel) ;
            title(' ');
            hold off;
            pause(time)
        end
        
    end
    
    methods (Access = private, Static)
        
        function nElems = computeNElems() % computes the number of elements that will be used
            nElems = [8 16 32 64 128 256 512 32768]; %last element is considered as exact solution
        end
        
        function Error = GaussQuadratureU(displacement1, displacement2, displacement1Exact, displacement2Exact,  coordx1, coordx2) %Computes the error using the GaussQuadrature for U
            xiG(1) = sqrt(3/5);
            xiG(2) = -sqrt(3/5);
            xiG(3) = 0;
            w(1) = 5/9;
            w(2) = 5/9;
            w(3) = 8/9;
            qFun1 = ((1-xiG)*(displacement1Exact - displacement1) + (1+xiG)*(displacement2Exact - displacement2)).^2; %See demostration of this equation in the report
            int = w*qFun1';
            he = coordx2 - coordx1;
            Error = 1/8*he*int;
        end
        
        function Error = GaussQuadratureDU(displacement1, displacement2, displacement1Exact, displacement2Exact,  coordx1, coordx2) %Computes the error using the GaussQuadrature for dU
            qFun1 = ((-displacement1Exact + displacement2Exact) - (-displacement1 + displacement2))^2;  %See demostration of this equation in the report
            he = coordx2 - coordx1;
            Error = 1/he*qFun1;
        end
        
    end
end

