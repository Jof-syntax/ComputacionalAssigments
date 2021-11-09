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
            Error = zeros(a, 2); % Matrix of errors
            for i = 1:1:a
                [Error(i, 1) Error(i, 2)]  = obj.computeErrorUdU(i);
            end
        end
        
        function [ErrorU ErrorDU] = computeErrorUdU(obj, currentNelem) % computes the error for 'U' and 'DU'
            nElems = obj.nElems;
            coordApprox = obj.displacements(currentNelem).COOR;
            dApprox = obj.displacements(currentNelem).displacement;
            ErrorU = 0;
            ErrorDU = 0;
            for i = 1:1:nElems(currentNelem) % for each element of the beam:
                coordx1 = coordApprox(i);
                coordx2 = coordApprox(i+1);
               
                displacementN1 = dApprox(i); %Displacement for approx solution
                displacementN2 = dApprox(i+1);
                
                ErrorU = ErrorU + obj.GaussQuadratureU(displacementN1, displacementN2,  coordx1, coordx2);
                ErrorDU = ErrorDU + obj.GaussQuadratureDU(displacementN1, displacementN2,  coordx1, coordx2);
            end
            %See demostration of these equations in the report in order to understand the process
            ErrorU = sqrt(ErrorU);
            ErrorDU = sqrt(ErrorDU);
        end
        
        function sizeElement = computeSizeElement(obj) %computes the size of the elements for each number of element
            nElems = obj.nElems;
            a = length(nElems);
            sizeElement = zeros(a, 1);
            for i = 1:1:a
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
            nElems = [8 16 32 64 128 256 512];
        end
        
        function Error = GaussQuadratureU(displacement1, displacement2, coordx1, coordx2) %Computes the error using the GaussQuadrature for U
            xiG(1) = 0.339981043584856;
            xiG(2) = -0.339981043584856;
            xiG(3) = 0.861136311594053;
            xiG(4) = -0.861136311594053;
            w(1) = 0.652145154862546;
            w(2) = 0.652145154862546;
            w(3) = 0.347854845137454;
            w(4) = 0.347854845137454;
            he = coordx2 - coordx1;
            uExactFunction = @(x) 0.01*cos(pi*x) + pi/100*sin(pi*x)+0.098696*x^2-0.02; %Exact solution
            Error = 0;
            for i = 1:length(w)
                N = 1/2*[1-xiG(i) 1+xiG(i)];
                UExact = uExactFunction(N*[coordx1; coordx2]);
                qFun1 = (UExact - 1/2*((1-xiG(i))*displacement1 + (1+xiG(i))*displacement2))^2; %See demostration of this equation in the report
                Error = Error + w(i)*qFun1;
            end
            Error = he/2*Error;
        end
        
        function Error = GaussQuadratureDU(displacement1, displacement2,  coordx1, coordx2) %Computes the error using the GaussQuadrature for dU
            xiG(1) = 0.339981043584856;
            xiG(2) = -0.339981043584856;
            xiG(3) = 0.861136311594053;
            xiG(4) = -0.861136311594053;
            w(1) = 0.652145154862546;
            w(2) = 0.652145154862546;
            w(3) = 0.347854845137454;
            w(4) = 0.347854845137454;
            DuExactFunction = @(x) 0.19739*x+0.098696*cos(pi*x)-pi/100*sin(pi*x); % Exact solution '
            he = coordx2 - coordx1;
            Error = 0;
            for i = 1:length(w)
                N = 1/2*[1-xiG(i) 1+xiG(i)];
                DUExact = DuExactFunction(N*[coordx1; coordx2]);
                qFun1 = (DUExact + (displacement1 - displacement2)/he)^2;  %See demostration of this equation in the report
                Error = Error + he/2*qFun1*w(i);
            end
        end
        
    end
end

