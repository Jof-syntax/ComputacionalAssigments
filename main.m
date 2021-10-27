%% ASSIGMENT 1 by Pau Tusell Tresserras & Jofre Geli Cerdan

% NOTE: By default, each plot will appear for 7 seconds. See function 'createTimeForPlot' to change this time.

classdef main < handle
    
    properties (Access = private)
        time
        Nelem
        data
    end
    
    methods (Access = public)
        
        function obj = main()
            obj.init();
            
            %%%%%%%%%%%%%%%%%%%% Part1 %%%%%%%%%%%%%%%%%%%%
            part2 = Part2(obj.data);
            part2.representSolution();
            
            part5 = Part5(part2.displacement, obj.data);
            
            part6 = Part6(obj.Nelem, obj.data);
            
            %%%%%%%%%%%%%%%%%%%% Part2 %%%%%%%%%%%%%%%%%%%%
            PART2 = AssigmentPART2(obj.data);
            
            %%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%
            obj.createPlot(part2, part5, part6, PART2);
        end
    end
    
    methods (Access = private)
        
        function init(obj) % Initializes all the functions needed to use this Class
            clc;
            close all;
            obj.createTimeForPlot();
            obj.createNelem();
            obj.createData();
        end
        
        function createNelem(obj) %Creates the number of elements used in 'Part6'
            obj.Nelem = 5;
        end
        
        function createTimeForPlot(obj) %Creates the time that the plots will appear
            obj.time = 7;
        end
        
        function createData(obj) %Geometric data of the problem
            obj.data.L = 1;
            obj.data.g = 0.01;
        end
        
        function createPlot(obj, part2, part5, part6, PART2) % Generates all the plots of the Assigment
            part2.plot(obj.time);
            part5.plot1(obj.time);
            part5.plot2(obj.time);
            part5.plot3(obj.time);
            part5.plot4(obj.time);
            part6.plot(obj.time);
            PART2.plot(obj.time);
        end
        
    end
    
end
