%% ASSIGMENT 1 by Pau Tusell Tresserras & Jofre Geli Cerdan

%NOTE: By default, each plot will appear for 7 seconds.

classdef main < handle
    
    properties (Access = private)
        time
        Nelem
        data
    end
    
    methods (Access = public)
        
        function obj = main()
            %%%%%%%%%%%%%%%%%%%% Part1 %%%%%%%%%%%%%%%%%%%%
            clc;
            close all;
            obj.createTimeForPlot();
            obj.createNelem();
            obj.createData();
            part2 = Part2(obj.data);
            %part2.representSolution();
            part5 = Part5(part2.displacement, obj.data);
            part6 = Part6(obj.Nelem, obj.data);
            obj.createPlot(part2, part5, part6);
            %%%%%%%%%%%%%%%%%%%% Part2 %%%%%%%%%%%%%%%%%%%%
            %PART2 = AssigmentPART2(obj.data);
        end
    end
    
    methods (Access = private)
        
        function createNelem(obj)
            obj.Nelem = 5;
        end
        
        function createTimeForPlot(obj)
            obj.time = 1;
        end
        
        function createData(obj)
            obj.data.L = 1;
            obj.data.g = 0.01;
        end
        
        function createPlot(obj, part2, part5, part6)
            %part2.plot(obj.time);
            %part5.plot1(obj.time);
            %part5.plot2(obj.time);
            %part5.plot3(obj.time);
            %part5.plot4(obj.time);
            part6.plot(obj.time);
        end
        
    end
end
