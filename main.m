%% ASSIGMENT 1 by Jofre Geli Cerdan

%NOTE: By default, each plot will appear for 7 seconds.

classdef main < handle
    
    properties (Access = private)
        time
    end
    
    methods (Access = public)
        
        function obj = main()
            clc;
            close all;
            obj.createTimeForPlot();
            time = obj.time;
            a = Part2();
            a.plot(time);
            b = Part5(a.solution);
            b.plot1(time);
            b.plot2(time);
            b.plot3(time);
            b.plot4(time);
        end
    end
    
    methods (Access = private)
        
        function createTimeForPlot(obj)
            obj.time = 7;
        end
        
    end
    
end
