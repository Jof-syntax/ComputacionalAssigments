%% ASSIGMENT 1 by Jofre Geli Cerdan

%NOTE: By default, each plot will appear for 7 seconds.

classdef main < handle
    
    properties (Access = private)
        time
        Nelem
    end
    
    methods (Access = public)
        
        function obj = main()
            %%%%%%%%%%%%%%%%%%%% Part1 %%%%%%%%%%%%%%%%%%%%
            clc;
            close all;
            obj.createTimeForPlot();
            obj.computeNelem();
            time = obj.time;
            %a = Part2();
            %a.plot(time);
            %b = Part5(a.solution);
%             b.plot1(time);
%             b.plot2(time);
%             b.plot3(time);
%             b.plot4(time);
            c = Part6(5);
            c.plot();
%             hold on;
%             c = Part6(10);
%             c.plot();
%             c = Part6(20);
%             c.plot();
%             c = Part6(40);
%             c.plot();
%             legend('Nelem = 5','Nelem = 10', 'Nelem = 20', 'Nelem = 40', 'Fontsize', 14)
%fer class part 7 que et grafiqui els 4 casos
              %%%%%%%%%%%%%%%%%%%% Part2 %%%%%%%%%%%%%%%%%%%%
        end
    end
    
    methods (Access = private)
        
        function computeNelem(obj)
            obj.Nelem = 40;
        end
        
        function createTimeForPlot(obj)
            obj.time = 7;
        end
        
    end
    
end
