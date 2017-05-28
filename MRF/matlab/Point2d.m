classdef Point2d < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x
        y 
    end
    
    methods
        function obj = Point2d(x,y,z)
            obj.x = x;
            obj.y = y;
        end       
    end
    
end

