classdef Point3d < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x
        y 
        z
    end
    
    methods
        function obj = Point3d(x,y,z)
            obj.x = x;
            obj.y = y;
            obj.z = z;
        end       
        
        function crossval = cross(p1, p2)
            crossval = Point3d(0,0,0);
            crossval.x = p1.y*p2.z - p1.z*p2.y;
            crossval.y = p1.z*p2.x - p1.x*p2.z;
            crossval.z = p1.x*p2.y - p1.y*p2.x;
        end
    end
    
end

