classdef Line3d < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        v
        p0
    end
    
    methods
        function obj = Line3d(v, p0)
            obj.v = Point3d(v.x,v.y,v.z);
            obj.p0 = Point3d(p0.x,p0.y,p0.z);
        end
        
        function p = getPoint(obj, lambda)
            p = Point3d(0,0,0);
            p.x = lambda * obj.v.x + obj.p0.x;
            p.y = lambda * obj.v.y + obj.p0.y;
            p.z = lambda * obj.v.z + obj.p0.z;
        end
        
        function p = setPoint(obj, p, lambda)
            p.x = lambda * obj.v.x + obj.p0.x;
            p.y = lambda * obj.v.y + obj.p0.y;
            p.z = lambda * obj.v.z + obj.p0.z;
        end

        function dist = distanceSquared(obj, pt) 
            num = obj.v.cross(Point3d(pt.x - obj.p0.x, ...
                pt.y - obj.p0.y, pt.z - obj.p0.z));
            dist = (num.x * num.x + num.y * num.y + num.z * num.z) / ...
                (obj.v.x * obj.v.x + obj.v.y * obj.v.y + obj.v.z * obj.v.z);
        end

    end
    
end

