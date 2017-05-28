classdef Line2d < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        v
        p0
    end
    methods (Static)
        function p = intersect( l1, l2)
            p = Point2d(0,0);
            ilambda = (l1.v.y * (l1.p0.x - l2.p0.x) - l1.v.x * (l1.p0.y - l2.p0.y)) / (l1.v.y * l2.v.x - l1.v.x * l2.v.y);
            p.x = l2.p0.x + ilambda * l2.v.x;
            p.y = l2.p0.y + ilambda * l2.v.y;
        end
    end
    methods
        function obj = Line2d(v, p0)
            obj.v.x = v.x;
            obj.v.y = v.y;
            obj.p0.x = p0.x;
            obj.p0.y = p0.y;
        end
        
        function p = getPoint(obj, lambda)
            p = Point2d(0,0);
            p.x = lambda * obj.v.x + obj.p0.x;
            p.y = lambda * obj.v.y + obj.p0.y;
        end
        
        function p = setPoint(obj, p, lambda)
            p.x = lambda * obj.v.x + obj.p0.x;
            p.y = lambda * obj.v.y + obj.p0.y;
        end
    end
    
end

