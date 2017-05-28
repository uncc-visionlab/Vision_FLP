classdef LineSegment2d < Line2d
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        l_start
        l_end
    end
    
    methods
        function obj = LineSegment2d(v, p0, l_start, l_end)
            obj@Line2d(v, p0);
            obj.l_start = l_start;
            obj.l_end = l_end;
        end
    end
    
end

