classdef Plane3d < Point3d
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        d
        ANGLE_THRESHOLD
        COPLANAR_COS_ANGLE_THRESHOLD
        PERPENDICULAR_SIN_ANGLE_THRESHOLD
        COPLANAR_COS_EXPECTED_DIHEDRAL_ANGLE
    end
    methods (Static)
        function test
            myplane3d = Plane3d(1,1,1,1);
        end
    end
    
    methods
        function obj = Plane3d(x,y,z,d)
            obj@Point3d(x,y,z);
            obj.ANGLE_THRESHOLD = 5;
            obj.COPLANAR_COS_ANGLE_THRESHOLD = cos(obj.ANGLE_THRESHOLD * pi / 180.); % degrees
            obj.PERPENDICULAR_SIN_ANGLE_THRESHOLD = sin(obj.ANGLE_THRESHOLD * pi / 180.); % degrees
            obj.COPLANAR_COS_EXPECTED_DIHEDRAL_ANGLE = 1.; % cos(0)
            obj.d = d;
        end
        
        function val = orthogonalDistanceSquared(obj, pt)
            val = obj.evaluate(pt) * obj.evaluate(pt);
        end
        
        function val = signedOrthogonalDistance(obj, pt)
            val = obj.evaluate(pt);
        end
        
        function [line, boolval] = intersect(planeA, planeB)
            line = Line3d(Point3d(0,0,0), Point3d(0,0,0));
            line.v = planeA.cross(planeB);
            %double detB = (planeA.x * planeB.y - planeB.x * planeA.y);
            %double detA = line.v.z;
            if (line.v.z == 0)
                boolval = false;
                return;
            end
            %std::cout << "detA " << detA << " detB " << detB << std::endl;
            line.p0.x = (-planeB.y * planeA.d + planeA.y * planeB.d) / line.v.z;
            line.p0.y = (planeB.x * planeA.d - planeA.x * planeB.d) / line.v.z;
            line.p0.z = 0;
            boolval = true;
        end
        
        function setCoeffs( obj, x,  y, z, d)
            obj.x = x;
            obj.y = y;
            obj.z = z;
            obj.d = d;
        end
        
        function scale(obj, scalef)
            obj.x = obj.x * scalef;
            obj.y = obj.y * scalef;
            obj.z = obj.z * scalef;
            obj.d = obj.d * scalef;
        end
        
        function val =  evaluate(obj, pt)
            val =  obj.x * pt.x + obj.y * pt.y + obj.z * pt.z + obj.d;
        end
        
        function setCoeffsFromPoints(obj, pt1, pt2, pt3)
            obj.x = (pt2.y - pt1.y)*(pt3.z - pt1.z)-(pt3.y - pt1.y)*(pt2.z - pt1.z);
            obj.y = (pt2.z - pt1.z)*(pt3.x - pt1.x)-(pt3.z - pt1.z)*(pt2.x - pt1.x);
            obj.z = (pt2.x - pt1.x)*(pt3.y - pt1.y)-(pt3.x - pt1.x)*(pt2.y - pt1.y);
            obj.d = -(obj.x * pt1.x + obj.y * pt1.y + obj.z * pt1.z);
        end
        
        function val = cosDihedralAngle(obj, test_plane)
            val = obj.x*test_plane.x + obj.y*test_plane.y + obj.z*test_plane.z;
        end
        
        function val = angleDistance(obj, planeA)
            val = obj.COPLANAR_COS_EXPECTED_DIHEDRAL_ANGLE - obj.cosDihedralAngle(planeA);
        end
        
        function boolval = epsilonEquals(obj, planeA, epsilon)
            if ~exist('epsilon', 'var')
                epsilon = obj.COPLANAR_COS_ANGLE_THRESHOLD;
            end
            
            boolval = obj.COPLANAR_COS_EXPECTED_DIHEDRAL_ANGLE - obj.cosDihedralAngle(planeA) ...
                < epsilon;
            
        end
        
        function boolval = epsilonPerpendicular(obj, planeA, epsilon)
            if ~exist('epsilon', 'var')
                epsilon = obj.PERPENDICULAR_SIN_ANGLE_THRESHOLD;
            end
            
            boolval = abs(obj.cosDihedralAngle(planeA)) < epsilon;
            
        end
        
        function interpolate(obj, alpha, planeA, planeB, pt)
            
            obj.x = alpha*planeA.x + (1 - alpha)*planeB.x;
            obj.y = alpha*planeA.x + (1 - alpha)*planeB.x;
            obj.z = alpha*planeA.x + (1 - alpha)*planeB.x;
            obj.d = -(obj.x * pt.x + obj.y * pt.y + obj.z * pt.z);
            
%             obj.x = alpha*planeA.x + (1 - alpha)*planeB.x;
%             obj.y = alpha*planeA.y + (1 - alpha)*planeB.y;
%             obj.z = alpha*planeA.z + (1 - alpha)*planeB.z;
%             obj.d = -(obj.x * pt.x + obj.y * pt.y + obj.z * pt.z);
            
        end
            
        function convertHessianNormalForm(obj) 
            normScale = 1.0 / sqrt(obj.x * obj.x + ...
                    obj.y * obj.y + obj.z * obj.z);
            obj.scale(normScale);
        end
    end
end

