/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   depth.h
 * Author: arwillis
 *
 * Created on February 28, 2017, 4:22 AM
 */

#ifndef DEPTH_H
#define DEPTH_H

#include <boost/shared_ptr.hpp>
#include <sstream>
#include <math.h>
//#ifdef __cplusplus
//extern "C" {
//#endif

    template<typename _Tpl> class Point3_ {
    public:
        Point3_() : x(0), y(0), z(0) {
        }

        Point3_(_Tpl _x, _Tpl _y, _Tpl _z) : x(_x), y(_y), z(_z) {
        }

        Point3_<_Tpl> clone() {
            return Point3_<_Tpl>(this->x, this->y, this->z);
        }
        
        friend std::ostream& operator<<(std::ostream& os, const Point3_<_Tpl>& p) {
            os << "[" << p.x << ", " << p.y << ", " << p.z << "]";
            return os;
        }        
        _Tpl x, y, z;
    };
    typedef Point3_<float> Point3f;
    typedef Point3_<double> Point3d;
    typedef Point3_<int> Point3i;

    template<typename _Tpl> class Plane3_ : public Point3_<_Tpl> {
    public:

//        static constexpr float ANGLE_THRESHOLD = 5.f;
//        static constexpr float COPLANAR_COS_ANGLE_THRESHOLD = cos(ANGLE_THRESHOLD * M_PI / 180.f); // degrees
//        static constexpr float PERPENDICULAR_SIN_ANGLE_THRESHOLD = sin(ANGLE_THRESHOLD * M_PI / 180.f); // degrees
//        static constexpr float COPLANAR_COS_EXPECTED_DIHEDRAL_ANGLE = 1.f; // cos(0)

        typedef boost::shared_ptr<Plane3_<_Tpl> > Ptr;
        typedef boost::shared_ptr<const Plane3_<_Tpl> > ConstPtr;

        Plane3_() : d(0), Point3_<_Tpl>() {
        }

        Plane3_(_Tpl _x, _Tpl _y, _Tpl _z, _Tpl _d) : d(_d), Point3_<_Tpl>(_x, _y, _z) {
        }

        Plane3_<_Tpl> clone() {
            return Plane3_<_Tpl>(this->x, this->y, this->z, this->d);
        }

        void alignPlanes(std::vector<Plane3_<_Tpl> > planesA, std::vector<Plane3_<_Tpl> > planesB) {

        }

        _Tpl orthogonalDistanceSquared(Point3_<_Tpl> pt) {
            return evaluate(pt) * evaluate(pt);
        }

        _Tpl signedOrthogonalDistance(Point3_<_Tpl> pt) {
            return evaluate(pt);
        }

//        bool intersect(Plane3_<_Tpl> planeB, Line3_<_Tpl>& line) {
//            line.v = this->cross(planeB);
//            //double detB = (planeA.x * planeB.y - planeB.x * planeA.y);
//            //double detA = line.v.z;
//            if (line.v.z == 0) {
//                return false;
//            }
//            //std::cout << "detA " << detA << " detB " << detB << std::endl;
//            line.p0.x = (_Tpl) (-planeB.y * this->d + this->y * planeB.d) / line.v.z;
//            line.p0.y = (_Tpl) (planeB.x * this->d - this->x * planeB.d) / line.v.z;
//            line.p0.z = 0;
//            return true;
//        }

        void setCoeffs(_Tpl _x, _Tpl _y, _Tpl _z, _Tpl _d) {
            this->x = _x;
            this->y = _y;
            this->z = _z;
            this->d = _d;
        }

        void scale(_Tpl scalef) {
            this->x *= scalef;
            this->y *= scalef;
            this->z *= scalef;
            this->d *= scalef;
        }

        _Tpl evaluate(const Point3_<_Tpl>& pt) {
            return this->x * pt.x + this->y * pt.y + this->z * pt.z + this->d;
        }

        void setCoeffs(const Point3_<_Tpl>& pt1, const Point3_<_Tpl>& pt2,
                const Point3_<_Tpl>& pt3) {
            this->x = (pt2.y - pt1.y)*(pt3.z - pt1.z)-(pt3.y - pt1.y)*(pt2.z - pt1.z);
            this->y = (pt2.z - pt1.z)*(pt3.x - pt1.x)-(pt3.z - pt1.z)*(pt2.x - pt1.x);
            this->z = (pt2.x - pt1.x)*(pt3.y - pt1.y)-(pt3.x - pt1.x)*(pt2.y - pt1.y);
            this->d = -(this->x * pt1.x + this->y * pt1.y + this->z * pt1.z);
        }

        _Tpl cosDihedralAngle(const Plane3_<_Tpl> test_plane) const {
            return this->x * test_plane.x + this->y * test_plane.y + this->z * test_plane.z;
        }

//        _Tpl angleDistance(Plane3_<_Tpl> planeA) const {
//            return (_Tpl) COPLANAR_COS_EXPECTED_DIHEDRAL_ANGLE - cosDihedralAngle(planeA);
//        }
//
//        bool epsilonEquals(Plane3_<_Tpl> planeA, float eps = COPLANAR_COS_ANGLE_THRESHOLD) const {
//            return (COPLANAR_COS_EXPECTED_DIHEDRAL_ANGLE - cosDihedralAngle(planeA)
//                    < COPLANAR_COS_ANGLE_THRESHOLD);
//        }

//        bool epsilonPerpendicular(Plane3_<_Tpl> planeA, float eps = PERPENDICULAR_SIN_ANGLE_THRESHOLD) const {
//            return (abs(cosDihedralAngle(planeA)) < eps);
//        }

        void interpolate(float alpha, Plane3_<_Tpl> planeA, Plane3_<_Tpl> planeB,
                Point3_<_Tpl> pt) {
            this->x = alpha * planeA.x + (1 - alpha) * planeB.x;
            this->y = alpha * planeA.x + (1 - alpha) * planeB.x;
            this->z = alpha * planeA.x + (1 - alpha) * planeB.x;
            this->d = -(this->x * pt.x + this->y * pt.y + this->z * pt.z);
        }

        void convertHessianNormalForm() {
            float normScale = 1.f / sqrt(this->x * this->x +
                    this->y * this->y + this->z * this->z);
            scale(normScale);
        }

        friend std::ostream& operator<<(std::ostream& os, const Plane3_<_Tpl>& p) {
            os << "[" << p.x << ", " << p.y << ", " << p.z
                    << ", " << p.d << "]";
            return os;
        }

        std::string toString() {
            std::ostringstream stringStream;
            stringStream << "(" << this->x << ", " << this->y
                    << ", " << this->z << ", " << this->d << ")";
            return stringStream.str();
        }

        _Tpl d;
    };
    typedef Plane3_<float> Plane3f;
    typedef Plane3_<double> Plane3d;
    typedef Plane3_<int> Plane3i;



//#ifdef __cplusplus
//}
//#endif

#endif /* DEPTH_H */

