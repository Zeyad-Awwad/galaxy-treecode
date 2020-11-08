
#ifndef Point_hpp
#define Point_hpp
#endif /* Point_hpp */

#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

class Point {
protected:
    
public:
    double m, x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    
    Point();
    Point(double m0, double x0, double y0, double z0, double vx0, double vy0, double vz0);
    
};

