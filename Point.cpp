
#include "Point.hpp"

Point::Point() {
    m = 0; x = 0; y = 0; z = 0;
    vx = 0; vy = 0; vz = 0;
    ax = 0; ay = 0; az = 0;
}

Point::Point(double m0, double x0, double y0, double z0, double vx0, double vy0, double vz0) {
    m = m0; x = x0; y = y0; z = z0;
    vx = vx0; vy = vy0; vz = vz0;
    ax = 0; ay = 0; az = 0;
}

