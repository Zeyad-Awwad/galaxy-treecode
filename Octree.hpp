
#ifndef Octree_hpp
#define Octree_hpp
#endif /* Octree_hpp */

#include "Node.hpp"
#include <vector>
#include "math.h"

class Octree {
protected:
    
public:
    Node* root;
    double max_opening;
    double eps;
    double dt;
    
    
    Octree(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double opening, double epsilon, double timestep);
    void clear( vector<Point> &points );
    bool add_point(Point* p);
    void add_point_vector( vector<Point> &points );
    void calculate_aggregates();
    void aggregate(Node* n, double &m0, double &x0, double &y0, double &z0);
    void calculate_all_accelerations();
    void calculate_accelerations(Node* n1, Node* n2);
    void update_accelerations(Node* n1, Node* n2);
    void initialize_leapfrog(Node* n);
    void update_positions(Node* n);
};
