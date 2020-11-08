
#ifndef Node_hpp
#define Node_hpp
#endif /* Node_hpp */

#include "Point.hpp"

class Node {
protected:

public:
    Point* point;
    Node* parent;
    Node* children[8];
    double xmin, xmax, ymin, ymax, zmin, zmax, s;
    bool is_aggregate;
    
    ~Node();
    Node();
    Node(Node* parent0, double xmin0, double ymin0, double zmin0, double xmax0, double ymax0, double zmax0);
    bool contains(Point* p);
    void split(Point* new_point);
    void insert_point(Point* p);
};
