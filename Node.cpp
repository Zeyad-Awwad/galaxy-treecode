
#include "Node.hpp"

Node::Node() {
    s = -1;
    point = nullptr;
    parent = nullptr;
    is_aggregate = false;
    xmin = -1; xmax = -1;
    ymin = -1; ymax = -1;
    zmin = -1; zmax = -1;
    for (int i=0; i<8; i++)
        children[i] = nullptr;
}

Node::Node(Node* parent0, double xmin0, double ymin0, double zmin0, double xmax0, double ymax0, double zmax0) {
    point = nullptr;
    parent = parent0;
    is_aggregate = false;
    xmin = xmin0; xmax = xmax0;
    ymin = ymin0; ymax = ymax0;
    zmin = zmin0; zmax = zmax0;
    s = xmax-xmin;
    if (ymax-ymin > s) s = ymax-ymin;
    if (zmax-zmin > s) s = zmax-zmin;
    for (int i=0; i<8; i++)
        children[i] = nullptr;
}

Node::~Node() {
    if (is_aggregate) {
        delete point;
        for (int i=0; i<8; i++) {
            if (children[i] != nullptr)
                delete children[i];
        }
    }
}

bool Node::contains(Point* p) {
    return p->x >= xmin && p->x < xmax && p->y >= ymin && p->y < ymax && p->z >= zmin && p->z < zmax;
}

void Node::split(Point* new_point) {
    is_aggregate = true;
    int i, j, k, index;
    double dx = (xmax-xmin)/2.0;
    double dy = (ymax-ymin)/2.0;
    double dz = (zmax-zmin)/2.0;
    
    for (k=0; k<2; k++)
        for (j=0; j<2; j++)
            for (i=0; i<2; i++)
            {
                index = 4*k + 2*j + i;                
                children[index] = new Node(this, xmin+i*dx, ymin+j*dy, zmin+k*dz, xmin+(i+1)*dx, ymin+(j+1)*dy, zmin+(k+1)*dz );
                children[index]->contains(point);
                if ( children[index]->contains(point) )
                    children[index]->insert_point(point);
                if ( children[index]->contains(new_point) )
                    children[index]->insert_point(new_point);
            }
    point = nullptr;
}

void Node::insert_point(Point* p) {
    if (point == nullptr)
        point = p;
    else
        split(p);
}
