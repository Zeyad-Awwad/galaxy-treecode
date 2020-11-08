
#include "Octree.hpp"

Octree::Octree(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, double opening, double epsilon, double timestep) {
    //Defines an octree data structure for cell-to-cell interactions 
    root = new Node(nullptr, xmin, ymin, zmin, zmax, ymax, zmax);
    max_opening = opening;
    eps = epsilon;
    dt = timestep;
}

void Octree::clear( vector<Point> &points ) {
    //Properly deletes all nodes and resets all accelerations
    Node* new_root = new Node(nullptr, root->xmin, root->ymin, root->zmin, root->zmax, root->ymax, root->zmax);
    delete root;
    root = new_root;
    for (int i=0; i<points.size(); i++) {
        points[i].ax = 0;
        points[i].ay = 0;
        points[i].az = 0;
    }
}

bool Octree::add_point(Point* p) {
    //Inserts a new object into a tree (as a node or leaf)
    if ( !(root->contains(p)) )
        return false;
    Node* current = root;
    while (current->is_aggregate)
        for (int i=0; i<8; i++)
            if ( current->children[i]->contains(p) ) {
                current = current->children[i];
                break;
            }
    current->insert_point(p);
    return true;
}

void Octree::add_point_vector( vector<Point> &points ) {
    //Inserts a vector of points into a tree
    for (int i=0; i<points.size(); i++) {
        add_point( &points[i] );
    }
}

void Octree::calculate_aggregates() {
    //Begins the aggregation of octree node properties from the root node
    double m = 0, x = 0, y = 0, z = 0;
    aggregate(root, m, x, y, z);
    root->point = new Point(m, x/m, y/m, z/m, 0, 0, 0);
}

void Octree::aggregate(Node* n, double &m0, double &x0, double &y0, double &z0) {
    //Aggregates the properties of all children for each octree node
    double m = 0, x = 0, y = 0, z = 0;
    if ( !(n->is_aggregate) && n->point != nullptr ) {
        m = n->point->m;     
        x = n->point->x;
        y = n->point->y;     
        z = n->point->z;
    }
    else if (n->is_aggregate) {
        for (int i=0; i<8; i++)
            aggregate(n->children[i], m, x, y, z);
        if ( m > 0 ) {
            x = x/m;
            y = y/m;
            z = z/m;
        }
        n->point = new Point(m, x, y, z, 0, 0, 0);
    }
    m0 += m;        
    x0 += m*x;
    y0 += m*y;        
    z0 += m*z;
}


void Octree::calculate_all_accelerations() {
    //Begins the cell-to-cell acceleration calculations from the children of the root node
    for (int i=0; i<8; i++) {
        for (int j=0; j<8; j++) {
            calculate_accelerations(root->children[i], root->children[j]);
        }
    }
}


void Octree::calculate_accelerations(Node* n1, Node* n2) {
    //Uses cell-to-cell interactions to estimate gravitational acceleration
    
    if ( n1->point == nullptr || n2->point == nullptr ) return;
    if ( n1 == n2 && !( n1->is_aggregate ) ) return;

    double d, dx, dy, dz, opening;


    if ( n1 == n2 )
        opening = max_opening + 1;
    else {
        dx = (n1->point->x) - (n2->point->x);
        dy = (n1->point->y) - (n2->point->y);
        dz = (n1->point->z) - (n2->point->z);
        d = sqrt( eps + dx*dx + dy*dy + dz*dz );
        
        if ( !(n1->is_aggregate || n2->is_aggregate) )
            opening = max_opening-1;
        else {
            if ( !(n1->is_aggregate) )
                opening = (n2->s)/d;
            else
                opening = (n1->s)/d; //( n1->s ) / ( d - (n1->s)/sqrt(3) );
        }
    }
    //Checks if the opening criteria is met to use cell-to-cell interactions
    //If not, recursively calculate acceleration between the child nodes 
    if ( opening < max_opening ) {
        double md3 = (n2->point->m)/(d*d*d);
        n1->point->ax -= dx*md3;
        n1->point->ay -= dy*md3;
        n1->point->az -= dz*md3;
    }
    else {
        if ( n1->is_aggregate == n2->is_aggregate ) {
            for (int i=0; i<8; i++) {
                for (int j= 0; j<8; j++) {
                    calculate_accelerations(n1->children[i], n2->children[j]);
                }
            }
        }
        else {
            if ( n1->is_aggregate ) {
                for (int i=0; i<8; i++)
                    calculate_accelerations(n1->children[i], n2);
            }
            if ( n2->is_aggregate ) {
                for (int i=0; i<8; i++)
                    calculate_accelerations(n1, n2->children[i]);
            }
        }
    }
}


void Octree::initialize_leapfrog(Node* n) {
    //Initializes the particle velocities for the leapfrog method
    if (n->is_aggregate) {
        Point* child;
        for (int i=0; i<8; i++) {
            child = n->children[i]->point;
            if (child != nullptr) {
                child->ax += n->point->ax;
                child->ay += n->point->ay;
                child->az += n->point->az;
                initialize_leapfrog(n->children[i]);
            }
        }
    }
    else if (n->point != nullptr) {
        n->point->vx -= (n->point->ax)*dt/2.0;
        n->point->vy -= (n->point->ay)*dt/2.0;
        n->point->vz -= (n->point->az)*dt/2.0;
    }
}

void Octree::update_positions(Node* n) {
    //Update the position of each leaf (i.e. particle) after recursively 
    //  accumulating the accelerations of all ancestor nodes
    if (n->is_aggregate) {
        Point* child;
        for (int i=0; i<8; i++) {
            child = n->children[i]->point;
            if (child != nullptr) {
                child->ax += n->point->ax;
                child->ay += n->point->ay;
                child->az += n->point->az;
                update_positions(n->children[i]);
            }
        }
    }
    else if (n->point != nullptr) {
        n->point->vx += (n->point->ax)*dt;
        n->point->vy += (n->point->ay)*dt;
        n->point->vz += (n->point->az)*dt;
        n->point->x += (n->point->vx)*dt;
        n->point->y += (n->point->vy)*dt;
        n->point->z += (n->point->vz)*dt;
    }
}
