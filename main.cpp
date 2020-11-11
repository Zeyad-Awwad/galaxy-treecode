
#include "Octree.hpp"

void print_positions( vector<Point> points, int N ) {
    for (int i=0; i<N; i++) {
        if (i > 0) cout << "\t";
        cout << points[i].x << "\t" << points[i].y;
    }
    cout << endl;
}


vector<Point> read_points(string filename)
{
    int i, lines;
    FILE * file;
    file = fopen(filename.c_str(), "r");
    fscanf(file, "Parameters:%d", &lines);
    vector<Point> points(lines);
    cout << "Reading " << lines << " lines" << endl;
    for (i=0; i<lines; i++)
    {
        fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &points[i].m, &points[i].x, &points[i].y, &points[i].z, &points[i].vx, &points[i].vy, &points[i].vz);
    }
    return points;
}

void write_points(string filename, vector<Point> points, double t)
{
    FILE * file;
    string key;
    file = fopen(filename.c_str(), "w");
    fprintf(file, "time=%lf", t);
    int i, N = (int) points.size();
    for (i = 0; i < N ; i++)
    {
        fprintf(file, "\n%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", points[i].x, points[i].y, points[i].z, points[i].vx, points[i].vy, points[i].vz, points[i].ax, points[i].ay, points[i].az);
    }
    fclose(file);
}

void accelerate_points(vector<Point> &points, bool initialized) {
    double d, dx, dy, dz;
    double eps = 0.000000;
    double dt = 0.0001;
    for (int i=0; i<points.size(); i++) {
        points[i].ax = 0;
        points[i].ay = 0;
        points[i].az = 0;
    }
    for (int i=0; i<points.size(); i++) {
        for (int j=i; j<points.size(); j++) {
            if (i == j) continue;
            dx = (points[i].x) - (points[j].x);
            dy = (points[i].y) - (points[j].y);
            dz = (points[i].z) - (points[j].z);
            d = sqrt( eps + dx*dx + dy*dy + dz*dz );
            double m1d3 = (points[i].m)/(d*d*d);
            double m2d3 = (points[j].m)/(d*d*d);
            points[i].ax -= dx*m2d3;
            points[i].ay -= dy*m2d3;
            points[i].az -= dz*m2d3;
            points[j].ax += dx*m1d3;
            points[j].ay += dy*m1d3;
            points[j].az += dz*m1d3;
        }
    }
    if (initialized) {
        for (int i=0; i<points.size(); i++) {
            points[i].vx += (points[i].ax)*dt;
            points[i].vy += (points[i].ay)*dt;
            points[i].vz += (points[i].az)*dt;
            points[i].x += (points[i].vx)*dt;
            points[i].y += (points[i].vy)*dt;
            points[i].z += (points[i].vz)*dt;
        }
    }
    else {
        for (int i=0; i<points.size(); i++) {
            points[i].vx -= (points[i].ax)*dt/2.0;
            points[i].vy -= (points[i].ay)*dt/2.0;
            points[i].vz -= (points[i].az)*dt/2.0;
        }
    }
}

int main(int argc, char **argv)
{
    int N = 4;
    double max_opening = 0.3;
    double eps = 0.025;
    double dt = 0.000001; 
    int interval = 20;
    Octree tree = Octree(-3, -3, -3, 4, 4, 4, max_opening, eps, dt);
    
    string directory = "/Users/Zeyad/Desktop/treecode/New/";
    string input_file = "data/galaxy.txt"; 
    vector<Point> points = read_points(directory + input_file );
    N = points.size();
    
    tree.add_point_vector( points );
    cout << "Calculating aggregates" << endl;
    tree.calculate_aggregates();
    Point* p = tree.root->point;
    cout << "Root has:\tm = " << p->m << "\tx = " << p->x << "\ty = " << p->y << "\tz = " << p->z << endl;
    
    tree.calculate_all_accelerations();
    tree.initialize_leapfrog(tree.root);
    tree.update_positions(tree.root);
    int counter = 0;
    for (int i=0; i<=1000; i++) {
        tree.clear(points);
        tree.add_point_vector(points);
        tree.calculate_aggregates();
        tree.calculate_all_accelerations();
        tree.update_positions(tree.root);
        if (i%interval == 0)
            cout << "Saving a snapshot at iteration" << i << endl;
            write_points(directory + "results/output_" + to_string(counter++) + ".txt", points, dt*i);
    }
    
    cout << "Finished running" << endl;
}
