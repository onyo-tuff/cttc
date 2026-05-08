#include "Mesh.h"


////// Public functions
int Mesh::node_number(int i, int j) {
    int n;
    n = i*M + j;
    return n;
}

std::vector<int> Mesh::ij_given_n(int n) {
    std::vector<int> indx = {0, 0};
    indx[0] = (int)(n/M);
    indx[1] = n - indx[0]*M;
    return indx;
}

// Read inputs
void Mesh::getInputs() {
    while (std::getline(input, strInput))
    {
        if (strInput == "Domain Width 'L' (m):") {
            std::getline(input, strInput);
            L = std::stod(strInput);
        }
        if (strInput == "Domain Height 'H' (m):") {
            std::getline(input, strInput);
            H = std::stod(strInput);
        }
        if (strInput == "Minimum horizontal element size 'hx_min' (m):") {
            std::getline(input, strInput);
            hx_min = std::stod(strInput);
        }
        if (strInput == "Maximum horizontal element size for x1 'hx_max1' (m):") {
            std::getline(input, strInput);
            hx_max1 = std::stod(strInput);
        }
        if (strInput == "Maximum horizontal element size for x2 'hx_max2' (m):") {
            std::getline(input, strInput);
            hx_max2 = std::stod(strInput);
        }
        if (strInput == "Maximum horizontal element size for x3 'hx_max3' (m):") {
            std::getline(input, strInput);
            hx_max3 = std::stod(strInput);
        }
        if (strInput == "Minimum vertical element size 'hy_min' (m):") {
            std::getline(input, strInput);
            hy_min = std::stod(strInput);
        }
        if (strInput == "Maximum vertical element size for y1 'hy_max1' (m):") {
            std::getline(input, strInput);
            hy_max1 = std::stod(strInput);
        }
        if (strInput == "Maximum vertical element size for y2 'hy_max2' (m):") {
            std::getline(input, strInput);
            hy_max2 = std::stod(strInput);
        }
        if (strInput == "Maximum vertical element size for y3 'hy_max3' (m):") {
            std::getline(input, strInput);
            hy_max3 = std::stod(strInput);
        }
        if (strInput == "Transition smoothness parameter along x1 'alpha_x1':") {
            std::getline(input, strInput);
            alpha_x1 = std::stod(strInput);
        }
        if (strInput == "Transition smoothness parameter along x2 'alpha_x2':") {
            std::getline(input, strInput);
            alpha_x2 = std::stod(strInput);
        }
        if (strInput == "Transition smoothness parameter along x3 'alpha_x3':") {
            std::getline(input, strInput);
            alpha_x3 = std::stod(strInput);
        }
        if (strInput == "Transition smoothness parameter along y1 'alpha_y1':") {
            std::getline(input, strInput);
            alpha_y1 = std::stod(strInput);
        }
        if (strInput == "Transition smoothness parameter along y2 'alpha_y2':") {
            std::getline(input, strInput);
            alpha_y2 = std::stod(strInput);
        }
        if (strInput == "Transition smoothness parameter along y3 'alpha_y3':") {
            std::getline(input, strInput);
            alpha_y3 = std::stod(strInput);
        }
        if (strInput == "x_r (m):") {
            std::getline(input, strInput);
            x_r = std::stod(strInput);
        }
        if (strInput == "y_r (m):") {
            std::getline(input, strInput);
            y_r = std::stod(strInput);
        }
        if (strInput == "Width of the obstacle 'L_r' (m):") {
            std::getline(input, strInput);
            L_r = std::stod(strInput);
        }
        if (strInput == "Height of the obstacle 'H_r' (m):") {
            std::getline(input, strInput);
            H_r = std::stod(strInput);
        }
    }
}

// Discretize functions
std::vector<double> Mesh::Discretize(double h_min, double h_max, double alpha, double x0, double xf) {
    std::vector<double> x = { 0 };
    double h, k, n=0, L=std::abs(xf-x0);
    if (alpha == 0) {
        while (x.back() < L) {
            n++;
            x.push_back(h_min*n);
        }
        k = L/x.back();
        for (size_t i = 0; i < x.size(); i++) { x[i] =k*x[i]; }
    }
    else {
        while (x.back() < L/2) {
            h = h_min + (h_max - h_min)*tanh(alpha*x.back());
            x.push_back(x.back() + h);
        }
        k = 0.5*L/x.back();
        for (size_t i = 0; i < x.size(); i++) { x[i] = k*x[i]; }
        std::vector<double> aux = x;
        aux.pop_back();
        for (size_t i = 0; i < aux.size(); i++) { aux[i] = L - x[i]; }
        std::reverse(aux.begin(), aux.end());
        x.insert(x.end(), aux.begin(), aux.end());
        for (size_t i = 0; i < x.size(); i++) { x[i] = x0 + (xf-x0)/L*x[i]; }
    }
    return x;
}

std::vector<double> Mesh::Discretize1(double h_min, double h_max, double alpha, double x0, double xf) {
    std::vector<double> x = { 0 };
    double h, L=std::abs(xf-x0), k;
    if (alpha == 0) {
        x = Discretize(h_min,0,0,x0,xf);
    }
    else {
        while (x.back() < L) {
            h = h_min + (h_max - h_min)*tanh(alpha*x.back());
            x.push_back(x.back() + h);
        }
        k = L/x.back();
        for (size_t i=0; i < x.size(); i++) { 
            x[i] = x0 + (xf-x0)/L*k*x[i]; 
        }
    }
    if (x0 > xf) { std::reverse(x.begin(), x.end()); }
    return x;
}

// Default constructor
Mesh::Mesh() {
    getInputs();
    Build();
    buildArrays();
    buildVxMesh();
    buildVyMesh();
    buildIdentifiers(false);
    N_vol = N*M;
    N_obstacle = 0;
}

void Mesh::Build() {
    x_coords = Discretize(hx_min, hx_max1, alpha_x1, 0, L);
    y_coords = Discretize(hy_min, hy_max1, alpha_y1, 0, H);
    N = static_cast<int>(x_coords.size()) - 1;
    M = static_cast<int>(y_coords.size()) - 1;
}

// Internal rectangle constructor
Mesh::Mesh(bool rectangle) {
    getInputs();
    Build1();
    buildArrays();
    buildVxMesh();
    buildVyMesh();
    buildIdentifiers(true);
    N_vol = N*M;
}

void Mesh::Build1() {
    std::vector<double> x1, x2, x3, y1, y2, y3;
    x1 = Discretize1(hx_min, hx_max1, alpha_x1, x_r, 0);
    x2 = Discretize(hx_min, hx_max2, alpha_x2, x_r, x_r+L_r);
    x3 = Discretize1(hx_min, hx_max3, alpha_x3, x_r+L_r, L);
    y1 = Discretize1(hy_min, hy_max1, alpha_y1, y_r, 0);
    y2 = Discretize(hy_min, hy_max2, alpha_y2, y_r, y_r+H_r);
    y3 = Discretize1(hy_min, hy_max3, alpha_y3, y_r+H_r, H);
    N1 = static_cast<int>(x1.size()) - 1;
    N2 = static_cast<int>(x2.size()) - 1;
    N3 = static_cast<int>(x3.size()) - 1;
    M1 = static_cast<int>(y1.size()) - 1;
    M2 = static_cast<int>(y2.size()) - 1;
    M3 = static_cast<int>(y3.size()) - 1;
    x1.pop_back();
    x2.pop_back();
    x_coords.insert(x_coords.end(), x1.begin(), x1.end());
    x_coords.insert(x_coords.end(), x2.begin(), x2.end());
    x_coords.insert(x_coords.end(), x3.begin(), x3.end());
    y1.pop_back();
    y2.pop_back();
    y_coords.insert(y_coords.end(), y1.begin(), y1.end());
    y_coords.insert(y_coords.end(), y2.begin(), y2.end());
    y_coords.insert(y_coords.end(), y3.begin(), y3.end());
    N = static_cast<int>(x_coords.size()) - 1;
    M = static_cast<int>(y_coords.size()) - 1;
}

// For both cases (after the initial build determining the x and y discretization)
void Mesh::buildArrays() {
    std::vector<std::vector<double>> x(5, std::vector<double>(2,0));
    std::vector<int> node(2,0);
    std::vector<double> VnS(5,0);
    n_nod = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            n_nod++;
            x[0] = Center(i,j); // center
            x[1] = West(i,j); // west
            x[2] = East(i,j); // east
            x[3] = South(i,j); // south
            x[4] = North(i,j); // north
            Coords.push_back(x);
            node[0] = i;
            node[1] = j;
            Connect.push_back(node);
            VnS[1] = std::abs(y_coords[j+1] - y_coords[j]); // S_w
            VnS[2] = VnS[1]; // S_e
            VnS[3] = std::abs(x_coords[i+1] - x_coords[i]); // S_s
            VnS[4] = VnS[3]; // S_n
            VnS[0] = VnS[1]*VnS[3]; // Vol
            Vols_n_Surfaces.push_back(VnS);
        }
    }
}

std::vector<double> Mesh::Center(int i, int j) {
    double x_c, y_c;
    x_c = 0.5*x_coords[i+1] + 0.5*x_coords[i];
    y_c = 0.5*y_coords[j+1] + 0.5*y_coords[j];
    std::vector<double> center = { x_c, y_c };
    return center;
}

std::vector<double> Mesh::West(int i, int j) {
    double x_w, y_w;
    x_w = x_coords[i];
    y_w = 0.5*y_coords[j+1] + 0.5*y_coords[j];
    std::vector<double> west = { x_w, y_w };
    return west;
}

std::vector<double> Mesh::East(int i, int j) {
    double x_e, y_e;
    x_e = x_coords[i+1];
    y_e = 0.5*y_coords[j+1] + 0.5*y_coords[j];
    std::vector<double> east = { x_e, y_e };
    return east;
}

std::vector<double> Mesh::South(int i, int j) {
    double x_s, y_s;
    x_s = 0.5*x_coords[i+1] + 0.5*x_coords[i];
    y_s = y_coords[j];
    std::vector<double> south = { x_s, y_s };
    return south;
}

std::vector<double> Mesh::North(int i, int j) {
    double x_n, y_n;
    x_n = 0.5*x_coords[i+1] + 0.5*x_coords[i];
    y_n = y_coords[j+1];
    std::vector<double> north = { x_n, y_n };
    return north;
}

bool Mesh::inObstacle(int i, int j) {
    std::vector<double> center = Center(i,j);
    bool inside = false;
    double x_r1 = x_r + L_r;
    double y_r1 = y_r + H_r;
    if ((center[0] > x_r) && (center[0] < x_r1)) {
        if ((center[1] > y_r) && (center[1] < y_r1)) {
            inside = true;
        }
    }
    return inside;
}

bool Mesh::inObstacleVx(int i, int j) {
    std::vector<double> x = Coords_Vx[i][j];
    bool inside = false;
    double x_r1 = x_r + L_r;
    double y_r1 = y_r + H_r;
    if ((x[0] > x_r) && (x[0] < x_r1)) {
        if ((x[1] > y_r) && (x[1] < y_r1)) {
            inside = true;
        }
    }
    return inside;
}

bool Mesh::inObstacleVy(int i, int j) {
    std::vector<double> x = Coords_Vy[i][j];
    bool inside = false;
    double x_r1 = x_r + L_r;
    double y_r1 = y_r + H_r;
    if ((x[0] > x_r) && (x[0] < x_r1)) {
        if ((x[1] > y_r) && (x[1] < y_r1)) {
            inside = true;
        }
    }
    return inside;
}

void Mesh::buildIdentifiers(bool yesObstacle) {
    std::vector<std::string> vec(1,"Interior");
    for (int n = 0; n < n_nod; n++) {
        Identif.push_back(vec);
    }      
    for (int i = 0; i <= N; i++) {
        vec.clear();
        for (int j = 0; j <= M; j++) {
            vec.push_back("Interior");
        }
        VxMesh_Ids.push_back(vec);
        VyMesh_Ids.push_back(vec);
    }            
    std::vector<int> elems;
    double x1, y1, x2, y2;
    std::string name;
    bool read = false;
    std::string incr;
    // Add boundary identifiers
    while (Boundary_Input >> strInput) {
        if (strInput == "----------------------------------") 
            read = true;
        if (strInput == "(" && read) {
            Boundary_Input >> strInput;
            x1 = std::stod(strInput);
            Boundary_Input >> strInput;
            y1 = std::stod(strInput);
            Boundary_Input >> strInput;
            x2 = std::stod(strInput);
            Boundary_Input >> strInput;
            y2 = std::stod(strInput);
            Boundary_Input >> incr;
            Boundary_Input >> strInput;
            Boundary_Input >> name;
            elems = findBoundaryElems(x1,y1,x2,y2,incr);
            buildVxVyIds(x1, y1, x2, y2, name);
            for (size_t i=0; i<elems.size(); i++) { 
                if (Identif[elems[i]][0] == "Interior")
                    Identif[elems[i]][0] = name;
                else
                    Identif[elems[i]].push_back(name);
            }
        }
    }
    // Add obstacle identifiers
    if (yesObstacle) {
    N_obstacle = 0;
    int n;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
           if (inObstacle(i,j)) {
		N_obstacle++;
                n = node_number(i,j);
                if (Identif[n][0] == "Interior")
                    Identif[n][0] = "Obstacle";
                else
                    Identif[n].push_back("Obstacle");
           }
        }
    }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            if (inObstacleVx(i,j)) {
                VxMesh_Ids[i][j] = "Obstacle";
            }
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            if (inObstacleVy(i,j)) {
                VyMesh_Ids[i][j] = "Obstacle";
            }
        }
    }
}

std::vector<int> Mesh::findBoundaryElems(double x1, double y1, double x2, double y2, std::string incr) {
    std::vector<int> i_bound, j_bound, elems;
    int k=0;
    bool stop = false;
    double dist1, dist2;
    // Case 1 -> vertical boundary
    if (x1 == x2) {
        dist1 = std::abs(x_coords[0]-x1);
        while (!stop) {
            k++;
            dist2 = std::abs(x_coords[k] - x1);
            if (k==N) { i_bound.push_back(k-1); stop = true; }
            if (dist2 > dist1) {
                if (incr == "+") { i_bound.push_back(k - 1); }
                if (incr == "-") { i_bound.push_back(k - 2); }
                stop = true;
            }
            dist1 = dist2;
        }
        stop = false;
        k = 0;
        dist1 = std::abs(y_coords[0]-y1);
        while (!stop) {
            k++; 
            dist2 = std::abs(y_coords[k] - y1);
            if (k==M) { j_bound.push_back(k-1); stop = true; }
            if (dist2 > dist1) { j_bound.push_back(k-1); stop = true; }
            dist1 = dist2;
        }
        stop = false;
        dist1 = std::abs(y_coords[k] - y2);
        while (!stop) {
            k++;
            j_bound.push_back(k-1);
            dist2 = std::abs(y_coords[k] - y2);
            if (k==M) { stop = true; }
            if (dist2 > dist1) { j_bound.pop_back(); stop = true; }
            dist1 = dist2;
        }   
    }
    // Case 2 -> horizontal boundary
    if (y1 == y2) {
        stop = false;
        k = 0;
        dist1 = std::abs(y_coords[0]-y1);
        while (!stop) {
            k++;
            dist2 = std::abs(y_coords[k] - y1);
            if (k==M) { j_bound.push_back(k-1); stop=true; }
            if (dist2 > dist1) {
                if (incr == "+") { j_bound.push_back(k-1); }
                if (incr == "-") { j_bound.push_back(k-2); }
                stop = true;
            }
            dist1 = dist2;
        }
        stop = false;
        k = 0;
        dist1 = std::abs(x_coords[0]-x1);
        while (!stop) {
            k++; 
            dist2 = std::abs(x_coords[k] - x1);
            if (k==N) { i_bound.push_back(k-1); stop = true; }
            if (dist2 > dist1) { i_bound.push_back(k-1); stop = true; }
            dist1 = dist2;
        }
        stop = false;
        dist1 = std::abs(x_coords[k]-x2);
        while (!stop) {
            k++;
            i_bound.push_back(k-1);
            dist2 = std::abs(x_coords[k] - x2);
            if (k==N) { stop = true; }
            if (dist2 > dist1) { i_bound.pop_back(); stop = true; }
            dist1 = dist2;
        }   
    }
    // Find elements corresponding to the boundary delimited by the x1,y1,x2,y2 given
    for (size_t a = 0; a < i_bound.size(); a++) {
        for (size_t b = 0; b < j_bound.size(); b++) {
            elems.push_back( node_number( i_bound[a], j_bound[b]));
        }
    }
    return elems;
}

bool Mesh::hasIdentifier(int i, int j, std::string Id) {
    bool has = false;
    int n;
    n = node_number(i,j);
    for (size_t k = 0; k < Identif[n].size(); k++) {
        if (Identif[n][k] == Id) { has = true; }
    }
    return has;
}

bool Mesh::hasIdentifier(int n, std::string Id) {
    bool has = false;
    for (size_t k = 0; k < Identif[n].size(); k++) {
        if (Identif[n][k] == Id) { has = true; }
    }
    return has;
}

void Mesh::buildVxMesh() {
    std::vector<double> a(2,0);
    std::vector<std::vector<double>> vec;
    for (int j = 0; j < M; j++) {
            vec.push_back(a);
    }
    for (int i = 0; i<= N; i++) {
        Coords_Vx.push_back(vec);
        SnV_Vx.push_back(vec);
    }        
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j < M; j++) {
            Coords_Vx[i][j][0] = x_coords[i];
            Coords_Vx[i][j][1] = 0.5*(y_coords[j+1] + y_coords[j]);
            SnV_Vx[i][j][0] = std::abs(y_coords[j+1] - y_coords[j]);
            if (i==0) 
                SnV_Vx[i][j][1] = SnV_Vx[i][j][0]*getCenter(i,j)[0];
            else if (i==N)
                SnV_Vx[i][j][1] = SnV_Vx[i][j][0]*(std::abs(L - getCenter(i-1,j)[0]));
            else
                SnV_Vx[i][j][1] = SnV_Vx[i][j][0]*(std::abs(getCenter(i,j)[0] - getCenter(i-1,j)[0]));
        }
    }
}

void Mesh::buildVyMesh() {
    std::vector<double> a(2,0);
    std::vector<std::vector<double>> vec;
    for (int j = 0; j <= M; j++) {
            vec.push_back(a);
    }
    for (int i = 0; i< N; i++) {
        Coords_Vy.push_back(vec);
        SnV_Vy.push_back(vec);
    }  

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= M; j++) {
            Coords_Vy[i][j][0] = 0.5*(x_coords[i+1] + x_coords[i]);
            Coords_Vy[i][j][1] = y_coords[j];
            SnV_Vy[i][j][0] = std::abs(x_coords[i + 1] - x_coords[i]);
            if (j==0)
                SnV_Vy[i][j][1] = SnV_Vy[i][j][0]*getCenter(i,j)[1];
            else if (j==M)
                SnV_Vy[i][j][1] = SnV_Vy[i][j][0]*(std::abs(H-getCenter(i,j-1)[1]));
            else
                SnV_Vy[i][j][1] = SnV_Vy[i][j][0]*(std::abs(getCenter(i,j)[1]-getCenter(i,j-1)[1]));
        }
    }
}

void Mesh::buildVxVyIds(double x1, double y1, double x2, double y2, std::string name) {
    double dist1, dist2;
    bool stop = false;
    int i = 0, j = 0;
    if (x1 == x2) {
        dist1 = std::abs(x_coords[0]-x1);
        while (!stop) {
            i++;
            dist2 = std::abs(x_coords[i]-x1);
            if (i==N) { stop = true; }
            if (dist2>dist1) { i = i-1; stop = true;}
            dist1 = dist2;
        }
        stop = false;
        dist1 = std::abs(y_coords[0]-y1);
        while(!stop) {
            j++;
            dist2 = std::abs(y_coords[j]-y1);
            if (j==M) { VxMesh_Ids[i][j] = name; stop = true;}
            if (dist2>dist1) { VxMesh_Ids[i][j] = name; stop = true; }
            dist1 = dist2;
        }
        stop = false;
        dist1 = std::abs(y_coords[j]-y2);
        while (!stop) {
            j++;
            dist2 = std::abs(y_coords[j]-y2);
            if (j==M) { VxMesh_Ids[i][j] = name; stop = true;}
            if (dist2>dist1) { stop = true; }
            else { VxMesh_Ids[i][j] = name; }
            dist1 = dist2;
        }
    }
    i = 0; j = 0;
    if (y1 == y2) {
        stop = false;
        dist1 = std::abs(y_coords[0]-y1);
        while (!stop) {
            j++;
            dist2 = std::abs(y_coords[j]-y1);
            if (j==M) { stop = true; }
            if (dist2>dist1) { j = j-1; stop = true;}
            dist1 = dist2;
        }
        stop = false;
        dist1 = std::abs(x_coords[0]-x1);
        while(!stop) {
            i++;
            dist2 = std::abs(x_coords[i]-x1);
            if (i==N) { VyMesh_Ids[i][j] = name; stop = true;}
            if (dist2>dist1) { VyMesh_Ids[i][j] = name; stop = true; }
            dist1 = dist2;
        }
        stop = false;
        dist1 = std::abs(x_coords[i]-x2);
        while (!stop) {
            i++;
            dist2 = std::abs(x_coords[i]-x2);
            if (i==N) { VyMesh_Ids[i][j] = name; stop = true;}
            if (dist2>dist1) { stop = true; }
            else { VyMesh_Ids[i][j] = name; }
            dist1 = dist2;
        }
    }
}
