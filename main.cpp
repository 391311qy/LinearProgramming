#include "LP.h"

// To run: 
// eigen_path = C:/C++libs/eigen-3.3.8/eigen-3.3.8/ for example
// g++ -I C:/C++libs/eigen-3.3.8/eigen-3.3.8/ Simplex.cpp Kamakar.cpp main.cpp -o solver

int main() {
    vector<vector<double>> A_s = {
        {3.0, 2.0, 1.0},
        {2.0, 5.0, 3.0},
    };
    vector<double> b_s = {10.0, 15.0};
    vector<double> c_s = {2.0, 3.0, 4.0};

    MatrixXd A(2,3);
    A << 3.0, 2.0, 1.0, 
        2.0, 5.0, 3.0;
    VectorXd b(2);
    b << 10, 15;
    VectorXd c(3);
    c << -3.0, -2.0, -4.0; 

    Simplex sln_spx;
    Kamakar sln;
    sln.solve(A, b, c);
    sln_spx.simplexTableaux(A_s, b_s, c_s);
    sln_spx.solve();
    return 0;
}