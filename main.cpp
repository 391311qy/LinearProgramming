#include "LP.h"
//////////////////////////////////////////////////////////////////////
////                     Genearal Linear Programing               ////
////                     max/min cTx                              ////
////                     s.t. ATx >= b                            ////
////                          eTx >= 1                            ////
//////////////////////////////////////////////////////////////////////

/* How to run: 
 This package uses Eigen. (see also README)
 eigen_path for example = C:/C++libs/eigen-3.3.8/eigen-3.3.8/ 
 g++ -I C:/C++libs/eigen-3.3.8/eigen-3.3.8/ Simplex.cpp Kamakar.cpp main.cpp -o solver
 */

int main() {
    // solving maximization problem using Simplex
    vector<vector<double>> A_s = {
        {3.0, 2.0, 1.0},
        {2.0, 5.0, 3.0},
    };
    vector<double> b_s = {10.0, 15.0};
    vector<double> c_s = {2.0, 3.0, 4.0};
    Simplex sln_spx;
    sln_spx.simplexTableaux(A_s, b_s, c_s);
    sln_spx.solve();

    // solving minimization problem using Simplex
    MatrixXd A(2,2);
    A << 3.0, 6.0, 
        3.0, 1.0;
    VectorXd b(2);
    b << 24, 9;
    VectorXd c(2);
    c << 2.0, 3.0; 
    Kamakar sln;
    sln.solve(A, b, c);
    
    return 0;
}