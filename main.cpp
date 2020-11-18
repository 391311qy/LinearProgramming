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
    // vector<vector<double>> A_s = {
    //     {3.0, 2.0, 1.0},
    //     {2.0, 5.0, 3.0},
    // };
    // vector<double> b_s = {10.0, 15.0};
    // vector<double> c_s = {2.0, 3.0, 4.0};
    // Simplex sln_spx;
    // sln_spx.simplexTableaux(A_s, b_s, c_s);
    // sln_spx.solve();

    // solving minimization problem using Karmarkar's algo
    // notice that A must be full rank (i.e. m > n)
    MatrixXd A(11,2);
    A << 0.0, 1, 0.2, 1, 0.4, 1, 0.6, 1, 0.8, 1, 
        1.0, 1, 1.2, 1, 1.4, 1, 1.6, 1, 1.8, 1, 2.0, 1;
    VectorXd b(11);
    b << 1.0, 1.01, 1.04, 1.09, 1.16, 1.25, 1.36, 1.49, 1.64, 1.81, 2.0; 
    VectorXd c(2);
    c<< 1.0 , 1.0; 
    VectorXd x0(2);
    x0<< 0.0, 0.0;
    Kamakar sln;
    sln.solve_default(A, b, c, x0);
    
    return 0;
}