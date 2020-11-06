/* The implementation of kamakar's algorithm
    eigen dependency:
    https://eigen.tuxfamily.org/dox/GettingStarted.html
    download and include the path to eigen folder
    g++ -I /path/to/eigen/ Kamakar.cpp -o km
    
    max cTx
    s.t. Ax <= b, x >= 0  
    usage: 
    simplexTableaux(A, b, c)
    solve()
*/
#include "LP.h"
#include "vector"
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

void Kamakar::solve(MatrixXf &A, VectorXf &b, VectorXf &c, VectorXf &x0) {
    k = 0;
    gamma = 0.98;
    VectorXf x = x0;
    
    while (true) {
        //need to add stopping condition
        if (k > 10) {cout<<x<<endl; return;}
        VectorXf v = b - A*x;
        MatrixXf Dv2 = v.cwiseProduct(v).asDiagonal();
        VectorXf hx = (A.transpose()*Dv2.inverse()*A).inverse()*c;
        VectorXf hv = -A*hx;
        if (hv.minCoeff() >= 0) {cout<<"unbounded"<<endl; return;}
        float alpha = gamma * (-v.cwiseProduct(hv.cwiseInverse())).minCoeff();
        x += alpha * hx;
        ++k;
    }

}

int main() {
    MatrixXf A(2,3);
    A << 3.0, 2.0, 1.0,
        2.0, 5.0, 3.0;
    VectorXf b(2);
    b << 10.0, 
        15.0;
    VectorXf c(3);
    c << 2.0, 
        3.0, 
        4.0;
    VectorXf x0(3);
    x0 << 0.0, 
        0.0, 
        0.0;
    Kamakar sln;
    sln.solve(A, b, c, x0);
    return 0;
}


