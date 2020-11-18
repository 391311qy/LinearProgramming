/* The implementation of kamakar's algorithm
    eigen dependency:
    https://eigen.tuxfamily.org/dox/GettingStarted.html
    download and include the path to eigen folder
    g++ -I /path/to/eigen/ Kamakar.cpp -o km
    
    min cTx
    s.t. Ax <= b, x >= 0  
    usage: 
    solve(A, b, c)
*/
#include "LP.h"
#include "vector"
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
using namespace std;
using namespace Eigen;

/*
The implementation follows the original Karmarkar algorithm, the problem 
 solved here is slightly different as 
 min cTx; 
 s.t. Ax = 0
     eTx = 1
      x >= 0
*/
void Kamakar::solve_original(MatrixXd &AT, VectorXd &c) {
    // initialization AT: n * m, here is AT by default.
    int m = AT.cols();
    int n = AT.rows();

    // learning rate 
    double beta = 0.25;
    double r = 1/sqrt(m*(m-1));
    double gamma = 0.5;

    // initialize e, x0, I
    RowVectorXd eT = RowVectorXd::Ones(m); // col identity (m*1) 
    MatrixXd I = eT.asDiagonal(); // D(0) (m*m)
    VectorXd x = pow(m,-1)*eT.transpose(); // initialize x0 (m*1)
    VectorXd x0 = pow(n,-1)*eT.transpose(); // initialize x0 (m*1)

    // iteration params
    int i = 0; 
    int max_iter = 10000;

    // iterative find optimum
    VectorXd lastx;
    while (c.transpose()*x >= 0 && i < max_iter) {
        cout<<x.transpose()<<endl;
        MatrixXd D = x.asDiagonal(); // D(i) (m*m)
        MatrixXd B(n + 1, m);B << AT*D, eT;// BT: (n+1*m)
        VectorXd cp = - (I - B.transpose()*((B*B.transpose()).inverse())*B) * D * c;
        VectorXd z = x0 + beta*r*cp/cp.norm();
        lastx = x;
        x = D*z/(eT*D*z); // update x via projective transformation
        ++i; // increment i
    }
    cout<<"minimized result is "<<lastx.transpose()*c<<endl;
    cout<<"  "<<endl;
}


/* 
revised version of karmarkar's algorithm 
maximize cTx
s.t. Ax + Dv*v_h = b
    v >= 0 
where v = Dv*v_h >= 0
    Dv = diag(v) at kth iteration
    v is a slack variable.
*/
void Kamakar::solve_default(MatrixXd &AT, VectorXd &b, VectorXd &c, VectorXd &x0) {
    int m = AT.cols();
    int n = AT.rows();
    
    RowVectorXd eT = RowVectorXd::Ones(m); // col identity (m*1) 
    MatrixXd I = eT.asDiagonal(); // D(0) (m*m)
    VectorXd x = x0; // initialize x0 (m*1) require to be interior.

    // iteration params
    int i = 0; 
    int max_iter = 100;
    double gamma = 0.1;
    
    while ( i < max_iter) {
        //need to add stopping condition
        VectorXd v = b - AT*x; // slack v
        MatrixXd Dv = v.asDiagonal();
        VectorXd hx = (AT.transpose()*(Dv*Dv).inverse()*AT).inverse()*c;
        VectorXd hv = -AT*hx;
        double alpha = INT_MAX;
        for (int i = 0 ; i <hv.size(); i++) {
            if (hv[i] >= 0) {cout<<"unbounded"<<endl ;return;}
            alpha = min (-v[i]/hv[i], alpha);
        }
        x += gamma * alpha * hx;
        ++i;
    }
    cout<< "x is "<<x.transpose()<<endl;
    cout<< "Maximized result is " << c.transpose()*x<<endl;

}


