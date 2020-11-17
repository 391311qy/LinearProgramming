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
using namespace std;
using namespace Eigen;

void Kamakar::solve(MatrixXd &AT, VectorXd &b, VectorXd &c) {
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
    while (c.transpose()*x <= 0 && i < max_iter) {
        MatrixXd D = x.asDiagonal(); // D(i) (m*m)
        MatrixXd B(n + 1, m);B << AT*D, eT;// BT: (n+1*m)
        VectorXd cp = - (I - B.transpose()*((B*B.transpose()).inverse())*B) * D * c;
        VectorXd z = x0 + beta*r*cp/cp.norm();
        x = D*z/(eT*D*z); // update x via projective transformation
        ++i; // increment i
    }
    cout<<"minimized result is "<<x.transpose()*c<<endl;
}


void Kamakar::solve_maximize(MatrixXd &AT, VectorXd &b, VectorXd &c) {
    int m = AT.cols();
    int n = AT.rows();
    
    RowVectorXd eT = RowVectorXd::Ones(m); // col identity (m*1) 
    MatrixXd I = eT.asDiagonal(); // D(0) (m*m)
    VectorXd x = pow(m, -1) *eT.transpose(); // initialize x0 (m*1)

    // iteration params
    int i = 0; 
    int max_iter = 100;
    double gamma;
    cout<< "what is learning rate ?" <<endl;
    cin >> gamma;
    
    while ( i < max_iter) {
        //need to add stopping condition
        VectorXd v = b - AT*x;
        MatrixXd Dv2 = v.cwiseProduct(v).asDiagonal();
        VectorXd hx = (AT.transpose()*Dv2.inverse()*AT).inverse()*c;
        VectorXd hv = -AT*hx;
        double alpha = 0;
        for (int i = 0 ; i <hv.size(); i++) {
            if (hv[i] >= 0) {cout<<"unbounded"<<endl ;return;}
            alpha = min (-v[i]/hv[i], alpha);
        }
        x += gamma * alpha * hx;
        cout<<x.transpose()<<endl;
        ++i;
    }
    cout << c.transpose()*x<<endl;

}


