/*  This is an implementation of solvers
    for linear programming methods
    
    maximize cTx;
    subjected to ATx <= b
                x >= 0

    A: p-by-n matrix
    x: n vector
    b: p vector
*/
#include <iostream>
#include <list>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

class Simplex {
    /* simplex 
        iteratively finding the pivot 
        useful resources:
        https://jeremykun.com/2014/12/01/linear-programming-and-the-simplex-algorithm/*/
    public:
        /*The simplex tabuleaux
          a = [A(M,N) I(M,M) b
               C(1,N) 0(1,M) 0]*/
        void simplexTableaux(vector<vector<double>> &A, vector<double> &b, vector<double> &c);

        /* pivot function 
           used to scale all elements except a specific loc (p,q)*/
        void pivot(int &p, int &q);

        /* solve simplex */
        void solve();

    private:
        vector<vector<double>> a;
        int M, N;
};

class Kamakar {
    /*Kamakar's internal point methods*/
    public:
        void solve(MatrixXf &A, VectorXf &b, VectorXf &c, VectorXf &x0);
    private:
        int k;
        float gamma;
};