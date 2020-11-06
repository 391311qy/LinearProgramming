/* The implementation of simplex algorithm
    max cTx
    s.t. Ax <= b, x >= 0  
    usage: 
    simplexTableaux(A, b, c)
    solve()
*/
#include "LP.h"
#include "vector"
#include <iostream>
using namespace std;

void Simplex::simplexTableaux(vector<vector<double>> &A, vector<double> &b, vector<double> &c) {
    M = b.size();
    N = c.size();
    a = vector<vector<double>>( M + 1 , vector<double> (N+M+1+1, 0)); 
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            a[i][j] = A[i][j];
        }
    }
    for (int j = N; j < N + M; j++) { a[j-N][j] = 1.0; } // allocate identity
    for (int k = 0; k < N;     k++) { a[M][k] = c[k];} // allocate c
    for (int j = 0; j < M;     j++) { a[j][M+N] = b[j];} // allocate b
    // print input
    cout<<"input tab"<<endl;
    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < M+N+1; j++) {
            cout<<a[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"______"<<endl;
}

void Simplex::pivot(int &p, int &q) {
    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= M + N; j++) {
            if (i != p && j != q) {
                a[i][j] -= a[p][j] * a[i][q] / a[p][q];
            }
        }
    }
    // zero out col q
    for (int i = 0; i <= M; i++) {if (i != p) a[i][q] = 0.0;}
    // scale row p
    for (int j = 0; j <= M + N; j++) {
        if (j != q) { a[p][j] /= a[p][q];}}
    a[p][q] = 1.0;
    
}

void Simplex::solve() {
    int cc = 0;
    while (true) {
        cc++;
        if (cc == 1000000) {break;}
        int p, q;
        //termination upon positive objective function
        for (q = 0; q < M + N; q++) { if (a[M][q] > 0) { break;}}
        if (q > M + N) {break;}
        for (p = 0; p < M; p++) { if (a[p][q] > 0) break;}
        // find p based on min ratio rule
        for (int i = p+1; i < M; i++) { 
            if (a[i][q] > 0) {
                if (a[i][M+N] / a[i][q] < a[p][M+N] / a[p][q]) {
                    p = i;
                }
            }
        }
        pivot(p, q);
    }
    // print result
    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < M+N+1; j++) {
            cout<<a[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"maximized value is "<<-a[M][M+N]<<endl;
}

int main() {
    /* Example:
        5x + 15y <= 480
        4x + 4y <= 160
        35x + 20 y <= 1190
        maximize 13x + 23y 
    */
    Simplex Sln;
    vector<vector<double>> A = {
        {3.0, 2.0, 1.0},
        {2.0, 5.0, 3.0},
    };
    vector<double> b = {10.0, 15.0};
    vector<double> c = {2.0, 3.0, 4.0};
    Sln.simplexTableaux(A, b, c);
    Sln.solve();
    return 0;
}
