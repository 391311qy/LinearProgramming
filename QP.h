/*  This is an implementation of solvers
    for quadric programming methods
    
    minimize 1/2 xTQx + cTx
    s.t.    Ax <= b

    A: m-by-n matrix
    Q: n-by-n matrix
    c: n vector
    x: n vector
    b: m vector
*/
#include <iostream>
#include <list>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

class ConjugateGradient {
    public:
    private:
};