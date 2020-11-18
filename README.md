# LinearProgramming

Algorithms on linear programming
[Simplex tabuleau form ](https://link.springer.com/chapter/10.1007/978-0-8176-4844-2_2);
[Karmarkar's interior point method](https://link.springer.com/article/10.1007/BF01587095)

eigen dependency:
    https://eigen.tuxfamily.org/dox/GettingStarted.html
download to eigen path and include when compiled with following cmd:

g++ -I ...(YourPath)/eigen-3.3.8/eigen-3.3.8/ Simplex.cpp Kamakar.cpp main.cpp -o solver
