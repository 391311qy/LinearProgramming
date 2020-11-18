# LinearProgramming

Algorithms on linear programming

[1]:[Simplex tabuleau form ](https://link.springer.com/chapter/10.1007/978-0-8176-4844-2_2)<br />
[2]:[Karmarkar's interior point method](https://link.springer.com/article/10.1007/BF01587095)<br />

[eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html) dependency is used for matrix manipulation.<br />
  download the package and extract to eigen path, include the path when compiled with following cmd:<br />
```
g++ -I ...(YourPath)/eigen-3.3.8/eigen-3.3.8/ Simplex.cpp Kamakar.cpp main.cpp -o solver
```
