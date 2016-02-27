#include <iostream>
#include <fstream>
using namespace std;

#include "matrix.hpp"

main(){
  srand48(time(NULL));
  // constructor
  matrix<float> A(4,4, "A");
  matrix<float> B(4,4, "B");

  // randomize entries
  A.randomize();
  B.randomize();

  // overloaded operator* for matrix-matrix multiplication
  matrix<float> C = A * B;

  // overloaded operator| for left matrix division
  matrix<float> D = A | C; // i.e. D = A\C in MATLAB
  matrix<float> E = B | D;

  // complicated expression with simple outcome
  matrix<float> F = (A+B) | ( (A+B) | ( (B+A)*A + (A+B)*B ) );

  // test the overloaded two-argument operator()
  F(4,1) = 2.;

  // stream out to standard out
  cout << A << endl;
  cout << B << endl;
  cout << C << endl;
  cout << D << endl;
  cout << E << endl;
  cout << F << endl;

}
