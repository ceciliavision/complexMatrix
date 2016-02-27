#include <iostream>
#include <fstream>
#include <complex>
using namespace std;

#include "matrixComplex.hpp"

main(){
  srand48(time(NULL));
  // constructor
  matrix<fcomplex> A(4,4, "A");
  matrix<fcomplex> B(4,4, "B");

  // randomize entries
  A.randomize();
  B.randomize();

  // overloaded operator* for matrix-matrix multiplication
  matrix<fcomplex> C = A * B;

  // overloaded operator| for left matrix division
  matrix<fcomplex> D = A | C; // i.e. D = A\C in MATLAB
  matrix<fcomplex> E = B | D;

  // complicated expression with simple outcome
  matrix<fcomplex> F = (A+B) | ( (A+B) | ( (B+A)*A + (A+B)*B ) );

  // test the overloaded two-argument operator()
  F(4,1) = fcomplex (2.,0.);

  // stream out to standard out
  cout << A << endl;
  cout << B << endl;
  cout << C << endl;
  cout << D << endl;
  cout << E << endl;
  cout << F << endl;

}
