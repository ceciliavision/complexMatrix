#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

extern "C"
{
  extern void dgesv_(int*,int*, double*,int*,int*,
		        double*,int*,int*);
  extern void sgesv_(int*, int*, float*, int*, int*,
		     float*, int*, int*);
  extern void dgemm_(char*,char*,int*,int*,int*,
		     double*, const double*,int*, const double*,int*,
		     double*,double*,int*);
  extern void sgemm_(char*, char*, int*, int*, int*,
		    float*, const float*, int*, const float*, int*,
		    float*, float*, int*);
};

template <class T>
class matrix{
private:
  int _Nrow;
  int _Ncolumn; 
  string _name; // name of the matrix for later reference
  T *_data;

public:
  // default constructor
  matrix(){
    cout << "In: default vector constructor" 
	 << endl;
    _Nrow = 0;
    _Ncolumn = 0;
    _name = "";
    _data = NULL;
  };

  //constructor with input size
  matrix(int row, int col, string str){
    _Nrow = row;
    _Ncolumn = col;
    _name = str;
    _data = new T[row*col];
    for(int n=0;n<row*col;++n)
      _data[n] = 0;
  };

  //copy constructor: 
  //Initializes a matrix using the dimensions and data from a second matrix
  matrix(const matrix &a){    
    _Nrow = a.row();
    _Ncolumn = a.col();

    _data = new T[_Nrow*_Ncolumn];
    for(int n=0;n<_Nrow*_Ncolumn;++n)
      _data[n] = a._data[n];
  };

  // destructor
  ~matrix(){
    cout << "In: ~matrix" << endl;
    if(_data != NULL)
      delete [] _data;
    
  };

  // randomize
  void randomize(){
    _name = "rand(" + this->name() + ")";
    for(int n=0;n<_Nrow*_Ncolumn;++n){
      _data[n] = drand48();
    }
  };

  // operator()
  // return a reference to an entry in the matrix at a requested row and column
  T &operator() (int row, int col) const{
    if(_data==NULL &&
       row>_Nrow || 
       row<=0 && col>_Ncolumn || col<=0){
      
      cout << "matrix::get warning matrix NULL" << endl;
      exit(-1);
    }

    return _data[_Nrow*(col-1)+(row-1)];
  };

  // operator +
  friend matrix operator+(const matrix &a, const matrix &b){
    matrix aPb(a._Nrow, a._Ncolumn, "("+a._name+"+"+b._name+")");
    int r = a.row();
    int c = a.col();
    for(int n=1;n<=r;++n){
      for(int m=1;m<=c;++m){
	aPb(n,m) = a(n,m)+b(n,m);
      }
    }
    return aPb;
  };

  int row() const{
    return _Nrow;
  };

  int col() const{ 
    return _Ncolumn;
  };

  string name() const{
    return _name;
  };

  template <class U>
  friend matrix<U> operator * (const matrix<U> &a, const matrix<U> &b);

  template <class U>
  friend matrix<U> operator | (const matrix<U> &a, const matrix<U> &b);
};

template<class T>
matrix<T> operator * (const matrix<T> &a, const matrix<T> &b);

template<class T>
matrix<T> operator | (const matrix<T> &a, const matrix<T> &b);


// operator *
template<>
matrix<float> operator * (const matrix<float> &a, const matrix<float> &b){
  matrix<float> aMb(a._Nrow, b._Ncolumn, "("+a._name+"*"+b._name+")");
  matrix<float> A(a);
  matrix<float> B(b);
  int M = a._Nrow;
  int N = b._Ncolumn;
  int K = a._Ncolumn;
  int LDA = M;
  int LDB = b._Nrow;
  int LDC = aMb._Nrow;
  float alpha = 1.0;
  float beta = 0.0;
  char transA = 'N';
  char transB = 'N';
    sgemm_(&transA, &transB,
	   &M, &N, &K,
	   &alpha, A._data, &LDA, 
	   B._data, &LDB, &beta,
	   aMb._data,
	   &LDC);
    return aMb;
};

template<>
matrix<double> operator * (const matrix<double> &a, const matrix<double> &b){
  matrix<double> aMb(a._Nrow, b._Ncolumn, "("+a._name+"*"+b._name+")");
  matrix<double> A(a);
  matrix<double> B(b);
  int M = a._Nrow;
  int N = b._Ncolumn;
  int K = a._Ncolumn;
  int LDA = M;
  int LDB = b._Nrow;
  int LDC = aMb._Nrow;
  double alpha = 1.0;
  double beta = 0.0;
  char transA = 'N';
  char transB = 'N';
    dgemm_(&transA, &transB,
	   &M, &N, &K,
	   &alpha, A._data, &LDA, 
	   B._data, &LDB, &beta,
	   aMb._data,
	   &LDC);
    return aMb;
};

// operator |
template<>
matrix<float> operator | (const matrix<float> &a, const matrix<float> &b){
  int N, NRHS,LDA,LDB,INFO;
  int *IPIV;
  matrix<float> A(a);
  matrix<float> B(b);
  B._name = "("+a.name()+"|"+b.name()+")";
  N = a.row();
  NRHS = b.col();
  LDA = N;
  LDB = NRHS;
  IPIV = (int*)calloc(a.row(),sizeof(a.row()));
  sgesv_(&N, &NRHS, 
	 A._data, &LDA, IPIV, 
	 B._data, &LDB, &INFO);
  free(IPIV);
  return B;
};

template<>
matrix<double> operator | (const matrix<double> &a, const matrix<double> &b){
  int N, NRHS,LDA,LDB,INFO;
  int *IPIV;
  matrix<double> A(a);
  matrix<double> B(b);
  B._name = "("+a.name()+"|"+b.name()+")";
  N = a.row();
  NRHS = b.col();
  LDA = N;
  LDB = NRHS;
  IPIV = (int*)calloc(a.row(),sizeof(a.row()));
  dgesv_(&N, &NRHS, 
	 A._data, &LDA, IPIV, 
	 B._data, &LDB, &INFO);
  free(IPIV);
  return B;
};
 
// outside of the class
template <class T>
ostream & operator<<(ostream &os, const matrix<T> &a){
  os << a.name() << endl;
  os << a.row() << " ";
  os << a.col() << endl;
  for(int n=1;n<=a.row();++n){
    for(int m=1;m<=a.col();++m){
      os << a(n,m) << " ";
    }
    os << endl;
  }
  return os;
};


