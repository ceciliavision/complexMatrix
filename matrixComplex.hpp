#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <complex>

typedef complex<double> dcomplex;
typedef complex<float> fcomplex;

using namespace std;

extern "C"
{
  extern void cgesv_(int*,int*, fcomplex*,int*,int*,
		        fcomplex*,int*,int*);
  extern void zgesv_(int*, int*, dcomplex*, int*, int*,
		     dcomplex*, int*, int*);
  extern void cgemm_(char*,char*,int*,int*,int*,
		     fcomplex*, const fcomplex*,int*, const fcomplex*,int*,
		     fcomplex*,fcomplex*,int*);
  extern void zgemm_(char*, char*, int*, int*, int*,
		    dcomplex*, const dcomplex*, int*, const dcomplex*, int*,
		    dcomplex*, dcomplex*, int*);
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
      _data[n] = (0,0);
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
      _data[n] = T (drand48(),drand48());
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
matrix<fcomplex> operator * (const matrix<fcomplex> &a, const matrix<fcomplex> &b){
  matrix<fcomplex> aMb(a._Nrow, b._Ncolumn, "("+a._name+"*"+b._name+")");
  matrix<fcomplex> A(a);
  matrix<fcomplex> B(b);
  int M = a._Nrow;
  int N = b._Ncolumn;
  int K = a._Ncolumn;
  int LDA = M;
  int LDB = b._Nrow;
  int LDC = aMb._Nrow;
  fcomplex alpha = fcomplex (1,0);
  fcomplex beta = fcomplex (0,0);
  char transA = 'N';
  char transB = 'N';
    cgemm_(&transA, &transB,
	   &M, &N, &K,
	   &alpha, A._data, &LDA, 
	   B._data, &LDB, &beta,
	   aMb._data,
	   &LDC);
    return aMb;
};

template<>
matrix<dcomplex> operator * (const matrix<dcomplex> &a, const matrix<dcomplex> &b){
  matrix<dcomplex> aMb(a._Nrow, b._Ncolumn, "("+a._name+"*"+b._name+")");
  matrix<dcomplex> A(a);
  matrix<dcomplex> B(b);
  int M = a._Nrow;
  int N = b._Ncolumn;
  int K = a._Ncolumn;
  int LDA = M;
  int LDB = b._Nrow;
  int LDC = aMb._Nrow;
  dcomplex alpha = dcomplex (1,0);
  dcomplex beta = dcomplex (0,0);
  char transA = 'N';
  char transB = 'N';
    zgemm_(&transA, &transB,
	   &M, &N, &K,
	   &alpha, A._data, &LDA, 
	   B._data, &LDB, &beta,
	   aMb._data,
	   &LDC);
    return aMb;
};

// operator |
template<>
matrix<fcomplex> operator | (const matrix<fcomplex> &a, const matrix<fcomplex> &b){
  int N, NRHS,LDA,LDB,INFO;
  int *IPIV;
  matrix<fcomplex> A(a);
  matrix<fcomplex> B(b);
  B._name = "("+a.name()+"|"+b.name()+")";
  N = a.row();
  NRHS = b.col();
  LDA = N;
  LDB = NRHS;
  IPIV = (int*)calloc(a.row(),sizeof(a.row()));
  cgesv_(&N, &NRHS, 
	 A._data, &LDA, IPIV, 
	 B._data, &LDB, &INFO);
  free(IPIV);
  return B;
};

template<>
matrix<dcomplex> operator | (const matrix<dcomplex> &a, const matrix<dcomplex> &b){
  int N, NRHS,LDA,LDB,INFO;
  int *IPIV;
  matrix<dcomplex> A(a);
  matrix<dcomplex> B(b);
  B._name = "("+a.name()+"|"+b.name()+")";
  N = a.row();
  NRHS = b.col();
  LDA = N;
  LDB = NRHS;
  IPIV = (int*)calloc(a.row(),sizeof(a.row()));
  zgesv_(&N, &NRHS, 
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


