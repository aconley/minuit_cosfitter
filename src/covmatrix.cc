//covmatrix.cc

#include <cmath>
#include <iostream>
#include <covmatrix.h>

#if USEMKL
#include <mkl_cblas.h>
#include <mkl_lapack.h>
#elif USEATLAS
extern "C" {
  #include <clapack.h>
  #include <cblas.h>
} 
#elif USEACCELERATE
#include <Accelerate/Accelerate.h>
#endif

//We avoid exceptions and most run-time error checking for 
//performance reasons, so be careful.

/*!
  \param[in] N Number of rows and columns
 */
void covMatrix::initialize(unsigned int N) {
  size_ = N*N;
  n_ = N;
  
  if ( N == 0 ) {
    v_ = NULL;
    row_ = NULL;
  } else {
    v_ = new double[size_]; 
    row_ = new double*[N];
    double* p = v_;              
    for (unsigned int i=0; i<N; i++) {
      row_[i] = p;
      p += N ;
    }      
  }
}

/*!
  \param[in] v Data to copy into internal array
*/
void covMatrix::copy(const double* v) {
  for (unsigned int i = 0; i < size_; ++i)
    v_[i] = v[i];
}

/*! 
  \param[in] val Sets all elements to this value
*/
void covMatrix::set(const double& val) {
  for (unsigned int i = 0; i < size_; ++i)
    v_[i] = val;
}

void covMatrix::destroy() {     
  /* do nothing, if no memory has been previously allocated */
  if (v_ == NULL) return ;
  
  /* if we are here, then matrix was previously allocated */
  if (v_ != NULL) delete [] (v_);     
  if (row_ != NULL) delete [] (row_);
  v_ = NULL;
  row_ = NULL;
  n_ = size_ = 0;
}

/*!
  This is internal because the user shouldn't be able to interface
  with the matrix when it's in lower triangular form.  Also, if you
  pass it something that isn't lower triangular, bad things will happen.
*/
void covMatrix::invertLowerTriangular() {
  double sum, *rowptr;
  for (unsigned int i = 0; i < n_; ++i) {
    rowptr = row_[i];
    //Do diagonal
    rowptr[i] = 1.0 / rowptr[i];
    for (unsigned int j = 0; j < i; ++j)
      row_[j][i] = 0.0;
    for (unsigned int j = i+1; j < n_; ++j) {
      rowptr = row_[j];
      sum = 0.0;
      for (unsigned int k = i; k < j; ++k)
	sum -= rowptr[k] * row_[k][i];
      rowptr[i] = sum / rowptr[j];
    }
  }
    
}

/*
  Factors the matrix into it's lower diagonal Cholesky form
  such that \f$A = L \cdot L^T\f$.  If the process fails,
  possibly because the matrix is not symmetric positive definite,
  then status will be non-zero on return.

  \param status 0 on success, -1 if non-spd
 */
void covMatrix::factorCholesky(int& status) {
  if (status) return;
  double d,s, *ptri, *ptrj;
  bool isspd;
  isspd = true;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    d = ptri[i];
    for (unsigned int j = 0; j < i; ++j)
      d -= ptri[j]*ptri[j];
    isspd = isspd && ( d > 0 );
    d = sqrt( d > 0.0 ? d : 0.0 );
    ptri[i] = d;
    for (unsigned int j = i+1; j < n_; ++j) {
      ptrj = row_[j];
      s = ptri[j];
      for (unsigned int k = 0; k < i; ++k)
	s -= ptri[k]*ptrj[k];
      ptrj[i] = ( d > 0 ? s/d : 0.0 );
    }
  }
  if ( isspd ) status = 0; else status = -1;
}
    
/*!
  Form the product of a lower triangular matrix and it's transpose,
  overwriting the current matrix.

  Note that the result is always symmetric.  This is internal
  because the user shouldn't interface with the lower triangular
  forms.
*/
void covMatrix::lowerTriangularTransposeMult() {

  double aii, sum, val;

  for (unsigned int i = 0; i < n_; ++i) {
    aii = v_[n_*i + i];
    if ( i < n_-1 ) {
      //Diagonal term
      sum = aii * aii;
      for ( unsigned int j = i+1; j < n_; ++j ) {
	val = v_[n_ * j + i];
	sum += val * val;
      }
      v_[ n_ * i + i ] = sum;

      //Do the lower part only, reflect later
      double *ptrk;
      for ( unsigned int j = 0; j < i; ++j) {
	//We have to handle the k=i case seperately,
	// since we already adjusted it's value
	ptrk = row_[i];
	sum = aii * ptrk[j];
	for ( unsigned int k = i+1; k < n_; ++k) {
	  ptrk = row_[k];
	  sum += ptrk[i] * ptrk[j];
	}
	v_[ n_ * i + j ] = sum;
      }
    } else {					
      //Bottom row, simple form
      double *rowptr;
      rowptr = row_[n_-1];
      for (unsigned int k = 0; k < n_; ++k)
	rowptr[k] *= aii;
    }
  }

  //Now reflect to upper half
  for (unsigned int i = 0; i < n_; ++i)
    for (unsigned int j = i+1; j < n_; ++j)
      v_[n_ * i + j] = v_[n_ * j + i];      
}

covMatrix::covMatrix() : n_(0), size_(0), v_(NULL), row_(NULL) { };

covMatrix::covMatrix(const covMatrix &A) {
  initialize(A.n_);
  copy(A.v_);
}

covMatrix::covMatrix(unsigned int M, double value) {
  initialize(M);
  set(value);
}

covMatrix::covMatrix(unsigned int N, double** arr) {
  initialize(N);
  unsigned int cntr;
  double *rowptr;
  for (unsigned int i = 0; i < N; ++i) {
    rowptr = arr[i];
    cntr = n_ * i;
    for (unsigned int j = 0; j < N; ++j) 
      v_[ cntr + j ] = rowptr[j];
  }
}

/*
  This operations occurs in place, i.e. when resizing to 
  a new matrix, original matrix elements
  are <b>NOT</b> retained.  Instead, one must explicit create
  a new matrix of this size and manually copy the elements, e.g.
  <pre>
  
  covMatrix B(N);
  
  unsigned int min_N = N < A.getNRows() ? N : A.getNRows();
  for (unsigned int i=0; i<=min_N; i++)
  for (unsigned int j=0; j<=min_N; j++)
  B[i][j] = A[i][j];
  
  A.destroy();
  </pre>
  
  \param[in] N new size
*/
covMatrix& covMatrix::resize(unsigned int N) {
  if ( n_ == N ) {
    return *this;
  }
  destroy();
  initialize(N);
  return *this;
}

covMatrix& covMatrix::resize(unsigned int N, double value) {
  if ( n_ == N ) {
    set(value);
    return *this;
  }
  destroy();
  initialize(N);
  set(value);
  return *this;
}

double covMatrix::getTotal() const {
  if (size_ == 0) return 0.0;
  double sum = v_[0];
  for (unsigned int i = 1; i < size_; ++i)
    sum += v_[i];
  return sum;
}

double covMatrix::getRowTotals( std::vector<double>& vec) const {
  if (size_ == 0) {
    vec.clear();
    return 0.0;
  }
  if (vec.size() != n_) vec.resize(n_);
  double sum, rowsum, *rowptr;
  sum = 0.0;
  for (unsigned int i = 0; i < n_; ++i) {
    rowptr = row_[i];
    rowsum = rowptr[0];
    for (unsigned int j = 1; j < n_; ++j)
      rowsum += rowptr[j];
    sum += rowsum;
    vec[i] = rowsum;
  }
  return sum;
}

covMatrix& covMatrix::operator=(const covMatrix& B) {
  //Self assignment check
  if (v_ == B.v_) return *this;
  if ( n_ == B.n_ ) //Don't need to realloc memory
    copy(B.v_);
  else {
    resize(B.n_);
    copy(B.v_);
  }
  return *this;
}

covMatrix& covMatrix::operator=(double val) {
  set(val);
  return *this;
}

std::vector<double> covMatrix::getDiag() const {
  std::vector<double> diag(n_);
  for (unsigned int i=0; i < n_; ++i)
    diag[i] = row_[i][i];
  return diag;
}

/*!
  Replaces the contents of the current matrix with
  \f$ A \cdot B \f$.

  \param[in] A One of the matricies
  \param[in] B The other matricix
  \param     status zero on successful calculation
*/
covMatrix& covMatrix::mult(const covMatrix&A, const covMatrix&B, int& status) {
  if ( status ) return *this;
  if ( v_ == A.v_ || v_ == B.v_ ) {
    //We can't handle this
    status = 2;
    return *this;
  }

  unsigned int n_A;
  n_A = A.getNRows();
  if ( n_A != B.getNRows() ) {
    status = 1;
    return *this;
  }
  if ( n_ != n_A ) resize( n_A );

#if USEMKL || USEATLAS || USEACCELERATE
  int in = static_cast<int>(n_);
  CBLAS_ORDER order = CblasRowMajor;
  CBLAS_TRANSPOSE trans = CblasNoTrans;
  cblas_dgemm(order,trans,trans,in,in,in,1.0,
	      A.v_,in,B.v_,in,0.0,v_,in );
#else
  double sum;
  const double* row_i;
  const double* col_k;
  for (unsigned int i = 0; i < n_; ++i) {
    row_i = A.row_[i];
    for (unsigned int k = 0; k < n_; ++k) {
      col_k = &(B[0][k]);
      sum = 0;
      for (unsigned int j=0; j< n_; ++j) {
	sum  += *row_i * *col_k;
	row_i++;
	col_k += n_;
      }
      row_[i][k] = sum; 
    }
  }
#endif

  return *this;
}

/*!
  Returns \f$ A \cdot \vec{b}\f$ where A is the current matrix

  \param[in] b The input vector
  \param     status zero on successful calculation

  The arithmetic is done in double precision internally.
*/
std::vector<float> covMatrix::mult( const std::vector<float> b,
				    int& status ) const {
  if ( status ) return std::vector<float>();
  unsigned int n_b;
  n_b = b.size();
  if ( n_b != n_ ) {
    status = 1;
    return std::vector<float>();
  }
  std::vector<float> tmp(n_);
  double sum, *rowptr;
  for (unsigned int i=0; i<n_; ++i) {
    sum = 0;
    rowptr = row_[i];
    for (unsigned int j=0; j<n_; ++j)
      sum +=  rowptr[j] * b[j];
    tmp[i] = static_cast<float>(sum); 
  }
  status = 0;
  return tmp;
}


/*!
  Returns \f$ A \cdot \vec{b}\f$, where A is the current matrix

  \param[in] b      The input vector
  \param     status zero on successful calculation, should be zero
                     on input.
*/
std::vector<double> covMatrix::mult( const std::vector<double> b,
				    int& status ) const {
  if ( status ) return std::vector<double>();
  unsigned int n_b;
  n_b = b.size();
  if ( n_b != n_ ) {
    status = 1;
    return std::vector<double>();
  }
  std::vector<double> tmp(n_);
  double sum, *rowptr;
  for (unsigned int i=0; i<n_; ++i) {
    sum = 0;
    rowptr = row_[i];
    for (unsigned int j=0; j<n_; ++j)
      sum +=  rowptr[j] * b[j];
    tmp[i] = sum; 
  }
  status = 0;
  return tmp;
}

covMatrix& covMatrix::scalarMultAndAdd( double scal, const covMatrix&A,
					int& status ) {
  if (status) return *this;
  if ( n_ != A.n_ ) {
    status = 1;
    return *this;
  }
  double *rowptr, *Arowptr;
  for (unsigned int i = 0; i < n_; ++i) {
    rowptr = row_[i];
    Arowptr = A.row_[i];
    for (unsigned int j = 0; j < n_; ++j)
      rowptr[j] += scal * Arowptr[j];
  }
  status = 0;
  return *this;
}

covMatrix& covMatrix::add(const covMatrix& A, int& status) {
  if (status) return *this;
  if ( v_ == A.v_ ) {
    //Self addition
    for (unsigned int i = 0; i < size_; ++i)
      v_[i] += v_[i];
    status = 0;
    return *this;
  }
  if ( A.getNRows() != n_ ) {
    status = 1;
    return *this;
  }
  for (unsigned int i = 0; i < size_; ++i)
    v_[i] += A.v_[i];
  status = 0;
  return *this;
}

covMatrix& covMatrix::operator*=(double val) {
  for (unsigned int i = 0; i < size_; ++i)
    v_[i] *= val;
  return *this;
}

/*!
  Invert using Cholesky decomposition
  \param status Returns 0 on success
*/
covMatrix& covMatrix::invert(int& status) {
  if (status) return *this;
  
#if USEMKL
  char uplo = 'U';
  int in = static_cast<int>(n_);
  dpotrf( &uplo, &in, v_, &in, &status );
  //std::cerr << "status from dpotrf: " << status << std::endl;
  if (status) return *this;
  dpotri( &uplo, &in, v_, &in, &status );
  //std::cerr << "status from dpotri: " << status << std::endl;
  //After this, only the lower half is correct.  Normally this
  // would be the upper half, since uplo='U', but CLAPACK
  // thinks we are in Column Major order
  //Anyways, reflect
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }
#elif USEATLAS
  ATLAS_ORDER Order = CblasRowMajor;
  ATLAS_UPLO Uplo = CblasLower;
  status = clapack_dpotrf( Order, Uplo, n_, v_, n_ );
  if (status) return *this;
  status = clapack_dpotri( Order, Uplo, n_, v_, n_ );
  if (status) return *this;
  //Only the lower half is correct -- reflect
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }
#elif USEACCELERATE
  char Uplo = 'U';
  __CLPK_integer N = static_cast<__CLPK_integer>(n_);
  __CLPK_integer STATUS = static_cast<__CLPK_integer>(status);
  dpotrf_(&Uplo, &N, v_, &N, &STATUS);
  status = static_cast<int>(STATUS);
  if (status) return *this;
  dpotri_(&Uplo, &N, v_, &N, &STATUS);
  status = static_cast<int>(STATUS);
  if (status) return *this;
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }  
#else
  status = 0;
  factorCholesky(status);
  if (status) return *this;
  invertLowerTriangular();
  lowerTriangularTransposeMult();
#endif
  return *this;
}

covMatrix& covMatrix::operator+=(const covMatrix& A) {
  int status;
  status = 0;
  return add( A, status );
}

std::ostream& operator<<(std::ostream &s, const covMatrix& A) {
  unsigned int N = A.getNRows();
  double const* rowptr;
  s << N << std::endl;
  for (unsigned int i=0; i<N; ++i) {
    rowptr = A[i];
    for (unsigned int j=0; j<N; ++j)
      s << rowptr[j] << " ";
    s << std::endl;
  }
  return s;
}

std::istream& operator>>(std::istream &s, covMatrix &A) {
  unsigned int N;
  s >> N;

  if ( !( N == A.getNRows() ) )
    A.resize(N);

  double *rowptr;
  for (unsigned int i=0; i<N; ++i) {
    rowptr = A[i];
    for (unsigned int j=0; j<N; ++j)
      s >> rowptr[j];
  }
  return s;
}

covMatrix operator*( const covMatrix&A, const covMatrix& B) {
  int status;
  status = 0;
  covMatrix tmp( A.getNRows() );
  return tmp.mult(A,B,status);
}

std::vector<double> operator*(const covMatrix&A, const std::vector<double>&b) {
  int status;
  status = 0;
  return A.mult( b, status );
}


