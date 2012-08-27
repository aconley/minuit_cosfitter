//simple_matrix.h

#ifndef __simple_matrix__
#define __simple_matrix__

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <assert.h>

#include <cosfitterexcept.h>


/*!
  \brief Code to handle Matricies

  The code is largely copied over from the TNT 3 beta,
  which is distributed without licensing restrictions.
  The author of TNT3 is Roldan Pozo.

*/
namespace simple_matrix {

  /*!
    \brief Matrix class
    
    It is stored in row major order, index starts at [0][0], and has
    copy by value semantics.
  */
  
  template <class T> class Matrix {
    
  private:
    unsigned int m_;       //!< Number of rows
    unsigned int n_;       //!< Number of columns
    unsigned int mn_;      //!< total size
    T* v_;        //!< holds elements                  
    T** row_;     //!< row pointers for efficient access
    
    /*!
      \brief internal helper function to create the array
      of row pointers
    */
    void initialize(unsigned int M, unsigned int N) {
      mn_ = M*N;
      m_ = M;
      n_ = N;
      
      v_ = new T[mn_]; 
      row_ = new T*[M];
      
      assert(v_  != NULL);
      assert(row_  != NULL);
      
      T* p = v_;              
      for (unsigned int i=0; i<M; i++) {
	row_[i] = p;
	p += N ;
	
      }      
    }

    /*!
      \brief Internal copying operator
    */
    void copy(const T* v) {
      for (unsigned int i = 0; i < mn_; ++i)
	v_[i] = v[i];
    }

    /*!
      \brief Sets all matrix entries to specified value
    */
    void set(const T& val) {
      for (unsigned int i=0; i< mn_; i++)
	v_[i] = val;
    }

    /*!
      \brief Free memory of array.
    */
    void destroy() {     
      /* do nothing, if no memory has been previously allocated */
      if (v_ == NULL) return ;
      
      /* if we are here, then matrix was previously allocated */
      if (v_ != NULL) delete [] (v_);     
      if (row_ != NULL) delete [] (row_);
      
    }
    
  public:
    
    int lbound() const { return 0;} //!< Return lower bound of indexing
    
    operator T**(){ return  row_; } //!< Gives access to row pointers
    operator const T**() const { return row_; } //!< Access to row pointers
    
    
    /*!
      \return the total number of items in matrix (M*N).
    */
    unsigned int size() const { return mn_; }
    
    /*!
      \brief Default constructor
    */
    Matrix() : m_(0), n_(0), mn_(0), v_(0), row_(0) {};
    
    /*!
      \brief Copy constructor
    */
    Matrix(const Matrix<T> &A) {
      initialize(A.m_, A.n_);
      copy(A.v_);
    }
    
    /*!
      \brief Create a MxN matrix, with each element assigned to the value 0.
      
      \param[in] M the number of rows
      \param[in] N the number of columns
      \param[in] value optional default value: 0 if not specified.
      
    */                                              
    Matrix(unsigned int M, unsigned int N, const T& value = T(0)) {
      initialize(M,N);
      set(value);
    }
    
    /*!
      \brief Make from c array

      Create an MxN matrix, filling in values (row-major order) from
      the list (C array) provided.
      
      \param[in] M the number of rows
      \param[in] N the number of columns
      \param[in] v  list (C array) of M*N values used to initialize matrix.
    */
    Matrix(unsigned int M, unsigned int N, const T* v) {
      initialize(M,N);
      copy(v);
    }
    
    /*!
      \brief Create matrix from string

      Create an MxN matrix, filling in values (row-major order) from
      a character string.
      
      \param[in] M the number of rows
      \param[in] N the number of columns
      \param[in] s  string of M*N values used to initialize matrix.
    */
    Matrix(unsigned int M, unsigned int N, const char *s) {
      initialize(M,N);
      std::istringstream ins(s);
      
      unsigned int i, j;
      
      for (i=0; i<M; i++)
	for (j=0; j<N; j++)
	  ins >> row_[i][j];
    }
    
    /*!
      \brief destructor
    */
    ~Matrix() {
      destroy();
    }
    
    /*!
      \brief Change size of matrix to MxN, reallocating memory if necessary.
      
      This operations occurs in place, i.e. when resizing to 
      a new matrix, original matrix elements
      are <b>NOT</b> retained.  Instead, one must explicit create
      a new matrix of this size and manually copy the elements, e.g.
      <pre>
      
      Matrix<double> B(M, N);
      
      unsigned int min_M = M < A.num_rows() ? M : A.num_rows();
      unsigned int min_N = N < A.num_cols() ? N : A.num_cols();
      for (unsigned int i=1; i<=min_M; i++)
      for (unsigned int j=1; j<=min_N; j++)
      B(i,j) = A(i,j);
      
      A.destroy();
      </pre>
      
      \param[in] M the number of rows of new size.
      \param[in] N the number of columns of new size.
    */
    Matrix<T>& resize(unsigned int M, unsigned int N) {
      if (num_rows() == M && num_cols() == N)
	return *this;
      
      destroy();
      initialize(M,N);
      
      return *this;
    }
    
    /*!
      Assign (copy) one matrix to another, e.g. A=B.  The
      contents of A are lost, and a new copy of B is created.  
      
      \param[in] B to matrix to be copied.
      
    */
    Matrix<T>& operator=(const Matrix<T> &B) {
      //Self assignment check
      if (v_ == B.v_)
	return *this;
      
      if (m_ == B.m_ && n_ == B.n_)      // no need to re-alloc
	copy(B.v_);
      else {
	destroy();
	initialize(B.m_, B.n_);
	copy(B.v_);
      }
      
      return *this;
    }
 
    /*!
      \brief Set all entries of matrix to scalar value
    */
    Matrix<T>& operator=(const T& scalar) {
      set(scalar); 
      return *this;
    }
    
    /*!
      \brief Get the extent along the specified dimension
      
      \param[in] d The dimension to get the extent along --
      either 1 or 2
    */
    unsigned int dim(unsigned int d) const {
      return (d==1) ? m_ : ((d==2) ? n_ : 0); 
    }
    
    unsigned int num_rows() const { return m_; } //!< Give the number of rows
    unsigned int num_cols() const { return n_; } //!< Give the number of columns
    
    /*! \brief Single subscripting operator */
    inline T* operator[](unsigned int i) {
      return row_[i];
    }
    
    /*! \brief const single subscripting operator */
    inline const T* operator[](unsigned int i) const {
      return row_[i];
    }
    
    /*! \brief Return the diagonal elements */
    std::vector<T> diag() const {
      unsigned int N = n_ < m_ ? n_ : m_;  //handle non-diagonal case
      std::vector<T> d(N);
      
      for (unsigned int i=0; i<N; i++)
	d[i] = row_[i][i];
      
      return d;
    }
    
  };
  
  //Input and output
  
  /*!
    \brief Output operator
  */
  template <class T>
    std::ostream& operator<<(std::ostream &s, const Matrix<T> &A) {

    unsigned int M=A.num_rows();
    unsigned int N=A.num_cols();
    
    for (unsigned int i=0; i<M; i++) {
      for (unsigned int j=0; j<N; j++)
	s << A[i][j] << " ";
      s << "\n";
    }
    
    return s;
  }
  
  
  /*!
    \brief Input operator
  */
  template <class T>
    std::istream& operator>>(std::istream &s, Matrix<T> &A)
    {
      
      unsigned int M, N;
      
      s >> M >> N;
      
      if ( !(M == A.num_rows() && N == A.num_cols() ))
	{
	  A.resize(M,N);
	}
      
      
      for (unsigned int i=0; i<M; i++)
	for (unsigned int j=0; j<N; j++)
	  s >>  A[i][j];
      
      return s;
    }
    
  /*!
    \brief Matrix-Matrix multiplication:  C = A * B.
    
    This is an optimizied (trinary) version of matrix multiply, where 
    the destination matrix has already been allocated.
    
    \param[in] A	matrix of size M x N.
    \param[in] B	matrix of size N x K.
    \param[out] C  the result A*B, of size M x K.  Must have been
    preallocated already.
    
    \return a reference to C, after multiplication.
  */
  template <class T>
    Matrix<T> & mult(Matrix<T>& C, const Matrix<T>  &A, 
		     const Matrix<T> &B) {
    
    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    unsigned int K = B.num_cols();
    
    if ( C.num_rows() != M )
      throw CosFitterExcept("matrix","mult",
			    "C array not right size",1);
    if ( C.num_cols() != K )
      throw CosFitterExcept("matrix","mult",
			    "C array not right size",1);
    if ( B.num_rows != N ) 
      throw CosFitterExcept("matrix","mult",
			    "Array sizes for mult not correct",2);
    
    T sum;
    
    const T* row_i;
    const T* col_k;
    
    for (unsigned int i=0; i<M; i++)
      for (unsigned int k=0; k<K; k++) {
	row_i  = &(A[i][0]);
	col_k  = &(B[0][k]);
	sum = 0;
	for (unsigned int j=0; j<N; j++) {
	  sum  += *row_i * *col_k;
	  row_i++;
	  col_k += K;
	}
	C[i][k] = sum; 
      }
    
    return C;
  }
  
  
  /*!
    \brief Matrix/matrix multiplication: compute A * B.
    
    \param[in] A matrix: left side operand  (size M x N).
    \param[in] B matrix: right side operand  (size N x K).
    
    \return A*B  (a new matirx of size M x K).
  */
  template <class T>
    inline Matrix<T> mult(const Matrix<T>  &A, const Matrix<T> &B) {
    
    
    if ( A.num_cols() != B.num_rows() )
      throw CosFitterExcept("matrix","mult",
			    "Array sizes for mult not correct",1);
    
    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    unsigned int K = B.num_cols();
    
    Matrix<T> tmp(M,K);
    
    mult(tmp, A, B);		// tmp = A*B
    
    return tmp;
  }
  
  /*!
    \brief Matrix/matrix multiplication: compute A * B.
    
    \param[in] A matrix: left side operand  (size M x N).
    \param[in] B matrix: right side operand  (size N x K).
    
    \return A*B  (a new matrix of size M x K).
  */
  template <class T>
    inline Matrix<T> operator*(const Matrix<T>  &A, const Matrix<T> &B) {
    return mult(A,B);
  }
  
  /*!
    \brief Matrix/vector multiplication: compute A * b.
    
    \param[in] A matrix: left side operand  (number of columns of A, 
       must match the number of elements in b.) size MxN
    \param[in] b vector: right side operand.
    \return \f$ A \cdot \vec{b}\f$ (a new vector of size M.)
  */
  
  template <class T>
    inline std::vector<T> mult(const Matrix<T>  &A, const std::vector<T> &b) {
    
    if ( A.num_cols() != b.size() )
      throw CosFitterExcept("matrix","mult",
			    "Array/vector sizes for mult not correct",1);
    
    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    
    std::vector<T> tmp(M);
    T sum;
    
    for (unsigned int i=0; i<M; i++) {
      sum = 0;
      for (unsigned int j=0; j<N; j++)
	sum +=  A[i][j] * b[j];
      
      tmp[i] = sum; 
    }
    
    return tmp;
  }
  
  /*!
    \brief Matrix/vector multiplication: compute A * b.
    
    \param[in] A matrix: left side operand  (number of columns of A must match 
    the number of elements in b.) Size MxN
    \param[in] b vector: right side operand.
    \return \f$ A \cdot \vec{b}\f$ (a new vector of size M.)
  */
  
  template <class T>
    inline std::vector<T> operator*(const Matrix<T>  &A, 
				    const std::vector<T> &b) {
    return mult(A,b);
  }
  
  /*!
    \brief Matrix scaling: multiply each element of A by scalar s.
    
    This creates a new copy of A.  To scale "in place",
    use *= or mult_eq().
    
    \param[in] A matrix: to be scaled. 
    \param[in] s scalar: multiplier.
    \return s*A, a new matrix with same size of A.
  */
  template <class T>
    inline Matrix<T> mult(const T& s, const Matrix<T> &A) {

    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    
    Matrix<T> R(M,N);
    for (unsigned int i=0; i<M; i++)
      for (unsigned int j=0; j<N; j++)
	R[i][j] = s * A[i][j];
    
    return R;
  }

  /*!
    \brief Matrix scaling: multiply each element of A by scalar s.
    
    Same as mult(A,s), as this is a commutative operation.
    
    This creates a new copy of A.  To scale "in place",
    use *= or mult_eq().
    
    \param[in] A matrix: to be scaled. 
    \param[in] s scalar: multiplier.
    \return s*A, a new matrix with same size of A.
  */
  template <class T>
    inline Matrix<T> mult(const Matrix<T> &A, const T& s) {
    
    return mult(s, A);
  }
  
  /*!
    \brief Matrix scale in place
    
    Matrix scale in-place, i.e. compute A *= s, where each element
    of A is multiplied (scaled) by the value s.
    
    \param[in] A matrix: to be scaled. 
    \param[in] s scalar: multiplier.
    \return A, after scaling.
  */
  template <class T>
    inline Matrix<T> mult_eq(const T& s, Matrix<T> &A) {

    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    
    for (unsigned int i=0; i<M; i++)
      for (unsigned int j=0; j<N; j++)
	A[i][j] *= s;
    
    return A;
  }

  /*!
    \brief Matrix scale in place.
  */
  template <class T>
    inline Matrix<T> mult_eq(Matrix<T> &A, const T&a) {
    return mult_eq(a, A);
  }
  
  /*!
    \brief Matrix-Matrix tranpose multiplication, i.e. compute tranpose(A)*B.
    
    This is more efficient than computing the tranpose(A) explicitly,
    and then multiplying, as the tranpose of A is never really constructed.
    
    \param[in] A  matrix: size M x N.
    \param[in] B	matrix: size M x K.
    \return a new matrix of size N x K equal to \f$ A^T \cdot B \f$
  */
  template <class T>
    inline Matrix<T> transpose_mult(const Matrix<T> &A, 
				    const Matrix<T> &B) {
    
    if ( A.num_cols() != B.num_rows() )
      throw CosFitterExcept("matrix","transpose_mult",
			    "Array sizes for mult not correct",1);
    
    unsigned int M = A.num_cols();
    unsigned int N = A.num_rows();
    unsigned int K = B.num_cols();
    
    Matrix<T> tmp(M,K);
    T sum;
    
    for (unsigned int i=0; i<N; i++)
      for (unsigned int k=0; k<K; k++) {
	sum = 0;
	for (unsigned int j=0; j<M; j++)
	  sum = sum +  A[j][i] * B[j][k];
	
	tmp[i][k] = sum; 
      }
    
    return tmp;
  }
  
  /*!
    \brief Matrix-Vector tranpose multiplication, i.e. compute tranpose(A)*b.
    
    This is more efficient than computing the tranpose(A) explicitly,
    and then multiplying, as the tranpose of A is not explicitly constructed.
    
    \param[in] A  Matrix: size M x N.
    \param[in] b	vector: size M.
    \return a new vector of size N equal to \f$ A^T \cdot \vec{b} \f$.
  */
  template <class T>
    inline std::vector<T> transpose_mult(const Matrix<T>  &A, 
					 const std::vector<T> &b) {
    
    if ( A.num_cols() != b.size() )
      throw CosFitterExcept("matrix","transpose_mult",
			    "Sizes for array/vec transpose mult not correct",1);
    
    unsigned int M = A.num_cols();
    unsigned int N = A.num_rows();
    
    std::vector<T> tmp(M);
    
    for (unsigned int i=0; i<M; i++) {
      T sum = 0;
      for (unsigned int j=0; j<N; j++)
	sum +=  A[j][i] * b[j];
      
      tmp[i] = sum; 
    }
    
    return tmp;
  }
  
  
  /*!
    \brief Matrix addition: compute A + B
    
    \param[in] A	matrix of size M x N.
    \param[out] B	matrix of size M x N.
    
    \return the sum A+B.
  */
  template <class T>
    Matrix<T> add(const Matrix<T> &A, const Matrix<T> &B) {
    
    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    
    if ( B.num_cols() != N )
      throw CosFitterExcept("matrix","add",
			    "Sizes for array add not correct",1);
    
    if ( B.num_rows() != M )
      throw CosFitterExcept("matrix","add",
			    "Sizes for array add not correct",1);
    
    
    Matrix<T> tmp(M,N);
    unsigned int i,j;
    
    for (i=0; i<M; i++)
      for (j=0; j<N; j++)
	tmp[i][j] = A[i][j] + B[i][j];
    
    return tmp;
  }
  
  /*!
    Matrix addition: compute A + B.
    
    
    \param[in] A	matrix of size M x N.
    \param[in] B	matrix of size M x N.
    
    \returns the sum A+B.
  */
  template <class T>
    inline Matrix<T> operator+(const Matrix<T> &A, const Matrix<T> &B) {
    return add(A,B);
  }
  
  
  /*!
    \brief Matrix subtraction : compute A - B.
    
    \param[in] A	matrix of size M x N.
    \param[in] B	matrix of size M x N.
    
    \return the result A-B.
  */
  template <class T>
    Matrix<T> minus(const Matrix <T>& A, const Matrix<T> &B) {
    
    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    
    if ( B.num_cols() != N )
      throw CosFitterExcept("matrix","minus",
			    "Sizes for array sub not correct",1);
    
    if ( B.num_rows() != M )
      throw CosFitterExcept("matrix","minus",
			    "Sizes for array sub not correct",1);
    
    Matrix<T> tmp(M,N);
    unsigned int i,j;
    
    for (i=0; i<M; i++)
      for (j=0; j<N; j++)
	tmp[i][j] = A[i][j] - B[i][j];
    
    return tmp;
    
  }
  
  /*!
    \brief Matrix subtraction : compute A - B.

    \param[in] A	matrix of size M x N.
    \param[in] B	matrix of size M x N.
    
    \return the result A-B.
  */
  template <class T>
    inline Matrix<T> operator-(const Matrix<T> &A, 
			       const Matrix<T> &B) {
    return minus(A,B);
  }
  
  
  /*!
    \brief Element by element multiplication
    
    Matrix element-by-elment multiplication: for each (i,j)
    compute A(i,j) * B(i,j).
    
    \param[in] A matrix of size M x N.
    \param[in] B matrix of size M x N.
    \return new matrix, where each (i,j) is A(i,j) * B(i,j);
    
  */
  template <class T>
    Matrix<T> mult_element(const Matrix<T> &A, const Matrix<T> &B) {

    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    
    if ( B.num_cols() != N )
      throw CosFitterExcept("matrix","mult_element",
			    "Sizes for arrays not correct",1);
    
    if ( B.num_rows() != M )
      throw CosFitterExcept("matrix","mult_element",
			    "Sizes for arrays not correct",1);
    
    Matrix<T> tmp(M,N);
    unsigned int i,j;
    
    for (i=0; i<M; i++)
      for (j=0; j<N; j++)
	tmp[i][j] = A[i][j] * B[i][j];
    
    return tmp;
  }
  
  /*!
    \brief In place element by element
    
    Matrix element-by-elment multiplication, in place: for each (i,j)
    compute A(i,j) = A(i,j) * B(i,j).
    
    \param A matrix of size M x N, modified on output
    \param[in] B matrix of size M x N.
  */
  template <class T>
    void mult_element_eq(const Matrix<T> &A, const Matrix<T> &B) {
    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    
    if ( B.num_cols() != N )
      throw CosFitterExcept("matrix","mult_element_eq",
			    "Sizes for arrays not correct",1);
    
    if ( B.num_rows() != M )
      throw CosFitterExcept("matrix","mult_element_eq",
			    "Sizes for arrays not correct",1);
    
    
    unsigned int i,j;
    
    for (i=0; i<M; i++)
      for (j=0; j<N; j++)
	A[i][j] *= B[i][j];
    
  }
  
  
  /*!
    \brief Compute Frobenius norm of matrix.  
    
    This is the square root of the sum of squares of each matrix entry, i.e.
    \f$
      \sqrt{ \sum_{i=1}^{N} \sum_{j=1}^{N} A_{ij}^2 }.
    \f$
    
    \param[in] A the matrix to compute its Frobeinus norm.
    \return the Frobenius norm of A.
  */
  template <class T>
    T norm(const Matrix<T> &A) {

    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
    
    T sum = 0.0;
    for (unsigned int i=1; i<=M; i++)
      for (unsigned int j=1; j<=N; j++)
	sum += A(i,j) * A(i,j);
    return sqrt(sum);
  }
  
  /*!
    
  \brief Matrix transpose
  
  
  \param[in] A matrix MxN
  \return new matrix of size N x M, where each (i,j) is A(j,i).
  */
  template <class T>
    Matrix<T> transpose(const Matrix<T> &A) {

    unsigned int M = A.num_rows();
    unsigned int N = A.num_cols();
  
    Matrix<T> S(N,M);
    unsigned int i, j;
    
    for (i=0; i<M; i++)
      for (j=0; j<N; j++)
	S[j][i] = A[i][j];
    
    return S;
  }

  /*!
    \brief Invert the lower triangular matrix L
    
    Invert the lower triangular part of the matrix, assuming
    the upper part is all zero.
    \return The inverse of the lower triagonal portion
  */
  template <class T>
    Matrix<T> invert_lower_triangular(const Matrix<T>& L) {
    
    unsigned int N = L.num_rows();
    
    Matrix<T> tmp(L);
    T sum;
    for (unsigned int i = 0; i < N; ++i) {
      tmp[i][i] = 1.0 / L[i][i];
      for (unsigned int j = 0; j < i; ++j)
	tmp[j][i] = static_cast<T>(0);
      for (unsigned int j = i+1; j < N; ++j) {
	sum = static_cast<T>(0);
	for (unsigned int k = i; k < j; ++k)
	  sum -= tmp[j][k]*tmp[k][i];
	tmp[j][i] = sum/L[j][j];
      }
    }
    
    return tmp;
  }

  /*!
    \brief Form the product L^T * L

    For the product of a lower triangular matrix and it's transpose.
    Note that the result is always symmetric.

    If the input matrix is not lower triangular, the upper portion
    is simply ignored.
    
    \return The product \f$L^T \cdot L\f$
  */
  template <class T>
    Matrix<T> lower_transpose_mult(const Matrix<T>& L) {
    
    unsigned int N = L.num_rows();
    
    if (L.num_cols() != N) 
      throw CosFitterExcept("matrix","lower_transpose_mult",
			    "Array not symmetric",1);

    Matrix<T> tmp(N,N,static_cast<T>(0));
    T sum;

    //Do the diagonal
    for (unsigned int i = 0; i < N; ++i) {
      sum = static_cast<T>(0);
      for (unsigned int j = i; j < N; ++j)
	sum += L[j][i] * L[j][i];
      tmp[i][i] = sum;
    }

    //Then we do the upper half
    for (unsigned int i = 0; i < N; ++i) 
      for (unsigned int j = i+1; j < N; ++j) {
	sum = static_cast<T>(0);
	for (unsigned int k = j; k < N; ++k)
	  sum += L[k][j] * L[k][i];
	tmp[i][j] = sum;
      }
    //And reflect
    for (unsigned int i = 0; i < N; ++i)
      for (unsigned int j = 0; j < i; ++j)
	tmp[i][j] = tmp[j][i];

    return tmp;
  }
  
}
#endif
