//covmatrix.h

#include <vector>
#include <fstream>

#ifndef __covmatrix__
#define __covmatrix__

/*!
  \defgroup linalg Linear Algebra
   
  \brief Code assosciated with linear algebra, primarily
    as it relates to covariance matricies.
*/

/*!
  \brief Code to handle covariance matricies
  \ingroup linalg

  More generally, this is specialized for handling
  double precision symmetric positive definite matricies.
  Based closely on the TNT3 beta code by Roldan Pozo.
  The code comes with a default inversion implementation using
  Cholesky decomposition (since covariance matricies are symmetric
  positive definite).  However, there are hooks during compilation
  for linking in LAPACK routines (ATLAS and MKL) which are highly
  recommended if you are thinking of using this code.

  Run-time error checking is mostly avoided for performance
  reasons, so be careful.  Some routines take a status variable,
  similar to the one used in cfitsio.  If a routine is passed
  a status variable that is non-zero, it will return immediately.
  If the routine fails somehow, it will set status to something
  that isn't 0.

*/
class covMatrix {
 private :
  unsigned int n_; //!< Number of rows and columns
  unsigned int size_; //!< Number of elements

  //The internal represenatation is as a single block rather than
  // the c-style because this makes everything allocate contiguously,
  // and because it is necessary for interfacing with external LAPACK
  // libraries
  //Note that the matrix is symmetric so row order versus column order
  // doesn't matter
  double *v_; //!< Holds elements
  double** row_;  //!< row pointers for efficient access
  
  /*! \brief Internal helper function for creating the array of row ptrs */
  void initialize(unsigned int M);

  /*! \brief Internal copying operator */
  void copy(const double* v);

  /*! \brief Sets all matrix elements to specified value */
  void set(const double&);

  /*! \brief Free memory of array */
  void destroy();

  //Routines having to do with the included inversion routine
  /*! \brief Factors matrix into Cholesky form in place */
  void factorCholesky(int& status);

  /*! \brief Inverts the matrix if it's in lower-triangular form*/
  void invertLowerTriangular();

  /*! \brief Multiply a lower triangular matrix by it's transpose */
  void lowerTriangularTransposeMult();

 public :

  /*! \brief Default constructor */
  covMatrix();
  
  /* ! \brief Copy constructor */
  covMatrix(const covMatrix&);

  /* ! \brief Constructor for given size (and value) */
  explicit covMatrix(unsigned int, double val = 0.0);

  /*! \brief Make from C array */
  explicit covMatrix(unsigned int, double**);

  /*! \brief Destructor */
  ~covMatrix() { destroy(); }

  /*! \brief Return the total number of elements */
  unsigned int getSize() const { return size_; }
  /*! \brief Return dimensions */
  unsigned int getNRows() const { return n_; }
  /*! \brief Return dimensions */
  unsigned int getNCols() const { return n_; }

  /*! \brief Change size of matrix, reallocating memory if necessary. */
  covMatrix& resize(unsigned int);
  /*! \brief Resize and set to value */
  covMatrix& resize(unsigned int,double);

  /*! \brief Clear out elements */
  void clear() { resize(0); }

  /*! \brief Single subscripting operator */
  inline double* operator[](unsigned int i) {
    return row_[i];
  }
  /*! \brief Single subscripting operator */
  inline double* operator[](int i) {
    return row_[i];
  }

  /*! \brief Single subscripting operator */
  inline double const* operator[](unsigned int i) const {
    return row_[i];
  }
  /*! \brief Single subscripting operator */
  inline double const* operator[](int i) const {
    return row_[i];
  }

  /*! \brief Get axis to raw data internals. Dangerous, but good for performance*/
  inline double* getData() { return v_; }
  /*! \brief Get axis to raw data internals. Dangerous, but good for performance*/
  inline double const* getData() const { return v_; }

  double** operator()() { return row_; } //!< Access to row pointers
  double const* const* operator()() const { return row_; } //!< Const access to row pointers

  /*! \brief Get the total of all the elements */
  double getTotal() const;
  /*! \brief Get the sum of the rows, also returning the total */
  double getRowTotals( std::vector<double>& ) const;

  covMatrix& operator=(const covMatrix& B); //!< Copy one matrix to another
  covMatrix& operator=(double val); //!< Set all elements to value

  /*! \brief Get the diagonal elements */
  std::vector<double> getDiag() const;

  /*! \brief Matrix/Matrix multiplication: C = A*B */
  covMatrix& mult(const covMatrix&A, const covMatrix&B, int& status);

  /*! \brief Matrix/Vector multiplication: \f$A \cdot \vec{b}\f$ */
  std::vector<float> mult( const std::vector<float> b, int& status ) const;

  /*! \brief Matrix/Vector multiplication: \f$A \cdot \vec{b}\f$ */
  std::vector<double> mult( const std::vector<double> b, int& status ) const;

  /*! \brief Add a scalar times another matrix to this one */
  covMatrix& scalarMultAndAdd( double scal, const covMatrix& A, int& status );

  /*! \brief Matrix addition */
  covMatrix& add(const covMatrix&A, int& status );

  /*! \brief Matrix scaling */
  covMatrix& operator*=(double);

  /*! \brief Matrix addition */
  covMatrix& operator+=(const covMatrix&);

  /*! \brief Invert in place*/
  covMatrix& invert(int& status);
};

/*! \brief Output operator */
std::ostream& operator<<(std::ostream &s, const covMatrix& A);
/*! \brief Input operator */
std::istream& operator>>(std::istream &s, covMatrix& A);

//Linear algebra
/*! \brief Matrix-Matrix multiplication: C = A*B */
covMatrix operator*( const covMatrix &A, const covMatrix& B);
/*! \brief Matrix-vector multiplication */
std::vector<double> operator*(const covMatrix& A, 
			      const std::vector<double>& b);
#endif
