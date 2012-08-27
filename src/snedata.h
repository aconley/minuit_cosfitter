#ifndef __snedata__
#define __snedata__

#include <vector>
#include <string>
#include <iterator>
#include <set>

#include <covmatrix.h>

/*!
  \defgroup snedata Supernova Data
  \brief Collection of information about supernovae

  These classes provide a mechanism for holding lists of SN parameters --
  that is, the data the cosmological parameters are fit to.
*/

/*!
  \brief Data structure specialized for holding supernova data.

  \ingroup snedata
*/
struct SNeDataEntry {
  std::string name; //!< Name of supernova
  float zcmb; //!< Redshift of supernova in CMB frame
  float zhel; //!< Heliocentric redshift of supernova
  float var_z; //!< Redshift variance
  float widthpar; //!< Width-luminosity parameter (stretch, delta m_15, etc.)
  float var_widthpar; //!< Variance of widthpar
  float mag; //!< Magnitude
  float var_mag; //!< Variance of magnitude
  float colourpar; //!< Colour related parameter ( E(B-V), A_V, etc. )
  float var_colourpar; //!< Variance of colourpar
  float cov_mag_widthpar; //!< Covariance between magnitude and widthpar
  float cov_mag_colourpar; //!< Covariance between magnitude and colourpar
  float cov_widthpar_colourpar; //!< Covariance between widthpar and colourpar
  unsigned int dataset; //!< Subset identifier

  /*! \brief Constructor */
  SNeDataEntry(std::string SNNAME="",float ZCMB=0,float ZHEL=0,float VARZ=0,
	       float WIDTHPAR=0, float VARWIDTHPAR=0,float MAG=0,
	       float VARMAG=0, float COLOURPAR=0, float VARCOLOURPAR=0, 
	       float COV_MAG_WIDTHPAR=0, float COV_MAG_COLOURPAR=0, 
	       float COV_WIDTHPAR_COLOURPAR=0, unsigned int DATASET=0) :
    name(SNNAME),zcmb(ZCMB),zhel(ZHEL),var_z(VARZ),widthpar(WIDTHPAR),
    var_widthpar(VARWIDTHPAR),mag(MAG),var_mag(VARMAG),colourpar(COLOURPAR),
    var_colourpar(VARCOLOURPAR),cov_mag_widthpar(COV_MAG_WIDTHPAR),
    cov_mag_colourpar(COV_MAG_COLOURPAR),
    cov_widthpar_colourpar(COV_WIDTHPAR_COLOURPAR),
    dataset(DATASET)
  { };

  /*! \brief Copy Constructor */
  SNeDataEntry(const SNeDataEntry& inval) {
    name = inval.name; zcmb = inval.zcmb; zhel = inval.zhel; 
    var_z = inval.var_z; widthpar = inval.widthpar; 
    var_widthpar = inval.var_widthpar; mag = inval.mag;
    var_mag = inval.var_mag; colourpar = inval.colourpar; 
    var_colourpar = inval.var_colourpar; 
    cov_mag_widthpar = inval.cov_mag_widthpar;
    cov_mag_colourpar = inval.cov_mag_colourpar;
    cov_widthpar_colourpar = inval.cov_widthpar_colourpar;
    dataset = inval.dataset;
  }

  /*! \brief Returns corrected magnitude */
  float getCorrectedMag(float,float) const;

  /*! \brief Returns error on corrected magnitude */
  float getCorrectedMagError(float,float) const;

  /*! \brief Returns variance of corrected magnitude */
  float getCorrectedMagVar(float,float) const;

  unsigned int getDataSet() const { return dataset; } //!< Returns dataset number

  bool isCovValid() const; //!< Check if covariance info is valid

  /*! \brief Copy operator */
  SNeDataEntry& operator=(const SNeDataEntry& inval) {
    if (this == &inval) return *this; //Self copy protection
    name = inval.name; zcmb = inval.zcmb; zhel = inval.zhel; 
    var_z = inval.var_z; widthpar = inval.widthpar; 
    var_widthpar = inval.var_widthpar; mag = inval.mag;
    var_mag = inval.var_mag; colourpar = inval.colourpar; 
    var_colourpar = inval.var_colourpar; 
    cov_mag_widthpar = inval.cov_mag_widthpar;
    cov_mag_colourpar = inval.cov_mag_colourpar;
    cov_widthpar_colourpar = inval.cov_widthpar_colourpar;
    dataset = inval.dataset;
    return *this;
  }

};


/*!
  \brief Sorting class for SNeDataEntry based on zcmb.
  \ingroup snedata
*/
class SNeSortByZcmb {
 public:
  /*! \brief Comparison operator on zcmb */
  int operator()(const SNeDataEntry& sn1, const SNeDataEntry& sn2) {
    return sn1.zcmb < sn2.zcmb;
  }
};


/*!
  \brief Collection of SNeDataEntry objects

  \ingroup snedata

  All data points are stored in the vector<SNeDataEntry> points, which you can
  access through the [] operator.  This class also provides serialization
  capability (i.e., it can write the data to a file and read it back).
*/
class SNeData {
 private:
  std::string Name; //!< Just a handy name, potentially useful for debugging
  std::vector<SNeDataEntry> points; //!< The data as SNeDataEntry's

  bool isdiagonal; //!< True if the errors are diagonal

  bool havemagcov;  //!< Do we have the magnitude covariance matrix?
  bool havewidthcov; //!< Do we have the width covariance matrix?
  bool havecolourcov; //!< Do we have the colour covariance matrix?
  bool havemagwidthcov; //!< Do we have the mag-width covariance matrix?
  bool havemagcolourcov; //!< Do we have the mag-colour covariance matrix?
  bool havewidthcolourcov; //!< Do we have the colour-width covariance matrix?

  covMatrix mag_covmatrix; //!< Stores the magnitude covariance matrix
  covMatrix width_covmatrix; //!< Stores the width covariance matrix
  covMatrix colour_covmatrix; //!< Stores the colour covariance matrix
  covMatrix magwidth_covmatrix; //!< Stores the mag-width covariance matrix
  covMatrix magcolour_covmatrix; //!< Stores the mag-colour covariance matrix
  covMatrix widthcolour_covmatrix; //!< Stores the width-colour covariance matrix

  void readCovData(const std::string& FileName,
		   covMatrix& covmat); //!< Reading function for cov matricies
  void writeCovData(const std::string& FileName,
		    const covMatrix& covmat) const; //!< Writing function for cov matricies


 public:
  
  typedef std::vector<SNeDataEntry>::iterator iterator;  //!< Forward iterator
  typedef std::vector<SNeDataEntry>::const_iterator const_iterator; //!< Forward iterator on constant
  typedef std::vector<SNeDataEntry>::reverse_iterator reverse_iterator; //!< Reverse iterator
  typedef std::vector<SNeDataEntry>::const_reverse_iterator const_reverse_iterator; //!< Reverse iterator on constant

  SNeData(); //!< Default constructor
  explicit SNeData(const std::string& name); //!< No data read
  SNeData(const std::string& name, const std::string& fileName); //!< Create SNeData object and read in data from file 

  void readData(const std::string&, bool verbose=false); //!< Read in the data from file 
  void readMagCovData(const std::string&); //!< Read in the mag covariance matrix from a file 
  void readWidthCovData(const std::string&); //!< Read in the width covariance matrix from a file
  void readColourCovData(const std::string&); //!< Read in the colour covariance matrix from a file
  void readMagWidthCovData(const std::string&); //!< Read in the mag-width covariance matrix from a file
  void readMagColourCovData(const std::string&); //!< Read in the mag-colour covariance matrix from a file
  void readWidthColourCovData(const std::string&); //!< Read in the width-colour covariance matrix from a file

  void writeData(const std::string&) const; //!< Write the data to a file

  void writeMagCovData(const std::string&) const; //!< Write the mag covariance matrix to a file
  void writeWidthCovData(const std::string&) const; //!< Write the width covariance matrix to a file
  void writeColourCovData(const std::string&) const; //!< Write the colour covariance matrix to a file
  void writeMagWidthCovData(const std::string&) const; //!< Write the mag-width covariance matrix to a file
  void writeMagColourCovData(const std::string&) const; //!< Write the mag-colour covariance matrix to a file
  void writeWidthColourCovData(const std::string&) const; //!< Write the width-colour covariance matrix to a file


  void setName(const std::string& name) { Name = name; } //!< Set name of object
  std::string getName() const { return Name; } //!< Return the name of the object

  bool areErrorsDiagonal() const { return isdiagonal; }  //!< Returns true if there are no covariance terms
  bool haveCovMatrix() const { return havemagcov || havewidthcov ||
      havecolourcov || havemagwidthcov || havemagcolourcov ||
      havewidthcolourcov; } //!< Do we have a full style covariance matrix
  bool haveMagCovMatrix() const { return havemagcov; } //!< Do we have a mag covariance matrix?
  bool haveWidthCovMatrix() const { return havewidthcov; } //!< Do we have a width covariance matrix
  bool haveColourCovMatrix() const { return havecolourcov; } //!< Do we have a colour covariance matrix
  bool haveMagWidthCovMatrix() const { return havemagwidthcov; } //!< Do we have a mag-width covariance matrix
  bool haveMagColourCovMatrix() const { return havemagcolourcov; } //!< Do we have a mag-colour covariance matrix
  bool haveWidthColourCovMatrix() const { return havewidthcolourcov; } //!< Do we have a width-colour covariance matrix

  //Getting covaraince matrix elements and the like
  covMatrix getMagCovMatrix() const; //!< Returns a copy of the mag covariance matrix
  const covMatrix& getMagCovMatrixRef() const { return mag_covmatrix; } //!< Returns a reference to the mag covariance matrix
  float getMagCovElement(unsigned int,unsigned int) const; //!< Returns the corresponding element of the mag covariance matrix
  unsigned int getMagCovMatrixNumRows() const { return mag_covmatrix.getNRows(); } //!< Number of rows in the mag covariance matrix

  //width cov matrix
  covMatrix getWidthCovMatrix() const; //!< Return a copy of the width covariance matrix
  const covMatrix& getWidthCovMatrixRef() const { return width_covmatrix; } //!< Return a reference to the width covariance matrix
  unsigned int getWidthCovMatrixNumRows() const { return width_covmatrix.getNRows(); } //!< Number of rows in the width covariance matrix

  //colour cov matrix
  covMatrix getColourCovMatrix() const; //!< Return a copy of the colour covariance matrix
  const covMatrix& getColourCovMatrixRef() const { return colour_covmatrix; } //!< Return a reference to the colour covariance matrix
  unsigned int getColourCovMatrixNumRows() const { return colour_covmatrix.getNRows(); } //!< Number of rows in the colour covariance matrix

  //mag-width
  covMatrix getMagWidthCovMatrix() const; //!< Returns a copy of the mag-width covariance matrix
  const covMatrix& getMagWidthCovMatrixRef() const { return magwidth_covmatrix;
} //!< Returns a reference to the mag-width covariance matrix
  unsigned int getMagWidthCovMatrixNumRows() const { return magwidth_covmatrix.getNRows(); } //!< Number of rows in the mag-width covariance matrix

  //mag-colour
  covMatrix getMagColourCovMatrix() const; //!< Returns a copy of the mag-colour covariance matrix
  const covMatrix& getMagColourCovMatrixRef() const { return magcolour_covmatrix; } //!< Returns a reference to the mag-colour covariance matrix
  unsigned int getMagColourCovMatrixNumRows() const { return magcolour_covmatrix.getNRows(); } //!< Number of rows in the mag-colour covariance matrix

  //width-colour
  covMatrix getWidthColourCovMatrix() const; //!< Returns a copy of the width-colour covariance matrix
  const covMatrix& getWidthColourCovMatrixRef() const { return widthcolour_covmatrix; } //!< Returns a reference to the width-colour covariance matrix
  unsigned int getWidthColourCovMatrixNumRows() const { return widthcolour_covmatrix.getNRows(); } //!< Number of rows in the width-colour covariance matrix

 //Combined covariance matrix, given alpha, beta
  covMatrix getCombinedCovMatrix(float,float) const; //!< Return a copy of the combined \f$\alpha, \beta\f$ and constant cov matrix
  float getCombinedCovElement(float,float,unsigned int,unsigned int) const; //!< Get element of combined covariance matrix


  //Iterators
  iterator begin() { return points.begin(); } //!< Returns read/write iterator pointing to first element
  const_iterator begin() const { return points.begin(); } //!< Returns read only iterator pointing to first element
  iterator end() { return points.end(); } //!< Returns read/write iterator pointing to one past the last element
  const_iterator end() const { return points.end(); } //!< Returns read only iterator pointing to one past the last element

  //Reverse iterators
  reverse_iterator rbegin() { return points.rbegin(); } //!< Returns read/write reverse iterator pointing to last element
  const_reverse_iterator rbegin() const { return points.rbegin(); } //!< Returns read only reverse iterator pointing to last element
  reverse_iterator rend() { return points.rend(); } //!< Returns read/write reverse iterator pointing to one before first element
  const_reverse_iterator rend() const { return points.rend(); } //!< Returns read only reverse iterator pointing to one before first element

  //Size/capacity stuff
  unsigned int capacity() const { return points.capacity(); } //!< Space allocated (but not necessarily filled) in points
  unsigned int size() const { return points.size(); } //!< Number of elements in points
  bool empty() const { return points.empty(); } //!< Is points empty?
  void resize(unsigned int newsize) { 
    points.resize(newsize); 
    if (havemagcov) mag_covmatrix.resize( newsize, newsize );
    if (havewidthcov) width_covmatrix.resize( newsize, newsize );
    if (havecolourcov) colour_covmatrix.resize( newsize, newsize );
    if (havemagwidthcov) magwidth_covmatrix.resize( newsize, newsize );
    if (havemagcolourcov) magcolour_covmatrix.resize( newsize, newsize );
    if (havewidthcolourcov) widthcolour_covmatrix.resize( newsize, newsize );
  } //!< Resize points
  void reserve(unsigned int newsize) { points.reserve(newsize); } //!< Allocate but don't initialize more space in points
  
  //sorting
  void zcmbsort(); //!< Sort supernovae by zcmb

  //concatenation
  SNeData& operator+=(const SNeData&); //!< Concatenate list to end of current one

  //indexing
  SNeData operator[](const std::vector<unsigned int>&) const; //!< Return new list indexed from old
  SNeData operator[](const std::vector<std::string>&) const; //!< Return new list indexed by supernova name
  SNeDataEntry& operator[](const int i) { return points[i]; } //!< Unchecked subscript
  const SNeDataEntry& operator[](const int i) const { return points[i]; } //!< Unchecked subscript
  SNeDataEntry& at(const int i) { return points.at(i); } //!<Checked subscript
  const SNeDataEntry& at(const int i) const { return points.at(i); } //!< Checked subscript

  //Removal
  void remove(const std::vector<int>&); //!< Removes specified elements in plac
  void remove(const std::vector<std::string>&, bool strict=true); //!< Removes specified elements in place using names
  SNeData copy_remove(const std::vector<int>&) const; //!< Returns a new list with indexed entries removed
  SNeData copy_remove(const std::vector<std::string>&, bool strict=true) const; //!< Returns a new list with entries with names in argument removed
  
  std::set<unsigned int> getDataSetList() const; //!< Returns the datasets present

  bool isCovValid() const; //!< Check if covariance info is valid for all SN

};

/*!
  \brief Sorting class for SNeDataEntry based on zcmb.

  \ingroup snedata

  This differs from SNeSortByZcmb in that it actually sorts
  an auxillary array based on the contents of a SN list
  instead of actually sorting the entries themselves.  This is
  useful when dealing with the covariance matrix.
*/
class SNeSortByZcmbIndex {
private:
  SNeData& data;  //!< SNeData to sort by
public :
  SNeSortByZcmbIndex( SNeData& d ) : data( d ) { }  //!< Constructor
  /*!
    \brief Comparison operator
   */
  int operator()( unsigned int i, unsigned int j ) {
    return data[i].zcmb < data[j].zcmb;
  }
};

#endif
