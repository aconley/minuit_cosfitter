#include <string>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <stdexcept>

#include <utility.h>
#include <cosfitterexcept.h>
#include <snedata.h>

using namespace std;

/*!
  If it is necessary to calculate a large number of these, it is more
  efficient to do it by hand than to use this function.

  \param[in] alpha Width (stretch) correction coefficient
  \param[in] beta Colour correction coefficient
  \returns Corrected mangitude value
*/
float SNeDataEntry::getCorrectedMag(float alpha, float beta) const {
  return static_cast<float>(mag + alpha * (widthpar - 1.0) - beta * colourpar);
}

/*!
  If it is necessary to calculate a large number of these, it is more
  efficient to do it by hand than to use this function.

  \param[in] alpha Width (stretch) correction coefficient
  \param[in] beta Colour correction coefficient
  \returns Error in corrected magnitude, not including peculiar velocity
    or intrinsic dispersion
*/
float SNeDataEntry::getCorrectedMagError(float alpha, 
                                         float beta) const {
  float var = var_mag 
    + (alpha*alpha*var_widthpar) 
    + (beta*beta*var_colourpar) 
    + (2.0*alpha*cov_mag_widthpar) 
    - (2.0*beta*cov_mag_colourpar) 
    - (2.0*alpha*beta*cov_widthpar_colourpar);
  return sqrt(var);
}

/*!
  If it is necessary to calculate a large number of these, it is more
  efficient to do it by hand than to use this function

  \param[in] alpha Width (stretch) correction coefficient
  \param[in] beta Colour correction coefficient
  \returns Variance corrected magnitude, not including peculiar velocity
    or intrinsic dispersion
*/
float SNeDataEntry::getCorrectedMagVar(float alpha, 
                                       float beta) const {
  return var_mag + alpha*alpha*var_widthpar + 
    beta*beta*var_colourpar + 2*alpha*cov_mag_widthpar -
    2*beta*cov_mag_colourpar - 2*alpha*beta*cov_widthpar_colourpar;
}

/*!
  Makes sure that all the defined correlation coefficients are
  less than 1.
 */
bool SNeDataEntry::isCovValid() const {
  if ( (cov_mag_widthpar > 0) && (var_mag > 0) && (var_widthpar > 0) &&
       ( fabs(cov_mag_widthpar/sqrt(var_mag*var_widthpar)) > 1 ) ) return false;
  if ( (cov_mag_colourpar > 0) && (var_mag > 0) && (var_colourpar > 0) &&
       ( fabs(cov_mag_colourpar/sqrt(var_mag*var_colourpar)) > 1) ) 
    return false;
  if ( (cov_widthpar_colourpar > 0) && (var_widthpar > 0) && 
       (var_colourpar > 0) &&
       ( fabs(cov_widthpar_colourpar/sqrt(var_widthpar*var_colourpar)) > 1) )
			      return false;
  return true;
}

SNeData::SNeData() : Name(""), isdiagonal( true ), havemagcov( false ),
		     havewidthcov( false ), havecolourcov( false ),
		     havemagwidthcov( false ), havemagcolourcov( false ),
		     havewidthcolourcov( false ) { }


SNeData::SNeData(const std::string& name) : Name(name), isdiagonal( true ),
					    havemagcov( false ), 
					    havewidthcov( false ),
					    havecolourcov( false ),
					    havemagwidthcov( false ),
					    havemagcolourcov( false ),
					    havewidthcolourcov( false ) { }

SNeData::SNeData(const std::string& name, 
		 const std::string &fileName) : Name(name), 
						isdiagonal( true ),
						havemagcov( false ), 
						havewidthcov( false ),
						havecolourcov( false ),
						havemagwidthcov( false ),
						havemagcolourcov( false ),
						havewidthcolourcov( false ) {
  readData(fileName,false);
}

void SNeData::zcmbsort() {
  if (points.size() == 0) return; //Nothing to do
  if ( isdiagonal ) {
    sort ( points.begin(), points.end(), SNeSortByZcmb() );
  } else {
    //More complex case -- we have a covariance matrix
    // to drag around
    unsigned int nsn = points.size();

    SNeSortByZcmbIndex srt( *this );
    std::vector<unsigned int> indxarr( nsn );
    for (unsigned int i = 0; i < nsn; ++i) indxarr[i] = i;
    sort( indxarr.begin(), indxarr.end(), srt );

    //Sort data points
    std::vector<SNeDataEntry> newpnts( nsn );
    for (unsigned int i = 0; i < nsn; ++i)
      newpnts[i] = points[ indxarr[i] ];
    points = newpnts;

    //The covariance matricies
    if ( havemagcov || havewidthcov || havecolourcov || havemagwidthcov ||
	 havemagcolourcov || havewidthcolourcov ) {
      covMatrix newcov( nsn, 0.0 ); //Allocate once, use many times
      if (havemagcov) {
	for (unsigned int i = 0; i < nsn; ++i)
	  for (unsigned int j = 0; j < nsn; ++j)
	    newcov[i][j] = mag_covmatrix[ indxarr[i] ][ indxarr[j] ];
	mag_covmatrix = newcov;
      } 
      if (havewidthcov) {
	for (unsigned int i = 0; i < nsn; ++i)
	  for (unsigned int j = 0; j < nsn; ++j)
	    newcov[i][j] = width_covmatrix[ indxarr[i] ][ indxarr[j] ];
	width_covmatrix = newcov;
      } 
      if (havecolourcov) {
	for (unsigned int i = 0; i < nsn; ++i)
	  for (unsigned int j = 0; j < nsn; ++j)
	    newcov[i][j] = colour_covmatrix[ indxarr[i] ][ indxarr[j] ];
	colour_covmatrix = newcov;
      } 
      if (havemagwidthcov) {
	for (unsigned int i = 0; i < nsn; ++i)
	  for (unsigned int j = 0; j < nsn; ++j)
	    newcov[i][j] = magwidth_covmatrix[ indxarr[i] ][ indxarr[j] ];
	magwidth_covmatrix = newcov;
      } 
      if (havemagcolourcov) {
	for (unsigned int i = 0; i < nsn; ++i)
	  for (unsigned int j = 0; j < nsn; ++j)
	    newcov[i][j] = magcolour_covmatrix[ indxarr[i] ][ indxarr[j] ];
	magcolour_covmatrix = newcov;
      } 
      if (havewidthcolourcov) {
	for (unsigned int i = 0; i < nsn; ++i)
	  for (unsigned int j = 0; j < nsn; ++j)
	    newcov[i][j] = widthcolour_covmatrix[ indxarr[i] ][ indxarr[j] ];
	widthcolour_covmatrix = newcov;
      } 
    }
  }
}

/*!
Each line in the data file should have the format:
snname zcmb zhel dz widthpar dwidthpar colourpar dcolourpar cov_mag_width cov_mag_colourpar cov_widthpar_colourpar [dataset]

snname is a string (quotes are not necessary, and in fact probably not
a good idea), but the others are all
floating point numbers.  There can't be any spaces in snname, or bad
things will happen.  The lines must have an entry for each of the
above parameters, although they can be set to zero.

dataset is optional, and is an integer.  This is to support different
 dispersions for different datasets, which is controlled via the
SNDISP option.  It is up to you to ensure that, if you specify datasets,
you specify a dispersion for each.

Lines that begin with # are ignored.

Any old data is eliminated, as is the covariance matrix if present.

Some examples can be found in the test subdirectory.

*/
void SNeData::readData(const std::string& FileName, bool verbose) { 
  ifstream fl(FileName.c_str());

  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Data file " << FileName << " not found";
    throw CosFitterExcept("SNeData","readData",errstrng.str(),2);
  }

  std::string line;
  std::vector<std::string> words;
  std::stringstream str("");
  SNeDataEntry sndata;

  float dz, dmag, dwidthpar, dcolourpar;
  unsigned int wsize;

  points.resize(0);

  //File reading loop
  while(fl) {
    getline(fl,line);
    if (line[0]=='#' || line[0]==';') continue; //Comment line
    utility::stringwords(line,words);
    wsize = words.size();
    if (wsize < 13) continue; //Not enough entries on line
    
   if (wsize > 14) {
      std::stringstream errstrng("");
      errstrng << "Too many entries on line: " << line;
      throw CosFitterExcept("SNeData","readData",
			    errstrng.str(),4);
    }

    sndata.name = words[0];

    //Yes, I'm sure there is a more elegant way to do this.
    str.str(words[1]); str.clear(); str >> sndata.zcmb;
    str.str(words[2]); str.clear(); str >> sndata.zhel;
    str.str(words[3]); str.clear(); str >> dz;
    str.str(words[4]); str.clear(); str >> sndata.mag;
    str.str(words[5]); str.clear(); str >> dmag;
    str.str(words[6]); str.clear(); str >> sndata.widthpar;
    str.str(words[7]); str.clear(); str >> dwidthpar;
    str.str(words[8]); str.clear(); str >> sndata.colourpar;
    str.str(words[9]); str.clear(); str >> dcolourpar;
    str.str(words[10]); str.clear(); str >> sndata.cov_mag_widthpar;
    str.str(words[11]); str.clear(); str >> sndata.cov_mag_colourpar;
    str.str(words[12]); str.clear(); str >> sndata.cov_widthpar_colourpar;

    if (wsize == 14) {
      str.str(words[13]); str.clear(); str >> sndata.dataset;
    }

    sndata.var_z = dz*dz;
    sndata.var_mag = dmag*dmag;
    sndata.var_widthpar = dwidthpar*dwidthpar;
    sndata.var_colourpar = dcolourpar*dcolourpar;

    points.push_back(sndata);
  }
  
  if (points.size() == 0) throw CosFitterExcept("SNeData","readData",
						"No data read",8);

  if (verbose) cout << "Read: " << points.size() << " lines from " << 
		 FileName << endl;

  fl.close();

  //Clean out old covmatrix if we are rereading
  if ( !isdiagonal ) {
    isdiagonal = true;
    if (havemagcov) mag_covmatrix.resize(0);
    havemagcov = false;
    if (havewidthcov) width_covmatrix.resize(0);
    havewidthcov = false;
    if (havecolourcov) colour_covmatrix.resize(0);
    havecolourcov = false;
    if (havemagwidthcov) magwidth_covmatrix.resize(0);
    havemagwidthcov = false;
    if (havemagcolourcov) magcolour_covmatrix.resize(0);
    havemagcolourcov = false;
    if (havewidthcolourcov) widthcolour_covmatrix.resize(0);
    havewidthcolourcov = false;
  }

  //Check covariance validity
  if (!isCovValid())
    throw CosFitterExcept("SNeData","readData",
			  "Found invalid correlation coefficient",16);


}

/*!
  Internal function to read in covariance matricies and 
  make sure we don't have strange inconsistencies.
  \param[in] FileName  File to read from
  \param[out] covmat   Matrix to load it into
*/
void SNeData::readCovData(const std::string& FileName,
			  covMatrix& covmat) {
  ifstream fl(FileName.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Data file " << FileName << " not found";
    throw CosFitterExcept("SNeData","readCovData",errstrng.str(),2);
  }
  fl >> covmat;
  fl.close();
}
			      

/*!
  Reads in the magnitude covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void SNeData::readMagCovData(const std::string& FileName) { 
  readCovData( FileName, mag_covmatrix );
  unsigned int n;
  n = mag_covmatrix.getNRows();
  if ( n != 0 ) {
    isdiagonal = false;
    havemagcov = true;
  }
}

/*!
  Reads in the width covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void SNeData::readWidthCovData(const std::string& FileName) { 
  readCovData( FileName, width_covmatrix );
  unsigned int n;
  n = width_covmatrix.getNRows();
  if ( n != 0 ) {
    isdiagonal = false;
    havewidthcov = true;
  }
}

/*!
  Reads in the colour covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void SNeData::readColourCovData(const std::string& FileName) { 
  readCovData( FileName, colour_covmatrix );
  unsigned int n;
  n = colour_covmatrix.getNRows();
  if ( n != 0 ) {
    isdiagonal = false;
    havecolourcov = true;
  }
}

/*!
  Reads in the mag-width covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void SNeData::readMagWidthCovData(const std::string& FileName) { 
  readCovData( FileName, magwidth_covmatrix );
  unsigned int n;
  n = magwidth_covmatrix.getNRows();
  if ( n != 0 ) {
    isdiagonal = false;
    havemagwidthcov = true;
  }
}


/*!
  Reads in the mag-colour covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void SNeData::readMagColourCovData(const std::string& FileName) { 
  readCovData( FileName, magcolour_covmatrix );
  unsigned int n;
  n = magcolour_covmatrix.getNRows();
  if ( n != 0 ) {
    isdiagonal = false;
    havemagcolourcov = true;
  }
}

/*!
  Reads in the width-colour covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void SNeData::readWidthColourCovData(const std::string& FileName) { 
  readCovData( FileName, widthcolour_covmatrix );
  unsigned int n;
  n = widthcolour_covmatrix.getNRows();
  if ( n != 0 ) {
    isdiagonal = false;
    havewidthcolourcov = true;
  }
}

/*!
  Uses the same format as SNeData::readData.
  \param[in] outfile File to write to
*/
void SNeData::writeData(const std::string& outfile) const {  
  //Since we are doing formatted output, and the c++ style formatted
  // output sucks out loud, use C style

  FILE *fp;
  fp = fopen(outfile.c_str(),"w");
  if (!fp) {
    std::stringstream errstrng("");
    errstrng <<  "Couldn't open file " << outfile << " for writing";
    throw CosFitterExcept("SNeData","writeFile",errstrng.str(),1);
  }

  std::string hdfmt="%-10s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s ";
  hdfmt += "%-11s %-11s %-11s";
  std::string fmt = "%-10s %8.6f %8.6f %8.6f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f ";
  fmt += "%11.8f %11.8f %11.8f";

  //Write the header
  fprintf(fp,hdfmt.c_str(),"name","zcmb","zhel","dz","mag","dmag",
	  "width","dwidth","colour","dcolour","cov_m_w","cov_m_c",
	  "cov_w_c");

  //And the data
  for (std::vector<SNeDataEntry>::const_iterator sn = points.begin(); 
       sn != points.end(); ++sn)
    fprintf(fp,fmt.c_str(),sn->name.c_str(),sn->zcmb,sn->zhel,sqrt(sn->var_z),
	    sn->mag,sqrt(sn->var_mag),sn->widthpar,sqrt(sn->var_widthpar),
	    sn->colourpar,sqrt(sn->var_colourpar),sn->cov_mag_widthpar,
	    sn->cov_mag_colourpar,sn->cov_widthpar_colourpar);

  fclose(fp);
}

/*!
  Writes a full covaraince matrix.  Right now this just uses the overloaded
  operator<< for covMatrix
  \param[in] FileName File to write to
  \param[in] covmat  Covariance matrix to write
*/
void SNeData::writeCovData(const std::string& FileName, 
			   const covMatrix& covmat) const {
  ofstream fl(FileName.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Can't open " << FileName << " to write";
    throw CosFitterExcept("SNeData","writeCovFile",errstrng.str(),2);
  }
  fl << covmat;
  fl.close();
}

/*!
  Writes the mag covariance matrix.  
*/
void SNeData::writeMagCovData(const std::string& FileName) const { 
  if (!havemagcov) return;
  writeCovData(FileName,mag_covmatrix);
}

/*!
  Writes the width covariance matrix.  
*/
void SNeData::writeWidthCovData(const std::string& FileName) const { 
  if (!havewidthcov) return;
  writeCovData(FileName,width_covmatrix);
}

/*!
  Writes the colour covariance matrix.  
*/
void SNeData::writeColourCovData(const std::string& FileName) const { 
  if (!havecolourcov) return;
  writeCovData(FileName,colour_covmatrix);
}

/*!
  Writes the mag-width covariance matrix.  
*/
void SNeData::writeMagWidthCovData(const std::string& FileName) const { 
  if (!havemagwidthcov) return;
  writeCovData(FileName,magwidth_covmatrix);
}

/*!
  Writes the mag-colour covariance matrix.  
*/
void SNeData::writeMagColourCovData(const std::string& FileName) const { 
  if (!havemagcolourcov) return;
  writeCovData(FileName,magcolour_covmatrix);
}

/*!
  Writes the width-colour covariance matrix.  
*/
void SNeData::writeWidthColourCovData(const std::string& FileName) const { 
  if (!havewidthcolourcov) return;
  writeCovData(FileName,widthcolour_covmatrix);
}

/*!
  Use getMagCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of mag_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix SNeData::getMagCovMatrix() const {
  if (!havemagcov) return covMatrix();
  return mag_covmatrix;
}

/*!
  Use getWidthCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of width_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix SNeData::getWidthCovMatrix() const {
  if (!havewidthcov) return covMatrix();
  return width_covmatrix;
}

/*!
  Use getColourCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of colour_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix SNeData::getColourCovMatrix() const {
  if (!havecolourcov) return covMatrix();
  return colour_covmatrix;
}

/*!
  Use getMagWidthCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of magwidth_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix SNeData::getMagWidthCovMatrix() const {
  if (!havemagwidthcov) return covMatrix();
  return magwidth_covmatrix;
}

/*!
  Use getMagColourCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of magcolour_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix SNeData::getMagColourCovMatrix() const {
  if (!havemagcolourcov) return covMatrix();
  return magcolour_covmatrix;
}

/*!
  Use getWidthColourCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of widthcolour_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix SNeData::getWidthColourCovMatrix() const {
  if (!havewidthcolourcov) return covMatrix();
  return widthcolour_covmatrix;
}


/*!
  An exception is thrown if the element is out of range (which
  includes the case where there is no mag covariance matrix).

  \return The specified element of the mag covariance matrix.  
*/
float SNeData::getMagCovElement(unsigned int i, 
				unsigned int j) const {
  if (!havemagcov) 
    throw CosFitterExcept("SNeData","getMagCovElement",
			  "No covariance matrix",1);

  unsigned int nrow = mag_covmatrix.getNRows();
  if ( i >= nrow || j >= nrow ) {
    std::stringstream errstr("");
    errstr << "Invalid index: " << i << " " << j << " size is: " <<
      nrow << " by " << nrow;
    throw CosFitterExcept("SNeData","getMagCovElement",errstr.str(),1);
  }

  return mag_covmatrix[i][j];
}

/*!
  \param[in] alpha  \f$\alpha\f$
  \param[in] beta   \f$\beta\f$
  \return The combined covariance matrix

*/
covMatrix SNeData::getCombinedCovMatrix(float alpha,float beta) const {
  covMatrix retval;
  int status = 0;
  if (havemagcov)  retval = mag_covmatrix; else 
    retval = covMatrix(points.size(),0.0);
  if (havewidthcov && alpha != 0.0) 
    retval.scalarMultAndAdd( alpha * alpha, width_covmatrix, status );
  if (havecolourcov && beta != 0.0 ) 
    retval.scalarMultAndAdd( beta * beta, colour_covmatrix, status );
  if (havemagwidthcov && alpha != 0.0) 
    retval.scalarMultAndAdd( 2.0 * alpha, magwidth_covmatrix, status );
  if (havemagcolourcov && beta != 0.0) 
    retval.scalarMultAndAdd( - 2.0 * beta, magcolour_covmatrix, status );
  if (havewidthcolourcov && alpha != 0.0 && beta != 0.0) 
    retval.scalarMultAndAdd( - 2.0 * alpha * beta, widthcolour_covmatrix, 
			     status );
  if (status) throw CosFitterExcept("snedata","getCombinedCovMatrix",
				    "Error combining cov matricies",1);

  return retval;
}

/*!
  \param[in] alpha \f$\alpha\f$
  \param[in] beta  \f$\beta\f$
  \param[in] i     Element index
  \param[in] j     Element index
  \return The corresponding element of the combined covariance matrix

 */
float SNeData::getCombinedCovElement(float alpha, float beta,
				     unsigned int i,unsigned int j) const 
{
  if (isdiagonal) return 0;
  float val;
  val = 0.0;
  if (havemagcov) val += mag_covmatrix[i][j];
  if (havewidthcov) val += alpha * width_covmatrix[i][j];
  if (havecolourcov)  val += beta * colour_covmatrix[i][j];
  if (havemagwidthcov) val += 2.0 * alpha * magwidth_covmatrix[i][j];
  if (havemagcolourcov) val -= 2.0 * beta * magcolour_covmatrix[i][j];
  if (havewidthcolourcov) val -= 2.0 * alpha * beta * 
    widthcolour_covmatrix[i][j];
  return val;
}


/*!
  Note that if a covariance matrix is present the cross terms
  between the old SN and the new ones will be left blank.

  \param[in] a List to append to current list
  \returns Reference to current list
*/
SNeData& SNeData::operator+=(const SNeData& a) {

  if (&a == this) {
    //Self addition.  It's an interesting question if
    // we should respect the const reference to a or not --
    // here we just ignore the request
    return *this;
  }
   
  unsigned int oldsize1 = points.size();
  unsigned int oldsize2 = a.points.size();
  unsigned int newsize = oldsize1 + oldsize2;
  points.reserve( newsize );
  for (SNeData::const_iterator i = a.begin(); i != a.end(); ++i)
    points.push_back( *i );
  
  if ( havemagcov || havewidthcov || havecolourcov || havemagwidthcov ||
       havemagcolourcov || havewidthcolourcov || a.havemagcov ||
       a.havewidthcov || a.havecolourcov || a.havemagwidthcov || 
       a.havemagcolourcov || a.havewidthcolourcov ) {
    covMatrix newcovmat( newsize );
    if ( havemagcov ) {
      if ( a.havemagcov ) {
	//Both have covariance matricies
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = mag_covmatrix[i][j];
	for (unsigned int i = 0; i < oldsize2; ++i) 
	  for (unsigned int j = 0; j < oldsize2; ++j) 
	    newcovmat[oldsize1+i][oldsize1+j] = a.mag_covmatrix[i][j];
      } else {
	//Only the old matrix is diagonal
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = mag_covmatrix[i][j];
      }
      isdiagonal = false;
      havemagcov = true;
      mag_covmatrix = newcovmat;
    } else {
      if ( a.havemagcov ) {
	//The new one has a covariance matrix
	//Only the old matrix is diagonal
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i+oldsize1][j+oldsize1] = a.mag_covmatrix[i][j];
	isdiagonal = false;
	havemagcov = true;
	mag_covmatrix = newcovmat;
      }
    }
    
    if ( havewidthcov ) {
      newcovmat = 0.0; //Wipe out cov vals to be safe
      if ( a.havewidthcov ) {
	//Both have width covariance matricies
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = width_covmatrix[i][j];
	for (unsigned int i = 0; i < oldsize2; ++i) 
	  for (unsigned int j = 0; j < oldsize2; ++j) 
	    newcovmat[oldsize1+i][oldsize1+j] = a.width_covmatrix[i][j];
      } else {
	//Only the old matrix is diagonal
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = width_covmatrix[i][j];
      }
      isdiagonal = false;
      havewidthcov = true;
      width_covmatrix = newcovmat;
    } else {
      if ( a.havewidthcov ) {
	//The new one has a covariance matrix
	//Only the old matrix is diagonal
	newcovmat = 0.0;
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i+oldsize1][j+oldsize1] = a.width_covmatrix[i][j];
	isdiagonal = false;
	havewidthcov = true;
	width_covmatrix = newcovmat;
      }
    }
    if ( havecolourcov ) {
      newcovmat = 0.0;
      if ( a.havecolourcov ) {
	//Both have colour covariance matricies
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = colour_covmatrix[i][j];
	for (unsigned int i = 0; i < oldsize2; ++i) 
	  for (unsigned int j = 0; j < oldsize2; ++j) 
	    newcovmat[oldsize1+i][oldsize1+j] = a.colour_covmatrix[i][j];
      } else {
	//Only the old matrix is diagonal
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = colour_covmatrix[i][j];
      }
      isdiagonal = false;
      havecolourcov = true;
      colour_covmatrix = newcovmat;
    } else {
      if ( a.havecolourcov ) {
	//The new one has a colour covariance matrix
	//Only the old matrix is diagonal
	newcovmat = 0.0;
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i+oldsize1][j+oldsize1] = a.colour_covmatrix[i][j];
	isdiagonal = false;
	havecolourcov = true;
	colour_covmatrix = newcovmat;
      }
    }
    if ( havemagwidthcov ) {
      newcovmat = 0.0;
      if ( a.havemagwidthcov ) {
	//Both have mag-width covariance matricies
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = magwidth_covmatrix[i][j];
	for (unsigned int i = 0; i < oldsize2; ++i) 
	  for (unsigned int j = 0; j < oldsize2; ++j) 
	    newcovmat[oldsize1+i][oldsize1+j] = a.magwidth_covmatrix[i][j];
      } else {
	//Only the old matrix is diagonal
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = magwidth_covmatrix[i][j];
      }
      isdiagonal = false;
      havemagwidthcov = true;
      magwidth_covmatrix = newcovmat;
    } else {
      if ( a.havemagwidthcov ) {
	//The new one has a covariance matrix
	//Only the old matrix is diagonal
	newcovmat = 0.0;
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i+oldsize1][j+oldsize1] = a.magwidth_covmatrix[i][j];
	isdiagonal = false;
	havemagwidthcov = true;
	magwidth_covmatrix = newcovmat;
      }
    }
    if ( havemagcolourcov ) {
      newcovmat = 0.0;
      if ( a.havemagcolourcov ) {
	//Both have mag-colour covariance matricies
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = magcolour_covmatrix[i][j];
	for (unsigned int i = 0; i < oldsize2; ++i) 
	  for (unsigned int j = 0; j < oldsize2; ++j) 
	    newcovmat[oldsize1+i][oldsize1+j] = a.magcolour_covmatrix[i][j];
      } else {
	//Only the old matrix is diagonal
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = magcolour_covmatrix[i][j];
      }
      isdiagonal = false;
      havemagcolourcov = true;
      magcolour_covmatrix = newcovmat;
    } else {
      if ( a.havemagcolourcov ) {
	//The new one has a covariance matrix
	//Only the old matrix is diagonal
	newcovmat = 0.0;
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i+oldsize1][j+oldsize1] = a.magcolour_covmatrix[i][j];
	isdiagonal = false;
	havemagcolourcov = true;
	magcolour_covmatrix = newcovmat;
      }
    }
    if ( havewidthcolourcov ) {
      newcovmat = 0.0;
      if ( a.havewidthcolourcov ) {
	//Both have width-colour covariance matricies
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = widthcolour_covmatrix[i][j];
	for (unsigned int i = 0; i < oldsize2; ++i) 
	  for (unsigned int j = 0; j < oldsize2; ++j) 
	    newcovmat[oldsize1+i][oldsize1+j] = a.widthcolour_covmatrix[i][j];
      } else {
	//Only the old matrix is diagonal
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i][j] = widthcolour_covmatrix[i][j];
      }
      isdiagonal = false;
      havewidthcolourcov = true;
      widthcolour_covmatrix = newcovmat;
    } else {
      if ( a.havewidthcolourcov ) {
	//The new one has a covariance matrix
	//Only the old matrix is diagonal
	newcovmat = 0.0;
	for (unsigned int i = 0; i < oldsize1; ++i) 
	  for (unsigned int j = 0; j < oldsize1; ++j) 
	    newcovmat[i+oldsize1][j+oldsize1] = a.widthcolour_covmatrix[i][j];
	isdiagonal = false;
	havewidthcolourcov = true;
	widthcolour_covmatrix = newcovmat;
      }
    }

  }

  return *this;
}


/*!
  \param[in] indxarr Array of indicies into old array.  
  \returns A new SNeData with only the specified objects present

  No bounds checking is performed, so it is up to the caller
  to ensure that the indicies in indxarr are valid.
*/
SNeData SNeData::operator[](const std::vector<unsigned int>& indxarr) const {
  unsigned int nelements = indxarr.size();

  if (nelements == 0) return SNeData(Name);

  SNeData retval(Name);

  retval.points.resize( nelements );
  
  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    retval.points[i++] = points[indxarr[j]];

  if ( havemagcov ) {
    retval.mag_covmatrix.resize( nelements );
    retval.havemagcov = true; retval.isdiagonal = false;
    for (unsigned int i = 0, j = 0; j < nelements; ++j) {
      for (unsigned int k = 0, h = 0; h < nelements; ++h)
	retval.mag_covmatrix[i][k++] = mag_covmatrix[indxarr[j]][indxarr[h]];
      ++i;
    }
  }
  if ( havewidthcov ) {
    retval.width_covmatrix.resize( nelements );
    retval.havewidthcov = true; retval.isdiagonal = false;
    for (unsigned int i = 0, j = 0; j < nelements; ++j) {
      for (unsigned int k = 0, h = 0; h < nelements; ++h)
	retval.width_covmatrix[i][k++] = 
	  width_covmatrix[indxarr[j]][indxarr[h]];
      ++i;
    }
  }
  if ( havecolourcov ) {
    retval.colour_covmatrix.resize( nelements );
    retval.havecolourcov = true; retval.isdiagonal = false;
    for (unsigned int i = 0, j = 0; j < nelements; ++j) {
      for (unsigned int k = 0, h = 0; h < nelements; ++h)
	retval.colour_covmatrix[i][k++] = 
	  colour_covmatrix[indxarr[j]][indxarr[h]];
      ++i;
    }
  }
  if ( havemagwidthcov ) {
    retval.magwidth_covmatrix.resize( nelements );
    retval.havemagwidthcov = true; retval.isdiagonal = false;
    for (unsigned int i = 0, j = 0; j < nelements; ++j) {
      for (unsigned int k = 0, h = 0; h < nelements; ++h)
	retval.magwidth_covmatrix[i][k++] = 
	  magwidth_covmatrix[indxarr[j]][indxarr[h]];
      ++i;
    }
  }
  if ( havemagcolourcov ) {
    retval.magcolour_covmatrix.resize( nelements );
    retval.havemagcolourcov = true; retval.isdiagonal = false;
    for (unsigned int i = 0, j = 0; j < nelements; ++j) {
      for (unsigned int k = 0, h = 0; h < nelements; ++h)
	retval.magcolour_covmatrix[i][k++] = 
	  magcolour_covmatrix[indxarr[j]][indxarr[h]];
      ++i;
    }
  }
  if ( havewidthcolourcov ) {
    retval.widthcolour_covmatrix.resize( nelements );
    retval.havewidthcolourcov = true; retval.isdiagonal = false;
    for (unsigned int i = 0, j = 0; j < nelements; ++j) {
      for (unsigned int k = 0, h = 0; h < nelements; ++h)
	retval.widthcolour_covmatrix[i][k++] = 
	  widthcolour_covmatrix[indxarr[j]][indxarr[h]];
      ++i;
    }
  }

  return retval;
}

/*!
  \param[in] namearr Vector of strings corresponding to names
    of SN
  \returns A new SNeData with only the specified objects present

  Error checking is done, and an out_of_range exception is thrown
  if any of the passed in names are not found.
*/
SNeData SNeData::operator[](const std::vector<std::string>& namearr) const {

  unsigned int nelements = namearr.size();

  if (nelements == 0) return SNeData(Name);  

  //The indexing is a lot more complex in this case, especially
  // because we want to retain the order requested
  std::map<std::string,unsigned int> name_index_map;
  unsigned int ncurrent = points.size();
  for (unsigned int i = 0; i < ncurrent; ++i)
    name_index_map[ points[i].name ] = i;

  std::map<std::string,unsigned int>::iterator endpoint = name_index_map.end();
  std::map<std::string,unsigned int>::iterator location;
  
  std::vector< unsigned int > indxarr;
  for (unsigned int i = 0, j = 0; i < nelements; ++i) {
    location = name_index_map.find( namearr[i] );
    if (location == endpoint) {
      //Not found
      std::stringstream s("");
      s << "Name " << namearr[i] << " not found in " << Name;
      throw out_of_range(s.str());
    }
    indxarr[j++] = location->second;
  }
  return (*this)[ indxarr ];
}

/*!
  \param[in] indxarr List of indicies to remove 

  Duplicate indicies are allowed, but no range checking is performed,
  so it is up to the caller to ensure that only valid indicies are passed.
*/
void SNeData::remove( const std::vector<int>& indxarr ) {

  unsigned int nelements = points.size();
  unsigned int nremove = indxarr.size();

  if (nremove == 0) return;
  
  std::vector<SNeDataEntry> newpoints(nelements-nremove);

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  std::vector<bool> included(nelements, true);

  for (unsigned int i =0; i < nremove; ++i)
    included[ indxarr[i] ] = false;
  
  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) newpoints[i++] = points[j];
  
  points = newpoints;

  if ( havemagcov || havewidthcov || havecolourcov || havemagwidthcov ||
       havemagcolourcov || havewidthcolourcov ) {
    //Now things get messy -- we also have to clean out the covariance
    // matrix
    covMatrix newcov(nelements-nremove,0.0);
    if ( havemagcov ) {
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = mag_covmatrix[i][j];
	  ++icntr;
	}
      mag_covmatrix = newcov;
    }
    if ( havewidthcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = width_covmatrix[i][j];
	  ++icntr;
	}
      width_covmatrix = newcov;
    }
    if ( havecolourcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = colour_covmatrix[i][j];
	  ++icntr;
	}
      colour_covmatrix = newcov;
    }
    if ( havemagwidthcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = magwidth_covmatrix[i][j];
	  ++icntr;
	}
      magwidth_covmatrix = newcov;
    }
    if ( havemagcolourcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = magcolour_covmatrix[i][j];
	  ++icntr;
	}
      magcolour_covmatrix = newcov;
    }
    if ( havewidthcolourcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = widthcolour_covmatrix[i][j];
	  ++icntr;
	}
      widthcolour_covmatrix = newcov;
    }
  }

}

/*!
  \param[in] namearr List of names to remove.
  \param[in] strict Complain if one of the passed in names is
     not found in the list.

  Duplicate indicies are allowed.  Invalid entries in namearr
  result in a out_of_range error being thrown unless strict
  is set to false
*/
void SNeData::remove( const std::vector<std::string>& namearr,
		      bool strict ) {

  unsigned int nelements = points.size();
  unsigned int nremove = namearr.size();

  if (nremove == 0) return; //Nothing to remove

  //The indexing is a lot more complex in this case, especially
  // because we want to retain the order requested, and because
  // we allow duplicates

  std::map<std::string,unsigned int> name_index_map;
  for (unsigned int i = 0; i < nelements; ++i)
    name_index_map[ points[i].name ] = i;

  std::map<std::string,unsigned int>::iterator endpoint = name_index_map.end();
  std::map<std::string,unsigned int>::iterator location;

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  std::vector<bool> included(nelements, true);

  for (std::vector<std::string>::const_iterator i = namearr.begin();
       i != namearr.end(); ++i) {
    location = name_index_map.find( *i );
    if (location == endpoint) {
      //Not found
      if (strict) {
	std::stringstream s("");
	s << "Name " << *i << " not found in " << Name;
	throw out_of_range(s.str());
      }
      //Otherwise just ignore
    } else included[ location->second ] = false;
  }

  unsigned int nremaining = count( included.begin(), included.end(), true );

  std::vector<SNeDataEntry> newpoints(nremaining);

  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) newpoints[i++] = points[j];

  points = newpoints;

  if ( havemagcov || havewidthcov || havecolourcov || havemagwidthcov ||
       havemagcolourcov || havewidthcolourcov ) {
    //Now things get messy -- we also have to clean out the covariance
    // matricies
    covMatrix newcov(nelements-nremove,0.0); //Allocate once, use many
    if ( havemagcov ) {
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = mag_covmatrix[i][j];
	  ++icntr;
	}
      mag_covmatrix = newcov;
    }
    if ( havewidthcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = width_covmatrix[i][j];
	  ++icntr;
	}
      width_covmatrix = newcov;
    }
    if ( havecolourcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = colour_covmatrix[i][j];
	  ++icntr;
	}
      colour_covmatrix = newcov;
    }
    if ( havemagwidthcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = magwidth_covmatrix[i][j];
	  ++icntr;
	}
      magwidth_covmatrix = newcov;
    }
    if ( havemagcolourcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = magcolour_covmatrix[i][j];
	  ++icntr;
	}
      magcolour_covmatrix = newcov;
    }
    if ( havewidthcolourcov ) {
      newcov = 0.0;
      unsigned int icntr = 0, jcntr;
      for (unsigned int i = 0; i < nelements; ++i)
	if (included[i]) {
	  jcntr=0;
	  for (unsigned int j = 0; j < nelements; ++j)
	    if (included[j])
	      newcov[icntr][jcntr++] = widthcolour_covmatrix[i][j];
	  ++icntr;
	}
      widthcolour_covmatrix = newcov;
    }
  }

}


/*!
  \param[in] indxarr List of indicies to remove from returned copy

  Duplicate indicies are allowed, but no range checking is performed,
  so it is up to the caller to ensure that only valid indicies are passed.

*/
SNeData SNeData::copy_remove( const std::vector<int>& indxarr ) const {

  SNeData retval(Name);

  unsigned int nelements = points.size();
  unsigned int nremove = indxarr.size();

  if (nremove == 0) return retval = *this; //Nothing removed

  retval.points.resize(nelements-nremove);

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  std::vector<bool> included(nelements, true);

  for (unsigned int i =0; i < nremove; ++i)
    included[ indxarr[i] ] = false;
  
  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) retval.points[i++] = points[j];

  if ( havemagcov ) {
    //Now things get messy -- we also have to clean out the covariance
    // matrix
    retval.mag_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.mag_covmatrix[icntr][jcntr++] = mag_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havemagcov = true; retval.isdiagonal = false;
  }
  if ( havewidthcov ) {
    retval.width_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.width_covmatrix[icntr][jcntr++] = width_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havewidthcov = true; retval.isdiagonal = false;
  }
  if ( havecolourcov ) {
    retval.colour_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.colour_covmatrix[icntr][jcntr++] = colour_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havecolourcov = true; retval.isdiagonal = false;
  }
  if ( havemagwidthcov ) {
    retval.magwidth_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.magwidth_covmatrix[icntr][jcntr++] = 
	      magwidth_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havemagwidthcov = true; retval.isdiagonal = false;
  }
  if ( havemagcolourcov ) {
    retval.magcolour_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.magcolour_covmatrix[icntr][jcntr++] = 
	      magcolour_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havemagcolourcov = true; retval.isdiagonal = false;
  }
  if ( havewidthcolourcov ) {
    retval.widthcolour_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.widthcolour_covmatrix[icntr][jcntr++] = 
	      widthcolour_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havewidthcolourcov = true; retval.isdiagonal = false;
  }

  return retval;
}

/*!
  \param[in] namearr List of names to remove from returned copy
  \param[in] strict Complain if one of the passed in names is
     not found in the list.
  \returns Copy of SNeData with specified entries removed

  Duplicate indicies are allowed.  Invalid entries in namearr
  result in a out_of_range error being thrown unless strict
  is set to false.
*/
SNeData SNeData::copy_remove( const std::vector<std::string>& namearr,
			      bool strict) const {

  SNeData retval(Name);

  unsigned int nelements = points.size();
  unsigned int nremove = namearr.size();

  if (nremove == 0) return retval = *this; //Nothing removed

  //The indexing is a lot more complex in this case, especially
  // because we want to retain the order requested, and because
  // we allow duplicates

  std::map<std::string,unsigned int> name_index_map;
  for (unsigned int i = 0; i < nelements; ++i)
    name_index_map[ points[i].name ] = i;

  std::map<std::string,unsigned int>::iterator endpoint = name_index_map.end();
  std::map<std::string,unsigned int>::iterator location;

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  std::vector<bool> included(nelements, true);

  for (std::vector<std::string>::const_iterator i = namearr.begin();
       i != namearr.end(); ++i) {
    location = name_index_map.find( *i );
    if (location == endpoint) {
      //Not found
      if (strict) {
	std::stringstream s("");
	s << "Name " << *i << " not found in " << Name;
	throw out_of_range(s.str());
      }
      //Otherwise just ignore this one
    } else included[ location->second ] = false;
  }

  unsigned int nremaining = count( included.begin(), included.end(), true );

  retval.points.resize(nremaining);

  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) retval.points[i++] = points[j];

  if ( havemagcov ) {
    retval.mag_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.mag_covmatrix[icntr][jcntr++] = mag_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havemagcov = true; retval.isdiagonal = false;
  }
  if ( havewidthcov ) {
    retval.width_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.width_covmatrix[icntr][jcntr++] = width_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havewidthcov = true; retval.isdiagonal = false;
  }
  if ( havecolourcov ) {
    retval.colour_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.colour_covmatrix[icntr][jcntr++] = colour_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havecolourcov = true; retval.isdiagonal = false;
  }
  if ( havemagwidthcov ) {
    retval.magwidth_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.magwidth_covmatrix[icntr][jcntr++] = 
	      magwidth_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havewidthcov = true; retval.isdiagonal = false;
  }
  if ( havemagcolourcov ) {
    retval.magcolour_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.magcolour_covmatrix[icntr][jcntr++] = 
	      magcolour_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havemagcolourcov = true; retval.isdiagonal = false;
  }
  if ( havewidthcolourcov ) {
    retval.widthcolour_covmatrix.resize(nelements-nremove );
    unsigned int icntr = 0, jcntr;
    for (unsigned int i = 0; i < nelements; ++i) {
      if (included[i]) {
	jcntr=0;
	for (unsigned int j = 0; j < nelements; ++j)
	  if (included[j])
	    retval.widthcolour_covmatrix[icntr][jcntr++] = 
	      widthcolour_covmatrix[i][j];
	++icntr;
      }
    }
    retval.havewidthcolourcov = true; retval.isdiagonal = false;
  }

  return retval;
}

std::set<unsigned int> SNeData::getDataSetList() const {
  std::set<unsigned int> retset;
  for (std::vector<SNeDataEntry>::const_iterator pos = points.begin();
       pos != points.end(); ++pos) 
    retset.insert( pos->dataset );
  return retset;
}

bool SNeData::isCovValid() const {
  if (points.size() == 0) return true;
  bool is_valid = true;
  for (std::vector<SNeDataEntry>::const_iterator it = points.begin();
       it != points.end(); ++it)
    is_valid &= it->isCovValid();
  return is_valid;
}
