#ifndef __chifunc__
#define __chifunc__

#include <vector>
#include <ostream>
#include <map>

#include <Minuit2/FCNBase.h>

#include <snedata.h>
#include <lumdist.h>
#include <utility.h>
#include <covmatrix.h>
#include <auxconstraint.h>

using namespace ROOT::Minuit2;

/*!
  \brief Abstract base class for evaluating the \f$\chi^2\f$ of SNe
*/
class chifunc_base : public FCNBase {
 protected:
  unsigned int nsn; //!< Number of supernova
  SNeData sne;      //!< Supernova data
  lumdist lmdist;   //!< Luminosity distance calculator

  mutable std::vector<double> dl; //!< Convenience vector for distances

  //Errors
  std::vector<double> pre_vars;  //!< Precalcluated error information

  utility::cosmo_fittype fittype;  //!< Type of fit being done

  std::map<unsigned int,double> intrinsicdisp; //!< Intrinsic dispersion in magnitudes
  double pecz;  //!< Peculiar velocity in redshift units
  double errdef;  //!< What errors to report -- see Minuit docs
  
  bool isprepped;  //!< Ready to evaluate \f$\chi^2\f$

  bool includebao; //!< Include the Baryon Acoustic Peak constraint
  bool includewmap; //!< Include WMAP 3rd year shift parameter constraint
  bool useEisenstein; //!< Use Eisenstein instead of Percival in BAO

  //Elements for auxilliary computations
  mutable auxconstraint::sdss_lrg_bao E05;
  mutable auxconstraint::bao_P09 P09;
  mutable auxconstraint::wmap7yr_dls WMAP7;
  mutable auxconstraint::wmap7_bao_P09 WMAP7_P09;

 public:

  chifunc_base( const SNeData& data, const std::map<unsigned int,double>& idisp, 
		double pz=0.001,utility::cosmo_fittype fit=utility::omegam_omegade); 

  double Up() const { return errdef; } //!< Return error defintion
  
  //This has to be run before the fit
  virtual void CalcPreErrs() = 0;  //!< Precalculate as many error quantities as possible

  virtual double GetRMS(const std::vector<double>&) const = 0; //!< Returns RMS around the fit

  //This is the important one
  virtual double operator()(const std::vector<double>&) const = 0; //!< Returns \f$\chi^2\f$. 
  virtual double EstimateScriptm(const std::vector<double>&) const = 0; //!< Estimates \f${\mathcal M}\f$

  void SetErrDef(double def) { errdef = def; } //!< Set error definition
  
  bool AreUsingKomatsuForm() const { return lmdist.usingKomatsuForm(); } //!< Are we using the Komatsu form for w(a)?
  void SetUseKomatsuForm(); //!< Use Komatsu form for w(a)
  void UnsetUseKomatsuForm(); //!< Use the Linder form for w(a)
  double GetAtrans() const { return lmdist.getAtrans(); } //!< Return \f$a_t\f$ for Komatsu w(a)
  void SetAtrans(double val); //!< Set \f$a_t\f$ for Komatsu form of w(a)

  void SetIncludeBAO(bool set) { includebao=set; } //!< Controls if BAO constraint included
  void SetIncludeWMAP(bool set) { includewmap=set; } //!< Controls if WMAP 7th year constraint is used
  void SetUseEisenstein(bool set) { useEisenstein=set; } //!< Controls if Eisensein '05 is used instead of Percival '09 for BAO

  virtual void print(std::ostream&) const; //!< Prints a summary of the data
  virtual void print(std::ostream& os, const std::vector<double>& par) const = 0; //!< More extedned printing function
};


/*!
  \brief Class for evaluating the \f$\chi^2\f$ of SNe including covariance 
     matricies
*/
class chifunc : public chifunc_base {
  
protected:
  mutable std::vector<double> working_vec; //!< Working variable


  mutable covMatrix invcovmatrix; //!< Holds current inverse cov matrix
  
  bool fixalphaset;  //!< Is \f$\alpha\f$ fixed for error propagation
  bool fixbetaset;   //!< Is \f$\beta\f$ fixed for error propagation
  double fixalphaval; //!< Value of fixed \f$\alpha\f$
  double fixbetaval; //!< Value of fixed \f$\beta\f$

  bool errsfixed;  //!< Are errors fixed?

  void computeInvCovMatrix(float,float) const; //!< Get inverse covariance matrix
  double GetChisq(double,double,double,double,double,double,double) const; //!< Non-diagonal case
  double GetDiagonalChisq(double,double,double,double,double,double,double) const; //!< Diagonal case

public:
  
  chifunc( const SNeData& data, const std::map<unsigned int,double>& idisp, 
	   double pz=0.001,utility::cosmo_fittype fit=utility::omegam_omegade); 
  void SetAlphaErrFix( double ); //!< Set \f$\alpha\f$ to fixed value for propagation
  void SetBetaErrFix( double ); //!< Set \f$\beta\f$ to fixed value for propagation

  void CalcPreErrs();  //!< Precalculate as many error quantities as possible

  double GetRMS(const std::vector<double>&) const; //!< Returns RMS around the fit

  void print(std::ostream& os, const std::vector<double>& par) const; //!< More extedned printing function

  double operator()(const std::vector<double>&) const; //!< Returns \f$\chi^2\f$. 
  double EstimateScriptm(const std::vector<double>&) const; //!< Estimates \f${\mathcal M}\f$

};

/*!
  \brief A version of chifunc that supports broken-linear alpha and beta
  relations, but not covariance matricies.
*/
class chifunc_multilin : public chifunc_base {
 protected:
  double widthcut; //!< Break point in alpha correction
  double colourcut; //!< Break point in colour correction
  double GetChisq(double,double,double,double,double,double, double, 
		  double,double) const; //!< Gets \f$\chi^2\f$
public:
  
  chifunc_multilin( const SNeData& data, const std::map<unsigned int,double>& idisp, 
		    double pz=0.001,utility::cosmo_fittype fit=utility::omegam_omegade); 

  double GetWidthCut() const { return widthcut; } //!< Returns widthcut
  double GetColourCut() const { return colourcut; } //!< Returns widthcut
  void SetWidthCut(double val) { widthcut=val; isprepped=false;} //!< Sets widthcut
  void SetColourCut(double val) { colourcut=val; isprepped=false;} //!< Sets colourcut

  void CalcPreErrs();  //!< Precalculate as many error quantities as possible

  double GetRMS(const std::vector<double>&) const; //!< Returns RMS around the fit

  double operator()(const std::vector<double>&) const; //!< Returns \f$\chi^2\f$. 
  double EstimateScriptm(const std::vector<double>&) const; //!< Estimates \f${\mathcal M}\f$
  
  void print(std::ostream& os, const std::vector<double>& par) const; //!< More extedned printing function

};


#endif
