#ifndef __auxconstraint__
#define __auxconstraint__

#include <Minuit2/FCNBase.h>
#include <gsl/gsl_integration.h>

#include <utility.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

using namespace ROOT::Minuit2;

/*!
  \defgroup other_const Other cosmological constraints
*/

/*!
  \brief Auxillary constraints

  \ingroup other_const

  Code for handling Baryon Acoustic Peak constraints from
  <A HREF="http://xxx.lanl.gov/abs/astro-ph/0501171">Eisenstein '05</A>.
  Specifically, these are based on the A parameter of that paper.
*/
namespace auxconstraint {
  
  //PDG (2004) value
  //const ogammah2 = 2.471e-5; //!< \f$\Omega_{\gamma} h^2\f$

  //Komatsu et al. value
  const double ogammah2 = 2.469e-5; //!< \f$\Omega_{\gamma} h^2\f$

  const double neff = 3.04; //!< Number of Neutrino species
  const double oradh2 = ogammah2*(1.0+0.2271*neff); //!< \f$\Omega_{r} h^2\f$

  double int_term_w( double z, void* params ); //!< 1/E term from Eisenstein, w=const
  double int_term_w0_wa( double z, void* params ); //!< 1/E term from Eisenstein, \f$w_0, w_a\f$ model
  double int_term_w0_waKom( double z, void* params ); //!< 1/E term from Eisenstein, \f$w_0, w_a\f$ model, Komatsu form

  double one_over_E_rad( double z, void * params );  //!< \f$1/E, w=-1\f$ with radiation
  double one_over_E_w_rad( double z, void * params ); //!< \f$1/E\f$ w=const, with radiation
  double one_over_E_w0wa_rad( double z, void * params ); //!< \f$1/E , w_0, w_a\f$, with radiation
  double one_over_E_w0waKom_rad( double z, void * params ); //!< \f$1/E , w_0, w_a\f$, Komatsu form with radiation


  double rsHovercint( double a, void * params); //<! \f$r_s H_0 / c\f$ integral
  double rsHovercint_w( double a, void * params); //<! \f$r_s H_0 / c\f$ integral with \f$w\f$
  double rsHovercint_w0wa( double a, void * params); //<! \f$r_s H_0 / c\f$ integral with \f$w_0, w_a\f$
  double rsHovercint_w0waKom( double a, void * params); //<! \f$r_s H_0 / c\f$ integral with \f$w_0, w_a\f$, Komatsu form

  bool isRebounding(double, double); //!< Check for rebounding Universe if w=-1

  //Can't be part of member class because we are interfacing to the GSL,
  // which doesn't handle member functions
  double wmap7_shift_chisq(const gsl_vector *v, void *params); //<! \f$\chi^2\f$ of WMAP7 shift parameters
  double percival_chisq(const gsl_vector *v, void *params); //<! \f$\chi^2\f$ of Percival BAO fits
  double percival_wmap7_chisq(const gsl_vector *v, void *params); //<! \f$\chi^2\f$ of Percival BAO + WMAP7 fits

  /*!
    \brief Helper class for some functions common to WMAP shift and BAO
  */
  class distance_helper {
  private:
    gsl_integration_workspace *work; //!< Workspace for numeric integratio
    bool useKomatsuForm; //!< Use Komatsu et al. (2008) form for \f$w\left(a\right)\f$
    double atrans; //!< Transition scale factor if Komatsu form used
  public:
    distance_helper(); //!< Constructor
    ~distance_helper(); //!< Destructor

    double GetAtrans() const { return atrans; } //!< Returns transition a
    void SetAtrans(double val) { atrans = val; } //!< Sets transition a
    bool UsingKomatsuForm() const { return useKomatsuForm; } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { useKomatsuForm=true; } //!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { useKomatsuForm = false; } //!< Turns off the Komatsu form

    double GetK1( double, double, double, double, double, double, double, bool& ) const; //<! Get \f$K_{1} = H_0 r_{s} \left( z  \right) / c\f$
    double GetK2( double, double, double, double, double, double, bool& ) const; //<! Get \f$K_{w} = \left(1 + z_{\star}\right) H_0 D_{A} \left( z \right) / c\f$
    double GetOneOverE(double, double, double, double, double, double, bool& ) const; //!< Get \f$1/E = H_0 / H\f$

    double GetZstar( double, double ) const; //<! \f$z_{\star}\f$
    double GetZdrag( double, double ) const; //<! \f$z_d\f$

  };

  /*!
    \brief Class for Baryon Acoustic Peak from the SDSS LRGs
  */
  class sdss_lrg_bao : public FCNBase {
  private:
    utility::cosmo_fittype fittype;  //!< Type of fit being done
    double ns; //!< Matter spectral index
    double fnu; //!< Neutrino mass fraction \f$\Omega_{\nu}/\Omega_m\f$
    double errdef; //!< What errors to report -- see Minuit docs
    bool useKomatsuForm; //!< Use the Komatsu for for w(a)?
    double atrans; //!< Transition scale factor \f$a_t\f$ for Komatsu form
  public:
    sdss_lrg_bao(utility::cosmo_fittype=utility::omegam_omegade,
		 double ns=0.98, double fnu=0.0); //!< Constructor
    void setFittype(utility::cosmo_fittype ft) { fittype=ft; } //!< Sets fit type
    double Up() const { return errdef; } //!< Return error defintion
    void SetNs( double val ) { ns = val; } //!< Set Matter spectral index
    void SetFnu( double val ) { fnu = val; } //!< Set Neutrino mass fraction
    void SetUseKomatsuForm() { useKomatsuForm=true; } //!< Use the Komatsu Form for w(a)
    void UnsetUseKomatsuForm() { useKomatsuForm=true; } //!< Don't use the Komatsu Form for w(a)
    void SetAtrans(double val) { atrans=val; } //!< Set \f$a_t\f$ for Komatsu form

    double aval(double om, double ode, double w=-1, 
		double wa=0.0) const; //!< Value of A parameter
    double aval(const std::vector<double>&) const; //!< Value of A parameter
    double operator()(double om, double ode, double w=-1,
		      double wa=0.0) const; //!< Returns \f$\chi^2\f$.
    double operator()(const std::vector<double>&) const; //!< Returns \f$\chi^2\f$.
    void SetErrDef(double def) { errdef = def; } //!< Set error definition
    
  };


  /*!
    \brief Class for Percival 09 BAO constraints
  */
  class bao_P09 : public FCNBase {
  private:
    utility::cosmo_fittype fittype;  //!< Type of fit being done
    double errdef; //!< What errors to report -- see Minuit docs
    gsl_vector * v; //!< Internal use vector for simplex minimization
    gsl_vector *ss; //!< Step sizes for minimization
    gsl_multimin_fminimizer *s;
    mutable gsl_multimin_function minex_func;
    static const gsl_multimin_fminimizer_type *T; //!< Type of fitter

    mutable distance_helper dhelp;
  public:
    bao_P09(utility::cosmo_fittype=utility::omegam_omegade); //!< Constructor
    ~bao_P09(); //!< Destructor
    void setFittype(utility::cosmo_fittype ft) { fittype=ft; } //!< Sets fit type
    double Up() const { return errdef; } //!< Return error defintion

    double GetAtrans() const { return dhelp.GetAtrans(); } //!< Returns transition a
    void SetAtrans(double val) { dhelp.SetAtrans(val); } //!< Sets transition a
    bool UsingKomatsuForm() const { return dhelp.UsingKomatsuForm(); } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { dhelp.SetUseKomatsuForm(); }//!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { dhelp.UnsetUseKomatsuForm(); } //!< Turns off the Komatsu form

    //Chi^2 code
    double operator()(double om, double ode, double w=-1,
		      double wa=0.0) const; //!< Returns \f$\chi^2\f$.
    double operator()(const std::vector<double>&) const; //!< Returns \f$\chi^2\f$.

    void SetErrDef(double def) { errdef = def; } //!< Set error definition
    
  };

  /*!
    \brief Class for WMAP 7th year distance to last scattering

    Uses the method described in Komatsu et al. (2009)
  */
  class wmap7yr_dls : public FCNBase {
  private:
    utility::cosmo_fittype fittype;  //!< Type of fit being done
    double errdef; //!< What errors to report -- see Minuit docs
    gsl_vector * v; //!< Internal use vector for simplex minimization
    gsl_vector *ss; //!< Step sizes for minimization
    gsl_multimin_fminimizer *s;
    mutable gsl_multimin_function minex_func;
    static const gsl_multimin_fminimizer_type *T; //!< Type of fitter

    mutable distance_helper dhelp;
  public:
    wmap7yr_dls(utility::cosmo_fittype=utility::omegam_omegade); //!< Constructor
    ~wmap7yr_dls(); //!< Destructor
    void setFittype(utility::cosmo_fittype ft) { fittype=ft; } //!< Sets fit type
    double Up() const { return errdef; } //!< Return error defintion

    double GetAtrans() const { return dhelp.GetAtrans(); } //!< Returns transition a
    void SetAtrans(double val) { dhelp.SetAtrans(val); } //!< Sets transition a
    bool UsingKomatsuForm() const { return dhelp.UsingKomatsuForm(); } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { dhelp.SetUseKomatsuForm(); }//!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { dhelp.UnsetUseKomatsuForm(); } //!< Turns off the Komatsu form

    //Chi^2 code
    double operator()(double om, double ode, double w=-1,
		      double wa=0.0) const; //!< Returns \f$\chi^2\f$.
    double operator()(const std::vector<double>&) const; //!< Returns \f$\chi^2\f$.

    void SetErrDef(double def) { errdef = def; } //!< Set error definition
    
  };


  /*!
    \brief Class for Percival 09 BAO + WMAP7 constraints
  */
  class wmap7_bao_P09 : public FCNBase {
  private:
    utility::cosmo_fittype fittype;  //!< Type of fit being done
    double errdef; //!< What errors to report -- see Minuit docs
    gsl_vector * v; //!< Internal use vector for simplex minimization
    gsl_vector *ss; //!< Step sizes for minimization
    gsl_multimin_fminimizer *s;
    mutable gsl_multimin_function minex_func;
    static const gsl_multimin_fminimizer_type *T; //!< Type of fitter

    mutable distance_helper dhelp;
  public:
    wmap7_bao_P09(utility::cosmo_fittype=utility::omegam_omegade); //!< Constructor
    ~wmap7_bao_P09(); //!< Destructor
    void setFittype(utility::cosmo_fittype ft) { fittype=ft; } //!< Sets fit type
    double Up() const { return errdef; } //!< Return error defintion

    double GetAtrans() const { return dhelp.GetAtrans(); } //!< Returns transition a
    void SetAtrans(double val) { dhelp.SetAtrans(val); } //!< Sets transition a
    bool UsingKomatsuForm() const { return dhelp.UsingKomatsuForm(); } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { dhelp.SetUseKomatsuForm(); }//!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { dhelp.UnsetUseKomatsuForm(); } //!< Turns off the Komatsu form

    //Chi^2 code
    double operator()(double om, double ode, double w=-1,
		      double wa=0.0) const; //!< Returns \f$\chi^2\f$.
    double operator()(const std::vector<double>&) const; //!< Returns \f$\chi^2\f$.

    void SetErrDef(double def) { errdef = def; } //!< Set error definition
    
  };


}


#endif
