#include<cmath>
#include<string>
#include<sstream>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

#include <auxconstraint.h>
#include <cosfitterexcept.h>

using namespace std;

/*!
  See equation (5) of Eisenstein '05, modified for non-flatness
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
    (\f$w\f$, \f$\Omega_m\f$,\f$\Omega_{DE}\f$)
  \returns 1/E term
*/
double auxconstraint::int_term_w( double z, void *params ) {
  double opz, w, om, ode, ok;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + z;
  w = in_params[0];
  om = in_params[1];
  ode = in_params[2];
  ok = 1.0 - om - ode;
  return 1.0/sqrt( (om*opz + ok)*opz*opz + ode*pow(opz,3.0*(1+w)) );
}

/*!
  See equation (5) of Eisenstein '05, modified for non-flatness
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
    (\f$w_0, w_a, \Omega_m, \Omega_{DE}\f$)
  \returns 1/E term
*/
double auxconstraint::int_term_w0_wa( double z, void *params ) {
  double opz, w0, wa, om, ode, ok;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + z;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  ok = 1.0 - om - ode;
  return 1.0/sqrt( (om*opz + ok)*opz*opz + 
		   ode*pow(opz,3.0*(1+w0+wa)) * exp(-3.0*wa*z/opz) );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters
    (\f$w_0\f$, \f$w_a\f$ \f$\Omega_m\f$,\f$\Omega_{DE}, atrans\f$)
  \returns 1/E term used in many cosmological expressions

  This uses the Komatsu et al. (2008) form for the expansion history
*/
double auxconstraint::int_term_w0_waKom( double z, void *params ) {
  double opz, w0, wa, om, ode, ok, atrans, val;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + z;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  atrans = in_params[4];
  ok = 1.0 - om - ode;
  double opat = 1.0 + atrans;
  val = opz*opz * ( opz*om + ok ) +
    ode * exp( -3.0*( z*wa/opz + (1+w0+opat*wa)*
                      log( (1/opz + atrans)/opat ) ) );
  return 1./sqrt(val);
}


/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
    (\f$\Omega_m\f$,\f$\Omega_{\Lambda}, \Omega_r\f$)
  \returns 1/E term used in many cosmological expressions
  This version includes \f$\Omega_r\f$
*/
double auxconstraint::one_over_E_rad( double z, void *params ) {
  double opz, om, ol, orad, ok;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + z;
  om = in_params[0];
  ol = in_params[1];
  orad = in_params[2];
  ok = 1.0 - om - ol - orad;
  return 1.0/sqrt( ((orad*opz+om)*opz + ok)*opz*opz + ol );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
    (\f$w\f$, \f$\Omega_m\f$,\f$\Omega_{DE}\f$,\f$\Omega_r\f$)
  \returns 1/E term used in many cosmological expressions
  This version includes \f$\Omega_r\f$
*/
double auxconstraint::one_over_E_w_rad( double z, void *params ) {
  double opz, w, om, ode, ok, orad;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + z;
  w = in_params[0];
  om = in_params[1];
  ode = in_params[2];
  orad = in_params[3];
  ok = 1.0 - om - ode - orad;
  return 1.0/sqrt( ((orad*opz+om)*opz + ok)*opz*opz + ode*pow(opz,3*(1+w)) );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
    (\f$w_0\f$, \f$w_a\f$ \f$\Omega_m\f$,\f$\Omega_{DE}\f$, \f$\Omega_r\f$)
  \returns 1/E term used in many cosmological expressions
  This version includes \f$\Omega_r\f$
*/
double auxconstraint::one_over_E_w0wa_rad( double z, void *params ) {
  double opz, w0, wa, om, ode, ok, orad;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + z;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  orad = in_params[4];
  ok = 1.0 - om - ode - orad;
  return 1.0/sqrt( ((orad*opz+om)*opz + ok)*opz*opz + 
		   ode*pow(opz,3*(1+w0+wa))*exp(-3 * wa * z/opz) );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters
    (\f$w_0\f$, \f$w_a\f$ \f$\Omega_m\f$,\f$\Omega_{DE}\f$, \f$\Omega_r,
    atrans\f$)
  \returns 1/E term used in many cosmological expressions
  This version includes \f$\Omega_r\f$, and is in the Komatsu form
  for w(a)
*/
double auxconstraint::one_over_E_w0waKom_rad( double z, void *params ) {
  double opz, w0, wa, om, ode, ok, orad, atrans,val;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + z;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  orad = in_params[4];
  atrans = in_params[5];
  ok = 1.0 - om - ode - orad;
  double opat = 1.0+atrans;
  val = ( (orad*opz + om)*opz + ok) * opz*opz +
    ode*exp( -3.0*( z*wa/opz + (1+w0+opat*wa)*log( (1/opz + atrans)/opat ) ) );
  return 1./sqrt(val);
}


/*!
  This is in the form for the GSL
  
  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
    (\f$\Omega_m , \Omega_{\Lambda} , \Omega_b h^2, h \f$)
  \returns The integral term for \f$r_s H_0 / c\f$  Does not include the
    \f$\sqrt(3)\f$ term, and is done in terms of \f$a\f$
*/
double auxconstraint::rsHovercint( double a, void *params ) {
  double om, ode, obh2, h, orad, ok;
  double *in_params = static_cast<double*>(params);

  om = in_params[0];
  ode = in_params[1];
  obh2 = in_params[2];
  h = in_params[3];

  //Includes 3.04 species of light Nuetrinos, value from PDG
  orad = oradh2 / (h*h); 

  ok = 1.0 - om - ode - orad;

  //a^2 * H/H0
  double a2E = sqrt( orad + a * (om + a*(ok + a*a*ode) ) );
  
  return 1.0 / ( a2E*sqrt(1.0 + 0.75*a*(obh2/ogammah2)) );
}

/*!
  This is in the form for the GSL
  
  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
    (\f$w, Omega_m , \Omega_{\Lambda} , \Omega_b h^2, h \f$)
  \returns The integral term for \f$r_s H_0 / c\f$  Does not include the
    \f$\sqrt(3)\f$ term, and is done in terms of \f$a\f$
*/
double auxconstraint::rsHovercint_w( double a, void *params ) {
  double w, om, ode, obh2, h, orad, ok;
  double *in_params = static_cast<double*>(params);
  
  w = in_params[0];
  om = in_params[1];
  ode = in_params[2];
  obh2 = in_params[3];
  h = in_params[4];

  //Includes 3.04 species of light Nuetrinos, value from PDG
  orad = oradh2 / (h*h); 

  ok = 1.0 - om - ode - orad;

  //a^2 * H/H0
  double a2E = sqrt( orad + a * (om + a*ok + ode*pow(a,-3.0*w)) );
  
  return 1.0 / ( a2E*sqrt(1.0 + 0.75*a*(obh2/ogammah2)) );
}

/*!
  This is in the form for the GSL
  
  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
    (\f$w_0, w_a, Omega_m , \Omega_{\Lambda} , \Omega_b h^2, h \f$)
  \returns The integral term for \f$r_s H_0 / c\f$  Does not include the
    \f$\sqrt(3)\f$ term, and is done in terms of \f$a\f$
*/
double auxconstraint::rsHovercint_w0wa( double a, void *params ) {
  double w0, wa, om, ode, obh2, h, orad, ok;
  double *in_params = static_cast<double*>(params);
  
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  obh2 = in_params[4];
  h = in_params[5];

  //Includes 3.04 species of light Nuetrinos, value from PDG
  orad = oradh2 / (h*h); 

  ok = 1.0 - om - ode - orad;

  //a^2 * H/H0
  double a2E = sqrt( orad + a * (om + a*ok + 
				ode*pow(a,-3.0*(w0+wa))*exp(-3.0*wa*(1-a)) ) );
  
  return 1.0 / ( a2E*sqrt(1.0 + 0.75*a*(obh2/ogammah2)) );
}

/*!
  This is in the form for the GSL

  \param[in] a Scale factor
  \param[in] params Cosmological parameters
    (\f$w_0, w_a, Omega_m , \Omega_{\Lambda} , \Omega_b h^2, h, atrans \f$)
  \returns The integral term for \f$r_s H_0 / c\f$  Does not include the
    \f$\sqrt(3)\f$ term, and is done in terms of \f$a\f$

  Uses the Komatsu et al. (2008) form for w(a)
*/
double auxconstraint::rsHovercint_w0waKom( double a, void *params ) {
  double w0, wa, om, ode, obh2, h, orad,  ok, atrans, opat;
  double *in_params = static_cast<double*>(params);

  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  obh2 = in_params[4];
  h = in_params[5];
  atrans = in_params[6];

  opat = 1.0 + atrans;

  //Includes 3.04 species of light Nuetrinos, value from PDG
  orad = oradh2 / (h*h);

  ok = 1.0 - om - ode - orad;

  //a^2 * H/H0
  double val = ode*exp( -3.0*( (1.0-a)*wa + (1+w0+opat*wa) *
                              log( (a + atrans)/opat ) ) );
  double a2E = sqrt( orad + a*(om + a*( ok + a*a*val ) ) );

  return 1.0 / ( a2E*sqrt(1.0 + 0.75*a*(obh2/ogammah2)) );
}


/*!
  Based on formulae from 
 <A HREF="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1992ARA%26A..30..499C&amp;db_key=AST&amp;data_type=HTML&amp;format=&amp;high=4270f9326909263">Carroll and Press, AARA 1993, 30: 499.</A>
  Ignores \f$\Omega_{rad}\f$.  This is only valid for
  the cosmological constant case (\f$w=-1\f$).
*/
bool auxconstraint::isRebounding(double om, double ol) {

  const double minom = 1e-4; //Smaller than this, use the small \om expansion
  if (om < minom) {
    return ol > 1.0 - om;
  }

  //Otherwise do the full calculation
  double arg = (1.0 - om)/om;
  if ( om >= 0.5 ) {
    return ol >= 4*om*pow( cos( acos( arg ) / 3.0 ), 3);
  } else {
    // cmath doesn't include inverse hyperbolic functions, so this
    // looks funny, but cosh^{-1}(z) = ln( z + sqrt(z+1)*sqrt(z-1) )
    return ol >= 4*om*pow( cosh( log( arg + sqrt(arg*arg-1) )/ 3.0),3);
  }
}


////////////////////////////////////////////////////////////////
//                   Distance Helper
////////////////////////////////////////////////////////////////
auxconstraint::distance_helper::distance_helper() :
  useKomatsuForm(false), atrans(0.1) {
  work = gsl_integration_workspace_alloc(1000); 
}

auxconstraint::distance_helper::~distance_helper() {
  gsl_integration_workspace_free(work);
}

/*!
  \param[in] z Redshift
  \param[in] obh2 \f$ \Omega_b h^2 \f$
  \param[in] om \f$ \Omega_m \f$
  \param[in] ode \f$ \Omega_{DE}\f$
  \param[in] w0 \f$w_{0}\f$
  \param[in] wa \f$w_{a}\f$
  \param[in] h \f$H_{0} = 100 h\f$ km/sec
  \param[out] success 1 if computation was successful, 0 if failed
  \returns \f$K_1 = H_{0} r_{s} \left( z \right) / c\f$

*/
double
auxconstraint::distance_helper::GetK1( double z, double obh2, double om,
				       double ode, double w0, double wa, 
				       double h, bool& success ) const {
  const double doubletol = 1e-5; //Used for various comso parameter tolerances
  const double prefac = 1.0 / sqrt(3.0);
  const double badval = 1e30;

  success = 0;

  double a = 1.0 / (1.0 + z);

  //Quick exit tests for bad inputs (bouncing, etc.)
  bool constw = fabs(wa) < doubletol;
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);
  if (cosmoconst && isRebounding( om, ode ) ) return badval;
  if (om < 0.0) return badval;

  //Set up integral
  gsl_function F;
  if ( cosmoconst ) {
    double params[4];
    params[0] = om;
    params[1] = ode;
    params[2] = obh2;
    params[3] = h;
    F.function = rsHovercint;
    F.params = &params;
  } else if (constw) {
    double params[5];
    params[0] = w0;
    params[1] = om;
    params[2] = ode;
    params[3] = obh2;
    params[4] = h;
    F.function = rsHovercint_w;
    F.params = &params;
  } else if (useKomatsuForm) {
    //w(a), Komatsu form
    double params[7];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = obh2;
    params[5] = h;
    params[6] = atrans;
    F.function = rsHovercint_w0waKom;
    F.params = &params;
  } else {
    //w(a), Linder form
    double params[6];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = obh2;
    params[5] = h;
    F.function = rsHovercint_w0wa;
    F.params = &params;
  }
  
  //This integral is more difficult than the BAO one, 
  // and requires adaptive integration
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  double intval, error;
  int st = gsl_integration_qag(&F,0,a,0,1e-7,1000,6,work,&intval,&error);
  gsl_set_error_handler( old_handler );

  //Failure if st != 0
  if (st) return badval;

  success = 1;
  return prefac * intval;
}


/*!
  \param[in] z Redshift
  \param[in] om \f$ \Omega_m \f$
  \param[in] ode \f$ \Omega_{DE}\f$
  \param[in] orad \f$ \Omega_r\f$, the density of radiation (not just photons)
  \param[in] w0 \f$w_{0}\f$
  \param[in] wa \f$w_{a}\f$
  \param[out] success 1 if computation was successful, 0 if failed
  \returns \f$K_2 = \left(1 + z \right) H_{0} 
  D_{A} \left( z_{\star} \right) / c\f$, where \f$D_{A}\f$ is
  the comoving angular diameter distance
*/
double
auxconstraint::distance_helper::GetK2( double z, double om,
				       double ode, double orad, double w0, 
				       double wa, bool& success ) const {
  const double doubletol = 1e-5; //Used for various comso parameter tolerance
  const double pi = 3.14159265358979323846264338327950288419716939937510582;
  const double badval = 1e30;

  success = 0;
  //Quick exit tests for bad inputs (bouncing, etc.)
  bool constw = fabs(wa) < doubletol;
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);
  if (cosmoconst && isRebounding( om, ode ) ) return badval;
  if (om < 0.0) return badval;

  //Set up integral
  gsl_function F;
  if ( cosmoconst ) {
    double params[3];
    params[0] = om;
    params[1] = ode;
    params[2] = orad;
    F.function = one_over_E_rad;
    F.params = &params;
  } else if (constw) {
    double params[4];
    params[0] = w0;
    params[1] = om;
    params[2] = ode;
    params[3] = orad;
    F.function = one_over_E_w_rad;
    F.params = &params;
  } else if (useKomatsuForm) {
    //w(a), Komatsu Form
    double params[6];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = orad;
    params[5] = atrans;
    F.function = one_over_E_w0waKom_rad;
    F.params = &params;
  } else {
    //w(a), Linder form
    double params[5];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = orad;
    F.function = one_over_E_w0wa_rad;
    F.params = &params;
  }
  
  //This integral is more difficult than the BAO one, 
  // and requires adaptive integration
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  double intval, error;
  int st = gsl_integration_qag(&F,0,z,0,1e-7,1000,6,work,&intval,&error);
  gsl_set_error_handler( old_handler );

  //Failure if st != 0
  if (st) return badval;
#if USEMKL
  if ( isnan(intval) ) return badval;
#else
  if (std::isnan(intval) ) return badval;
#endif

  double sqrtk, ok;
  double finalval;
  ok = 1.0 - om - ode - orad;
  if (fabs(ok) < doubletol) {
    finalval = intval;
  } else if (ok > 0.0) {
    //Open
    sqrtk = sqrt( fabs( ok ) );
    finalval = sinh( sqrtk * intval )/sqrtk;
  } else {
    //Closed
    sqrtk = sqrt( fabs( ok ) );

    //proper distance is zero or negative at this redshift
    double compval = pi / sqrtk;
    if (intval >= compval) return badval;

    finalval = sin( sqrtk * intval )/sqrtk;
  }

  success = 1;
  return finalval;
}

/*!
  \param[in] z Redshift
  \param[in] om \f$ \Omega_m \f$
  \param[in] ode \f$ \Omega_{DE}\f$
  \param[in] orad \f$ \Omega_r\f$, the density of radiation (not just photons)
  \param[in] w0 \f$w_{0}\f$
  \param[in] wa \f$w_{a}\f$
  \param[out] success 1 if computation was successful, 0 if failed
  \returns \f$1/E = H_{0} / H\left( z \right)\f$
*/
double
auxconstraint::distance_helper::GetOneOverE( double z, double om,
					     double ode, double orad, 
					     double w0, double wa, 
					     bool& success ) const {
  const double doubletol = 1e-5; //Used for various comso parameter tolerance
  const double badval = 1e30;

  success = 0;
  //Quick exit tests for bad inputs (bouncing, etc.)
  bool constw = fabs(wa) < doubletol;
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);
  if (cosmoconst && isRebounding( om, ode ) ) return badval;
  if (om < 0.0) return badval;

  //Set up integral
  if ( cosmoconst ) {
    double params[3];
    params[0] = om;
    params[1] = ode;
    params[2] = orad;
    success = 1;
    return one_over_E_rad( z, params );
  } else if (constw) {
    double params[4];
    params[0] = w0;
    params[1] = om;
    params[2] = ode;
    params[3] = orad;
    success = 1;
    return one_over_E_w_rad( z, params );
  } else if (useKomatsuForm) {
    //w(a), Komatsu Form
    double params[6];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = orad;
    params[5] = atrans;
    success = 1;
    return one_over_E_w0waKom_rad( z, params );
  } else {
    //w(a), Linder form
    double params[5];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = orad;
    success = 1;
    return one_over_E_w0wa_rad( z, params );
  }

  return badval;
}



/*!
  \param[in] obh2  \f$\Omega_b h^2\f$
  \param[in] omh2  \f$\Omega_m h^2\f$
  \returns The value of \f$z_{\star}\f$, the decopuling redshift

  Fitting forumula from Hu & Sugiyama (1996)
*/
double
auxconstraint::distance_helper::GetZstar( double obh2, double omh2 ) const {
  double g1 = 0.0783*pow(obh2,-0.238) / ( 1.0 + 39.5*pow(obh2,0.763) );
  double g2 = 0.560 / ( 1.0 + 21.1 * pow(obh2,1.81)  );
  return 1048.0*(1 + 0.00124*pow(obh2, -0.738))*(1+g1*pow(omh2,g2));
}

/*!
  \param[in] obh2  \f$\Omega_b h^2\f$
  \param[in] omh2  \f$\Omega_m h^2\f$
  \returns The value of \f$z_{d}\f$, the drag redshift

  Fitting formula from Eisenstein & Hu (1998)
*/
double
auxconstraint::distance_helper::GetZdrag( double obh2, double omh2 ) const {
  double b1 = 0.313*pow(omh2,-0.419)*(1.0 + 0.607*pow(omh2,0.674));
  double b2 = 0.238*pow(omh2,0.223);
  return 1291.0*pow(omh2,0.251) * (1.0 + b1*pow(obh2,b2)) / 
    (1.0+0.659*pow(omh2,0.828));
}



/*!
  \param[in] v Vector of \f$h, \Omega_b h^2\f$
  \param[in] params Other parameters, very indirectly
  \returns \f$\chi^2\f$
*/
double auxconstraint::wmap7_shift_chisq(const gsl_vector *v, void *params) {

  //Values for weak H prior to avoid edge effects
  const double minh = 0.2, maxh = 1.0, varh = 0.01;

  const double pi = 3.14159265358979323846264338327950288419716939937510582;

  //Values from Komatsu et al. (2010)
  const double wmap7_params[3] = { 302.09, 1.725, 1091.3 };

  //Return on failure
  const double badval = 1e50;

  double h, h2, obh2;
  h = gsl_vector_get(v,0);
  obh2 = gsl_vector_get(v,1);
  h2 = h*h;
  
  void **vptr = static_cast<void**>(params);

  distance_helper *const dh = static_cast<distance_helper *const>(vptr[0]);
  
  double* in_params = static_cast<double*>(vptr[1]);
  double om, ode, w0, wa;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  
  double omh2;
  omh2 = om * h2;

  double orad;
  orad = oradh2 / h2;

  //Bad value test
  if ( om < 0 || obh2 < 0 || oradh2 < 0 || h < 0 || obh2/h2 > om ) 
    return badval;

  double zs;
  zs = dh->GetZstar( obh2, omh2 );

  double K1, K2;
  bool success_K1, success_K2;
  K1=dh->GetK1( zs, obh2, om, ode, w0, wa, h, success_K1);
  if (! success_K1) return badval; //Failure
  if (K1 == 0) return badval; //Failure
  K2=dh->GetK2( zs, om, ode, orad, w0, wa, success_K2);
  if (! success_K2) return badval; //Failure

  double la, R;
  la = pi * K2 / K1;
  R  = sqrt( om )*K2;

  double delta_params[3]; //Between params and WMAP fits
  delta_params[0] = la - wmap7_params[0];
  delta_params[1] = R - wmap7_params[1];
  delta_params[2] = zs - wmap7_params[2];

  double invcov[3][3];
  invcov[0][0] = 2.305;
  invcov[0][1] = invcov[1][0] = 29.698;
  invcov[0][2] = invcov[2][0] = -1.333;
  invcov[1][1] = 6825.270;
  invcov[1][2] = invcov[2][1] = -113.180;
  invcov[2][2] = 3.414;

  double chisq;
  chisq = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      chisq += delta_params[i] * invcov[i][j] * delta_params[j];

  //Add weak H prior
  if ( h < minh ) chisq += ( minh - h ) * ( minh - h ) / varh;
  if ( h > maxh ) chisq += ( h - maxh ) * ( h - maxh ) / varh;

  return chisq;
}

/*!
  \param[in] v Vector of \f$h, \Omega_b h^2\f$
  \param[in] params Other parameters, very indirectly
  \returns \f$\chi^2\f$
*/
double auxconstraint::percival_chisq(const gsl_vector *v, void *params) {

  const unsigned int nz = 2;
  const double percival_z[nz] = { 0.2, 0.35 };
  const double percival_params[nz] = { 0.1905, 0.1097 };
  const double badval = 1e50;

  double h, h2, obh2;
  h = gsl_vector_get(v,0);
  obh2 = gsl_vector_get(v,1);
  h2 = h*h;
  
  void **vptr = static_cast<void**>(params);

  distance_helper *const dh = static_cast<distance_helper *const>(vptr[0]);
  
  double* in_params = static_cast<double*>(vptr[1]);
  double om, ode, w0, wa;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  
  double omh2;
  omh2 = om * h2;

  double orad;
  orad = oradh2 / h2;

  //Bad value test
  if ( om < 0 || obh2 < 0 || h < 0 || obh2 / h2 > om) 
    return badval;

  double zdrag;
  zdrag = dh->GetZdrag( obh2, omh2 );

  //Get H0 rs(zd)/c
  double K1, K2_0, K2_1;
  bool success_K1, success_K2;
  K1=dh->GetK1( zdrag, obh2, om, ode, w0, wa, h, success_K1 );
  if (! success_K1) return badval; //Failure

  //We need two values of K2, one for each z
  K2_0=dh->GetK2( percival_z[0], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_0 == 0) return badval; //Failure
  K2_1=dh->GetK2( percival_z[1], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_1 == 0) return badval; //Failure

  //Really, we hold H0/c CV = [K2^2 z/E]^1/3
  double DV_0, DV_1;
  bool success_DV;
  DV_0 = pow( K2_0 * K2_0 * percival_z[0] * 
	      dh->GetOneOverE( percival_z[0], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_0 == 0 ) return badval;  //Failure
  DV_1 = pow( K2_1 * K2_1 * percival_z[1] * 
	      dh->GetOneOverE( percival_z[1], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_1 == 0 ) return badval;  //Failure
  
  double delta_params[2]; //Difference between values and Percival fits
  delta_params[0] = K1/DV_0 - percival_params[0];
  delta_params[1] = K1/DV_1 - percival_params[1];

  double invcov[2][2];
  invcov[0][0] = 30124.0;
  invcov[0][1] = invcov[1][0] = -17227.0;
  invcov[1][1] = 86977.0;

  double chisq;
  chisq = 0.0;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      chisq += delta_params[i] * invcov[i][j] * delta_params[j];

  return chisq;
}

/*!
  \param[in] v Vector of \f$h, \Omega_b h^2\f$
  \param[in] params Other parameters, very indirectly
  \returns \f$\chi^2\f$ of Percival BAO and WMAP7
*/
double auxconstraint::percival_wmap7_chisq(const gsl_vector *v, void *params) {

  const unsigned int nz = 2;
  const double percival_z[nz] = { 0.2, 0.35 };
  const double percival_params[nz] = { 0.1905, 0.1097 };
  const double pi = 3.14159265358979323846264338327950288419716939937510582;
  const double wmap7_params[3] = { 302.09, 1.725, 1091.3 };
  const double badval = 1e50;

  double h, h2, obh2;
  h = gsl_vector_get(v,0);
  obh2 = gsl_vector_get(v,1);
  h2 = h*h;
  
  void **vptr = static_cast<void**>(params);

  distance_helper *const dh = static_cast<distance_helper *const>(vptr[0]);
  
  double* in_params = static_cast<double*>(vptr[1]);
  double om, ode, w0, wa;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  
  double omh2;
  omh2 = om * h2;

  double orad;
  orad = oradh2 / h2;

  //Bad value test
  if ( om < 0 || obh2 < 0 || h < 0 || obh2 / h2 > om) 
    return badval;

  double zdrag, zstar;
  zdrag = dh->GetZdrag( obh2, omh2 );
  zstar = dh->GetZstar( obh2, omh2 );

  //Get H0 rs(zd)/c
  double K1_wmap, K1_bao, K2_wmap, K2_bao_0, K2_bao_1;
  bool success_K1, success_K2;
  K1_bao=dh->GetK1( zdrag, obh2, om, ode, w0, wa, h, success_K1 );
  if (! success_K1) return badval; //Failure
  K1_wmap=dh->GetK1( zstar, obh2, om, ode, w0, wa, h, success_K1 );
  if (! success_K1 || K1_wmap == 0.0) return badval; //Failure

  K2_wmap=dh->GetK2( zstar, om, ode, orad, w0, wa, success_K2 );
  if (! success_K2) return badval; //Failure
  K2_bao_0=dh->GetK2( percival_z[0], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_bao_0 == 0) return badval; //Failure
  K2_bao_1=dh->GetK2( percival_z[1], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_bao_1 == 0) return badval; //Failure

  //WMAP shift params
  double la, R;
  la = pi * K2_wmap / K1_wmap;
  R  = sqrt( om )*K2_wmap;

  //Really, we hold H0/c CV = [K2^2 z/E]^1/3
  double DV_0, DV_1;
  bool success_DV;
  DV_0 = pow( K2_bao_0 * K2_bao_0 * percival_z[0] * 
	      dh->GetOneOverE( percival_z[0], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_0 == 0 ) return badval;  //Failure
  DV_1 = pow( K2_bao_1 * K2_bao_1 * percival_z[1] * 
	      dh->GetOneOverE( percival_z[1], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_1 == 0 ) return badval;  //Failure

  double delta_params_bao[2]; //Difference between values and Percival fits
  delta_params_bao[0] = K1_bao/DV_0 - percival_params[0];
  delta_params_bao[1] = K1_bao/DV_1 - percival_params[1];

  double invcov_bao[2][2];
  invcov_bao[0][0] = 30124.0;
  invcov_bao[0][1] = invcov_bao[1][0] = -17227.0;
  invcov_bao[1][1] = 86977.0;

  double chisq;
  chisq = 0.0;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      chisq += delta_params_bao[i] * invcov_bao[i][j] * delta_params_bao[j];

  double delta_params_wmap[3]; //Between params and WMAP fits
  delta_params_wmap[0] = la - wmap7_params[0];
  delta_params_wmap[1] = R - wmap7_params[1];
  delta_params_wmap[2] = zstar - wmap7_params[2];

  //Table 10 of comatsu
  double invcov_wmap[3][3];
  invcov_wmap[0][0] = 2.305;
  invcov_wmap[0][1] = invcov_wmap[1][0] = 29.698;
  invcov_wmap[0][2] = invcov_wmap[2][0] = -1.333;
  invcov_wmap[1][1] = 6825.270;
  invcov_wmap[1][2] = invcov_wmap[2][1] = -113.180;
  invcov_wmap[2][2] = 3.414;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      chisq += delta_params_wmap[i] * invcov_wmap[i][j] * 
	delta_params_wmap[j];

  return chisq;
}


///////////////////////////////////////////////////////////////////////

auxconstraint::sdss_lrg_bao::sdss_lrg_bao(utility::cosmo_fittype ft, 
			   double n, double f) :
  fittype(ft), ns(n), fnu(f), errdef(1.0), useKomatsuForm(false), atrans(0.01)
{
}

/*!
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \param[in] w    \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \returns Eisenstein A parameter

  Uses the value of A given in the paper, possibly modified for
  the neutrino mass fraction as in Goobar et al. (2006).
  Note that an expansion is used for the neutrino fraction, so if
  the value is large the expansion is probably invalid
 */
double 
auxconstraint::sdss_lrg_bao::aval(double om, double ode, double w,
				  double wa) const {

  const double doubletol = 1e-5; //Used for various comso parameter tolerances

  double zval = 0.35; //Eisenstein value

  double ok = 1.0 - om - ode;
  if (om < 0.0) return 10000;

  bool cosmoconst = (!useKomatsuForm) && (fabs(w + 1) < doubletol) 
    && ( fabs(wa) < doubletol );

  //Bouncing universe test -- hard to have the BAO constraint,
  // which assumes things about the CMB
  if (cosmoconst && isRebounding( om, ode ) ) return 1e20;

  //Set up usage of GSL integrator
  gsl_function F;
  size_t neval;
  double intval, error, termval;
  int st;
  if (useKomatsuForm) {
    double params[5];
    params[0] = w; params[1] = wa, params[2] = om; params[3] = ode;
    params[4] = atrans;
    termval = int_term_w0_waKom( zval, params );
    F.function = int_term_w0_waKom;
    F.params = &params;
    st = gsl_integration_qng(&F,0,zval,0,1e-7,&intval,&error,&neval);
  } else {
    double params[4];
    params[0] = w; params[1] = wa, params[2] = om; params[3] = ode;
    termval = int_term_w0_wa( zval, params );
    F.function = int_term_w0_wa;
    F.params = &params;
    st = gsl_integration_qng(&F,0,zval,0,1e-7,&intval,&error,&neval);
  }

  double curra;
  if (st) {
    //Integral didn't evaluate
    return 10000;
  } else {
    //Recall that f1 is 1/E, not E from Eisenstein
    if ( fabs(ok) < doubletol ) {
      //Treat as flat
      curra = sqrt( om ) * pow( termval, 1.0/3.0 ) * 
	pow( intval / zval, 2.0/3.0 );
    } else if ( ok < 0.0 ) {
      //Closed Universe
      curra = sqrt( om ) * 
	pow( termval / fabs( ok ), 1.0/3.0 ) *
	pow( 1.0/zval * sin( sqrt(fabs(ok)) * intval ), 2.0/3.0 );
    } else {
      //Open Universe
      curra = sqrt( om ) * 
	pow( termval / ok,1.0/3.0 ) *
	pow( 1.0/zval * sinh( sqrt(ok) * intval ), 2.0/3.0 );
    }
  }
  return curra;
}

/*!
  \param[in] par  Params in order \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \returns Eisenstein A at parameter values
 */
double 
auxconstraint::sdss_lrg_bao::aval(const std::vector<double>& par) const {
  double om, ode, w, wa, alpha, beta;
  utility::parse_param_values_nosm( par, fittype, om, ode, w, wa, alpha, 
				    beta );
  return aval(om,ode,w,wa);
}

/*!
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \param[in] w   \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \returns \f$\chi^2\f$ at parameter values

  Uses the value of A given in the paper, possibly modified for
  the neutrino mass fraction as in Goobar et al. (2006).
  Note that an expansion is used for the neutrino fraction, so if
  the value is large the expansion is probably invalid
 */
double 
auxconstraint::sdss_lrg_bao::operator()(double om, double ode, 
					double w, double wa) const 
{

  const double doubletol = 1e-5; //Used for various comso parameter tolerances

  double eis_aval = 0.469; //Eisenstein value
  double eis_aval_error = 0.017; //Eisenstein value

  //Modify for different scalar index
  if (fabs(ns - 0.98) > doubletol) eis_aval *= pow( ns/0.98, -0.35 );

  //And neutrino mass
  if (fabs(fnu) > doubletol) eis_aval *= (1.0 + 0.94*fnu);

  double chisq = ( eis_aval-aval(om,ode,w,wa) )/eis_aval_error;
  return chisq*chisq;
}
/*!
  \param[in] par  Params in order \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \returns \f$\chi^2\f$ at parameter values
 */
double 
auxconstraint::sdss_lrg_bao::operator()(const std::vector<double>& par) const {
  double om, ode, w, wa, alpha, beta;
  utility::parse_param_values_nosm( par, fittype, om, ode, w, wa, alpha, 
				    beta );
  return (*this)(om,ode,w,wa);
}

////////////////////////////////////////////////////////////////////////
//                         bao_P09                                    //
////////////////////////////////////////////////////////////////////////

auxconstraint::bao_P09::bao_P09(utility::cosmo_fittype ft) :
  fittype(ft), errdef(1.0)
{
  v = gsl_vector_alloc(2);
  ss = gsl_vector_alloc(2);
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc(T, 2);
}

auxconstraint::bao_P09::~bao_P09() {
  gsl_vector_free(v);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
}

double 
auxconstraint::bao_P09::operator()(double om, double ode, 
				       double w0, double wa) const
{
  const size_t maxiter = 1000;
  const double badval = 1e30;
  const double sizetol1 = 1e-3, sizetol2 = 1e-4;

  if (om <= 0) return badval;  //Can't handle no matter universe

  //Do an initial evalution of K2 to make sure we aren't in a bad region
  // using default values.  Basically, this is to test for bad closed
  // universes with 0 or negative D_A
  bool K2_test_success;
  dhelp.GetK2( 1090.0, om, ode, 5e-5, w0, wa, K2_test_success );
  if (!K2_test_success) return badval;

  //Set initial conditions for h, obh2 to WMAP7 values
  gsl_vector_set(v, 0, 0.724);
  gsl_vector_set(v, 1, 0.02268);

  //Set initial step sizes
  gsl_vector_set(ss, 0, 0.05);
  gsl_vector_set(ss, 1, 0.004);
  
  //Make par vector
  double params[4];
  params[0] = w0;
  params[1] = wa;
  params[2] = om;
  params[3] = ode;
  
  //We have to jump through some serious hoops to pass
  // in an instance of the class as a parameter
  void **varr;
  varr = new void*[2];
  varr[0] = static_cast<void*>(&dhelp); //This is evil, but such is c++->c
  varr[1] = static_cast<void*>(params);
  void *vptr;
  vptr = static_cast<void*>(varr);

  //Setup function
  minex_func.n = 2;
  minex_func.f = &percival_chisq;
  minex_func.params = vptr;

  gsl_multimin_fminimizer_set(s, &minex_func, v, ss);

  int status;
  size_t iter = 0;
  double size;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
           
    if (status) 
      break;
     
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, sizetol1);
  } while (status == GSL_CONTINUE && iter < maxiter);

  //Give large chisq if minimz didn't converge
  if (status != GSL_SUCCESS) {
    return badval;
  } else {
    //Following the recommendation of NumRec, if we succeeded
    // we restart the minimization at the current guess

    //reset initial step sizes to half old values
    gsl_vector_set(ss, 0, 0.025);
    gsl_vector_set(ss, 1, 0.002);

    gsl_multimin_fminimizer_set(s, &minex_func, s->x, ss);
    iter = static_cast<size_t>(0);
    do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status) 
	break;
      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, sizetol2);
    } while (status == GSL_CONTINUE && iter < maxiter);
    //Here we use even if the second didn't converge
  }

  //Get chisq at minimum
  double chisq;
  chisq = percival_chisq( s->x, varr );

  //Clean up
  delete[] varr;

  return chisq;

}

/*!
  \param[in] par  Params in order \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \returns \f$\chi^2\f$ at parameter values
 */
double 
auxconstraint::bao_P09::operator()(const std::vector<double>& par) const {
  double om, ode, w, wa, alpha, beta;
  utility::parse_param_values_nosm( par, fittype, om, ode, w, wa, alpha, 
				    beta );
  return (*this)(om,ode,w,wa);
}

////////////////////////////////////////////////////////////////////////
//                         wmap7yr_dls                                //
////////////////////////////////////////////////////////////////////////

auxconstraint::wmap7yr_dls::wmap7yr_dls(utility::cosmo_fittype ft) :
  fittype(ft), errdef(1.0)
{
  v = gsl_vector_alloc(2);
  ss = gsl_vector_alloc(2);
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc(T, 2);
}

auxconstraint::wmap7yr_dls::~wmap7yr_dls() {
  gsl_vector_free(v);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
}

double 
auxconstraint::wmap7yr_dls::operator()(double om, double ode, 
				       double w0, double wa) const
{
  const size_t maxiter = 500;
  const double badval = 1e30;
  const double sizetol1 = 1e-4, sizetol2 = 1e-6;
  if (om <= 0) return badval;  //Can't handle no matter universe

  //Do an initial evalution of K2 to make sure we aren't in a bad region
  // using default values.  Basically, this is to test for bad closed
  // universes with 0 or negative D_A
  bool K2_test_success;
  dhelp.GetK2( 1090.0, om, ode, 5e-5, w0, wa, K2_test_success );
  if (!K2_test_success) return badval;

  //Set initial conditions for h, obh2 to WMAP7 values
  gsl_vector_set(v, 0, 0.724);
  gsl_vector_set(v, 1, 0.02268);

  //Set initial step sizes
  gsl_vector_set(ss, 0, 0.05);
  gsl_vector_set(ss, 1, 0.004);
  
  //Make par vector
  double params[4];
  params[0] = w0;
  params[1] = wa;
  params[2] = om;
  params[3] = ode;
  
  //We have to jump through some serious hoops to pass
  // in an instance of the class as a parameter
  void **varr;
  varr = new void*[2];
  varr[0] = static_cast<void*>(&dhelp); //This is evil, but such is c++->c
  varr[1] = static_cast<void*>(params);
  void *vptr;
  vptr = static_cast<void*>(varr);

  //Setup function
  minex_func.n = 2;
  minex_func.f = &wmap7_shift_chisq;
  minex_func.params = vptr;

  gsl_multimin_fminimizer_set(s, &minex_func, v, ss);

  int status;
  size_t iter = 0;
  double size;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
           
    if (status) 
      break;
     
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, sizetol1);
  } while (status == GSL_CONTINUE && iter < maxiter);

  //Give large chisq if minimz didn't converge
  if (status != GSL_SUCCESS) {
    return badval;
  } else {
    //Following the recommendation of NumRec, if we succeeded
    // we restart the minimization at the current guess

    //reset initial step sizes to half old values
    gsl_vector_set(ss, 0, 0.025);
    gsl_vector_set(ss, 1, 0.002);

    gsl_multimin_fminimizer_set(s, &minex_func, s->x, ss);
    iter = static_cast<size_t>(0);
    do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status) 
	break;
      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, sizetol2);
    } while (status == GSL_CONTINUE && iter < maxiter);
    //Here we use even if the second didn't converge
  }

  //Get chisq at minimum
  double chisq;
  chisq = wmap7_shift_chisq( s->x, varr );

  //Clean up
  delete[] varr;

  return chisq;

}

/*!
  \param[in] par  Params in order \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \returns \f$\chi^2\f$ at parameter values
 */
double 
auxconstraint::wmap7yr_dls::operator()(const std::vector<double>& par) const {
  double om, ode, w, wa, alpha, beta;
  utility::parse_param_values_nosm( par, fittype, om, ode, w, wa, alpha, 
				    beta );
  return (*this)(om,ode,w,wa);
}

////////////////////////////////////////////////////////////////////////
//                         wmap7_bao_P09                              //
////////////////////////////////////////////////////////////////////////

auxconstraint::wmap7_bao_P09::wmap7_bao_P09(utility::cosmo_fittype ft) :
  fittype(ft), errdef(1.0)
{
  v = gsl_vector_alloc(2);
  ss = gsl_vector_alloc(2);
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc(T, 2);
}

auxconstraint::wmap7_bao_P09::~wmap7_bao_P09() {
  gsl_vector_free(v);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
}

double 
auxconstraint::wmap7_bao_P09::operator()(double om, double ode, 
				       double w0, double wa) const
{
  const size_t maxiter = 500;
  const double badval = 1e30;
  const double sizetol1 = 1e-4, sizetol2 = 1e-6;

  if (om <= 0) return badval;  //Can't handle no matter universe

  //Do an initial evalution of K2 to make sure we aren't in a bad region
  // using default values.  Basically, this is to test for bad closed
  // universes with 0 or negative D_A
  bool K2_test_success;
  dhelp.GetK2( 1090.0, om, ode, 5e-5, w0, wa, K2_test_success );
  if (!K2_test_success) return badval;

  //Set initial conditions for h, obh2 to WMAP7 values
  gsl_vector_set(v, 0, 0.724);
  gsl_vector_set(v, 1, 0.02268);

  //Set initial step sizes
  gsl_vector_set(ss, 0, 0.05);
  gsl_vector_set(ss, 1, 0.004);
  
  //Make par vector
  double params[4];
  params[0] = w0;
  params[1] = wa;
  params[2] = om;
  params[3] = ode;
  
  //We have to jump through some serious hoops to pass
  // in an instance of the class as a parameter
  void **varr;
  varr = new void*[2];
  varr[0] = static_cast<void*>(&dhelp); //This is evil, but such is c++->c
  varr[1] = static_cast<void*>(params);
  void *vptr;
  vptr = static_cast<void*>(varr);

  //Setup function
  minex_func.n = 2;
  minex_func.f = &percival_wmap7_chisq;
  minex_func.params = vptr;

  gsl_multimin_fminimizer_set(s, &minex_func, v, ss);

  int status;
  size_t iter = 0;
  double size;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
           
    if (status) 
      break;
     
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, sizetol1);
  } while (status == GSL_CONTINUE && iter < maxiter);

  //Give large chisq if minimz didn't converge
  if (status != GSL_SUCCESS) {
    return badval;
  } else {
    //Following the recommendation of NumRec, if we succeeded
    // we restart the minimization at the current guess

    //reset initial step sizes to half old values
    gsl_vector_set(ss, 0, 0.025);
    gsl_vector_set(ss, 1, 0.002);

    gsl_multimin_fminimizer_set(s, &minex_func, s->x, ss);
    iter = static_cast<size_t>(0);
    do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status) 
	break;
      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, sizetol2);
    } while (status == GSL_CONTINUE && iter < maxiter);
    //Here we use even if the second didn't converge
  }

  //Get chisq at minimum
  double chisq;
  chisq = percival_wmap7_chisq( s->x, varr );

  //Clean up
  delete[] varr;

  return chisq;

}

/*!
  \param[in] par  Params in order \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \returns \f$\chi^2\f$ at parameter values
 */
double 
auxconstraint::wmap7_bao_P09::operator()(const std::vector<double>& par) const {
  double om, ode, w, wa, alpha, beta;
  utility::parse_param_values_nosm( par, fittype, om, ode, w, wa, alpha, 
				    beta );
  return (*this)(om,ode,w,wa);
}
