#include<Minuit2/FunctionMinimum.h>
#include<Minuit2/MnUserParameterState.h>
#include<Minuit2/MnPrint.h>
#include<Minuit2/MnMigrad.h>
#include<Minuit2/MnMinos.h>
#include<Minuit2/MnContours.h>
#include<Minuit2/MnMachinePrecision.h>
#include<Minuit2/MnPlot.h>
#include<Minuit2/MinosError.h>
#include<Minuit2/MnUserCovariance.h>
#include<Minuit2/ContoursError.h> 

#include<iostream>
#include<iomanip>
#include<string>
#include<sstream>
#include<vector>
#include<map>
#include<getopt.h>

#include <chifunc.h>
#include <auxconstraint.h>
#include <cosfitterexcept.h>
#include <utility.h>

using namespace ROOT::Minuit2;

int main(int argc, char **argv) {

  bool fixalphaerrset, fixbetaerrset, fixalpha, fixbeta, verbose;
  bool flatfit, wfit, wafit, includebao, includewmap, docontours;
  bool have_mag_cov, have_width_cov, have_colour_cov;
  bool have_magwidth_cov, have_magcolour_cov, have_widthcolour_cov;
  bool useKomatsuForm, useEisenstein;
  std::string mag_cov_file, width_cov_file, colour_cov_file;
  std::string magwidth_cov_file, magcolour_cov_file, widthcolour_cov_file;
  std::string datafile;
  double pecz, alpha, beta, fixalphaerrval, fixbetaerrval;
  double equiv_width_cut;
  double atrans, om, ode, w0, wa;
  bool fixom, fixode, fixw0, fixwa;
  std::map<unsigned int,double> intrinsicdisp;
  int ncontours;

  //Defaults
  docontours = false;
  fixalpha = fixbeta = false;
  verbose = false;
  fixalphaerrset = false; alpha = fixalphaerrval = 1.5;
  fixbetaerrset = false; beta = fixbetaerrval = 2.3;
  have_mag_cov=false;have_width_cov=false;have_colour_cov=false;
  have_magwidth_cov=false;have_magcolour_cov=false;have_widthcolour_cov=false;
  flatfit = false;
  wfit = false;
  wafit = false;
  useKomatsuForm = false;
  useEisenstein = false;
  atrans = 0.01;
  includebao = false;
  includewmap = false;
  pecz = 0.001;
  intrinsicdisp[0] = 0.17;
  om = 0.25;
  ode = 0.75;
  w0 = -1.0;
  wa = 0.0;
  fixom = fixode = fixw0 = fixwa = false;
  ncontours = 40;
  equiv_width_cut = 0.0;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"contours",no_argument,0,'c'},
    {"flat",no_argument,0,'f'},
    {"wfit",no_argument,0,'w'},
    {"wafit",no_argument,0,'5'},
    {"usekomatsu",no_argument,0,'k'},
    {"verbose",no_argument,0,'v'},
    {"intrinsicdisp",required_argument,0,'i'},
    {"fixalphaerr",required_argument,0,'1'},
    {"fixbetaerr",required_argument,0,'2'},
    {"includebao",no_argument,0,'3'},
    {"includewmap",no_argument,0,'4'},
    {"useeisenstein",no_argument,0,'E'},
    {"alpha",required_argument,0,'a'},
    {"beta",required_argument,0,'b'},
    {"omega_m",required_argument,0,'&'},
    {"omega_de",required_argument,0,'*'},
    {"w0",required_argument,0,'('},
    {"wa",required_argument,0,')'},
    {"ncontours",required_argument,0,'n'},
    {"pecz",required_argument,0,'p'},
    {"atrans",required_argument,0,'6'},
    {"equiv_widthcut",required_argument,0,'e'},
    {"magcov",required_argument,0,'!'},
    {"widthcov",required_argument,0,'@'},
    {"colourcov",required_argument,0,'#'},
    {"magwidthcov",required_argument,0,'$'},
    {"magcolourcov",required_argument,0,'%'},
    {"widthcolourcov",required_argument,0,'^'},
    {"help", no_argument, 0, 'h'},
    {0,0,0,0}
  };

  while ( ( c = getopt_long(argc,argv,
			    "hcEfwv345ki:1:2:a:b:e:p:6:!:@:#:$:%:^:&:*:(:):n:",
			    long_options, &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tcosfitter-twoscriptm -- Performs a Minuit cosmology fit."
		<< std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tcosfitter-twoscriptm [ -i | --intrinsicdisp=DATASET IDISP ] " 
		<< std::endl;
      std::cerr << "\t [ -e, --equivwidthcut ] [ --omega_m OMEGA_M ] "
		<< "[ --omega_de OMEGA_DE ]" << std::endl;
      std::cerr << "\t [ --w0 w0 ] [ --wa WA ] [ -f, --flat ]"
		" [ -a, --alpha ALPHA ]" << std::endl;
      std::cerr << "\t [ -b, --beta BETA ] [ -w, --wfit ] [ --wafit ] "
		<< std::endl;
      std::cerr << "\t [ -k, --usekomatsu ] [ --atrans=ATRANS ]" << std::endl;
      std::cerr << "\t [ --includebao ] [ --includewmap ] "
		<< "[ -E, --useeisenstein ]" << std::endl;
      std::cerr << "\t [ -p, --pecz PECZ ] [ --fixalphaerr FIXAL ] "
		<< "[ --fixbetaerr FIXBETA ]" << std::endl;
      std::cerr << "\t [ --magcov=MAGCOVFILE ] [ --widthcov=WIDTHCOVFILE ]"
		<< std::endl;
      std::cerr << "\t [ --colourcov=COLOURCOVFILE ]"
		<< " [ --magwidthcov=MAGWIDTHCOVFILE ]" << std::endl;
      std::cerr << "\t [ --magcolourcov=MAGCOLOURCOVFILE ]" << std::endl;
      std::cerr << "\t [ --widthcolourcov=WIDTHCOLOURCOVFILE ]" 
		<< " [ -v, --verbose ]" << std::endl;
      std::cerr << "\t [ -c, --contours ] [ -n, --ncontours NCONTOURS ]"
		<< " datafile" << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEstimates the cosmological and nuisance parameters" <<
	" from a set of" << std::endl;
      std::cerr << "\tsupernova data.  The input datafile should be in the "
		<< std::endl;
      std::cerr << "\tsimple_cosfitter format" << std::endl;
      std::cerr << "OPTIONS:" << std::endl;
      std::cerr << "\t-f, --flat" << std::endl;
      std::cerr << "\t\tDo fit assuming flat Universe." << std::endl;
      std::cerr << "\t-w, --wfit" << std::endl;
      std::cerr << "\t\tFit for w." << std::endl;
      std::cerr << "\t--wafit" << std::endl;
      std::cerr << "\t\tFit for wa (requires a flat universe)." << std::endl;
      std::cerr << "\t-k, --usekomatsu" << std::endl;
      std::cerr << "\t\tUse the Komatsu et al. (2008) form for w(a), rather"
		<< std::endl;
      std::cerr << "\t\tthan the Linder (2003) one.  Ignored unless --wafit." 
		<< std::endl;
      std::cerr << "\t-i, --intrinsicdisp IDISP" << std::endl;
      std::cerr << "\t\tThe intrinsic dispersion, in mags.  This takes two"
		<< "arguments" << std::endl;
      std::cerr << "\t\tthe data set number followed by the value.  The " 
		<< "default " << std::endl;
      std::cerr << "\t\tis 0 0.17 (i.e., dataset 0 has 0.17 mag dispersion)" 
		<< std::endl;
      std::cerr << "\t-e, --equiv_widthcut EQUIV_WIDTHCUT" << std::endl;
      std::cerr << "\t\tCut on equivalent width parameter.  Values below this"
		<< std::endl;
      std::cerr << "\t\tare matched to scriptm1, values above this to"
		<< " scriptm2." << std::endl;
      std::cerr << "\t\t(def: 0.0)" << std::endl;
      std::cerr << "\t--omega_m OMEGA_M" << std::endl;
      std::cerr << "\t\tDon't fit for omega_m, but fix it at this value."
		<< std::endl;
      std::cerr << "\t--omega_de OMEGA_DE" << std::endl;
      std::cerr << "\t\tDon't fit for omega_de, but fix it at this value."
		<< std::endl;
      std::cerr << "\t--w0 W0" << std::endl;
      std::cerr << "\t\tDon't fit for w_0, but fix it at this value."
		<< std::endl;
      std::cerr << "\t--wa WA" << std::endl;
      std::cerr << "\t\tDon't fit for w_a, but fix it at this value."
		<< std::endl;
      std::cerr << "\t-a, --alpha ALPHA" << std::endl;
      std::cerr << "\t\tDon't fit for alpha, but instead fix it at the " <<
	std::endl;
      std::cerr << "\t\tgiven value." << std::endl;
      std::cerr << "\t-b, --beta BETA" << std::endl;
      std::cerr << "\t\tDon't fit for beta, but instead fix it at the " <<
	std::endl;
      std::cerr << "\t\tgiven value." << std::endl;
      std::cerr << "\t--fixalphaerr FIXAL" << std::endl;
      std::cerr << "\t\tFix alpha at the specified value for error " << 
	" propagation." << std::endl;
      std::cerr << "\t--fixbetaerr FIXBETA" << std::endl;
      std::cerr << "\t\tFix beta at the specified value for error " 
		<< " propagation." << std::endl;
      std::cerr << "\t-p, --pecz PECZ" << std::endl;
      std::cerr << "\t\tPeculiar velocity in redshift units. Default: 0.001"
		<< std::endl;
      std::cerr << "\t--atrans=ATRANS" << std::endl;
      std::cerr << "\t\tTransition scale factor if Komatsu form for w(a) is"
		<< std::endl;
      std::cerr << "\t\tbeing used (def: 0.01)." << std::endl;
      std::cerr << "\t--includebao" << std::endl;
      std::cerr << "\t\tInclude the Percival et al. (2009) Baryon Acoustic"
		<< std::endl;
      std::cerr << "\t\tOscillations prior.  Also see --useeisenstein"
		<< std::endl;
      std::cerr << "\t--includewmap" << std::endl;
      std::cerr << "\t\tInclude the Komatsu (2010) WMAP 7th year"
		<< std::endl;
      std::cerr << "\t\tdistance to last scattering constraints." << std::endl;
      std::cerr << "\t-E, --useeisenstein" << std::endl;
      std::cerr << "\t\tUse the Eisenstein '05 BAO constraints instead of"
		<< " Percival 09." << std::endl;
      std::cerr << "\t-c, --contours" << std::endl;
      std::cerr << "\t\tDraw (with primitive text graphics) contours." 
		<< std::endl;
      std::cerr << "\t-n, --ncontours NCONTOURS" << std::endl;
      std::cerr << "\t\tNumber of contours points to make (def: 40)."
		<< std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tActivate verbose mode." << std::endl;
      return 0;
      break;
    case 'a' :
      alpha = atof(optarg);
      fixalpha = true;
      break;
    case 'b' :
      beta = atof(optarg);
      fixbeta = true;
      break;					
    case 'c' :
      docontours = true;
      break;
    case 'n' :
      ncontours = atoi(optarg);
      break;
    case 'e' :
      equiv_width_cut = atof(optarg);
      break;
    case '1' :
      fixalphaerrset = true;
      fixalphaerrval = atof(optarg);
      break;
    case '2' :
      fixbetaerrset = true;
      fixbetaerrval = atof(optarg);
      break;
    case '3' :
      includebao = true;
      break;
    case 'E' :
      useEisenstein = true;
      break;
    case '4' :
      includewmap = true;
      break;
    case '5' :
      //Can't fit W(a) without also fitting w
      wfit = true;
      wafit = true;
      break;
    case 'k' :
      useKomatsuForm = true;
      break;
    case '6' :
      atrans = atof(optarg);
      break;
    case 'f' :
      flatfit = true;
      break;
    case 'w' :
      wfit = true;
      break;
    case 'i' :
      if (optind == argc) {
	std::cerr << "Missing second arg for -i" << std::endl;
	return 4;
      }
      intrinsicdisp[static_cast<unsigned int>(atoi( optarg ))] = 
	atof( argv[optind++] );
      break;
    case 'p' :
      pecz = atof(optarg);
      break;
    case '!' :
      have_mag_cov = true;
      mag_cov_file = std::string( optarg );
      break;
    case '@' :
      have_width_cov = true;
      width_cov_file = std::string( optarg );
      break;
    case '#' :
      have_colour_cov = true;
      colour_cov_file = std::string( optarg );
      break;
    case '$' :
      have_magwidth_cov = true;
      magwidth_cov_file = std::string( optarg );
      break;
    case '%' :
      have_magcolour_cov = true;
      magcolour_cov_file = std::string( optarg );
      break;
    case '^' :
      have_widthcolour_cov = true;
      widthcolour_cov_file = std::string( optarg );
      break;
    case '&' :
      fixom = true;
      om = atof(optarg);
      break;
    case '*' :
      fixode = true;
      flatfit = false; //This may be overwritten later.
      ode = atof(optarg);
      break;
    case '(' :
      fixw0 = true;
      wfit = true;
      w0 = atof(optarg);
      break;
    case ')' :
      fixwa = true;
      wafit = true;
      wa = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if ( optind == argc ) {
    std::cerr << "Required argument datafile not provided" << std::endl;
    return 1;
  }
  if (useEisenstein && (!includebao)) {
    std::cout << "Warning: You are trying to use Eisenstein '05 but haven't"
	      << " set --includebao" << std::endl;
    std::cout << "BAO constraints will not be included" << std::endl;
  }

  datafile = std::string( argv[optind] );

  if (verbose) {
    std::cout << "Running fit to datafile: " << datafile 
	      << std::endl;
    std::cout << "Intrinsic dispersion values: " << std::endl;
    for ( std::map<unsigned int, double>::const_iterator 
	    it=intrinsicdisp.begin(); it != intrinsicdisp.end(); ++it) 
      std::cout << " Dataset: " << it->first << " value: "
		<< it->second << std::endl;
  }
  
  SNeData sne;
  try {
    sne.setName("Supernova");
    sne.readData( datafile );

    if (have_mag_cov) sne.readMagCovData( mag_cov_file );
    if (have_width_cov) sne.readWidthCovData( width_cov_file );
    if (have_colour_cov) sne.readColourCovData( colour_cov_file );
    if (have_magwidth_cov) sne.readMagWidthCovData( magwidth_cov_file );
    if (have_magcolour_cov) sne.readMagColourCovData( magcolour_cov_file );
    if (have_widthcolour_cov) 
      sne.readWidthColourCovData( widthcolour_cov_file );
    sne.zcmbsort();

  } catch (const CosFitterExcept& ex) {
    std::cerr << "Error performing fit" << std::endl;
    std::cerr << ex << std::endl;
    std::cerr << "Aborting" << std::endl;
    return 2;
  }
  
  //Determine the fit type
  utility::cosmo_fittype ft;
  if (flatfit) {
    if (wfit) {
      if (wafit) ft = utility::flat_omegam_w0_wa; 
      else ft = utility::flat_omegam_w;
    } else ft = utility::flat_omegam;
  } else if (wfit) {
    if (wafit) {
      std::cerr << "Error: omega_m, omega_de, w0, wa fit not supported"
		<< std::endl;
      return 4;
    }
    ft = utility::omegam_omegade_w;
  } else ft = utility::omegam_omegade;
  
  /*
  if (wafit) {
    //Right now they aren't working well
    std::cerr << "w_a fits not yet implemented" << std::endl;
    return 8;
  }
  */
  
  chifunc fFCN( sne, intrinsicdisp, pecz, ft, equiv_width_cut );
  if (fixalphaerrset) fFCN.SetAlphaErrFix( fixalphaerrval );
  if (fixbetaerrset) fFCN.SetBetaErrFix( fixbetaerrval );
  if (useKomatsuForm) {
    fFCN.SetUseKomatsuForm();
    fFCN.SetAtrans( atrans );
  }
  if (includebao) fFCN.SetIncludeBAO( true );
  if (includebao && useEisenstein) fFCN.SetUseEisenstein( true );
  if (includewmap) fFCN.SetIncludeWMAP( true );

  
  //Initialize the errors 
  try {
    fFCN.CalcPreErrs();
  } catch (const CosFitterExcept& ex) {
    std::cerr << "Error performing fit" << std::endl;
    std::cerr << ex << std::endl;
    std::cerr << "Aborting" << std::endl;
    return 2;
  }
  
  //Set up parameters
  double alpha_err=0.1;
  double beta_err=0.1;
  double om_err=0.1;
  double ode_err=0.1;
  double w0_err = 0.1;
  double wa_err = 0.5;
  double scriptm, scriptm_err = 0.1;
  
  MnUserParameters upar;
  upar.Add("\\Omega_m",om,om_err);
  if (fixom) upar.Fix("\\Omega_m"); else upar.SetLowerLimit("\\Omega_m",0.0);
  if (!flatfit) {
    upar.Add("\\Omega_DE",ode,ode_err);
    if (fixode) upar.Fix("\\Omega_DE");
  }
  if (wfit) {
    upar.Add("w_0",w0,w0_err);
    if (fixw0) upar.Fix("w_0");
  }
  if (wafit) {
    upar.Add("w_a",wa,wa_err);
    if (fixwa) upar.Fix("w_a");
  }
  upar.Add("\\alpha",alpha,alpha_err);
  upar.Add("\\beta",beta,beta_err);
  if (fixalpha) upar.Fix("\\alpha");
  if (fixbeta) upar.Fix("\\beta");
  scriptm = fFCN.EstimateScriptm( upar.Params() );
  upar.Add("\\scriptM_1",scriptm,scriptm_err);    
  upar.Add("\\scriptM_2",scriptm,scriptm_err);    
  
  //Create the minimizer
  MnMachinePrecision mprec;
  mprec.SetPrecision(1e-5);
  MnMigrad migrad( fFCN, upar);
  
  int nparams = upar.Parameters().size();
  
 //Figure free params
  std::vector< std::string > free_params;
  unsigned int nfree = 0;
  for (int i = 0; i < nparams; ++i)
    if (!upar.Parameter(i).IsFixed()) 
      free_params.push_back( upar.Name(i) );
  nfree = free_params.size();


  //And minimize
  std::cout << "Starting minimization" << std::endl;
  try {
    FunctionMinimum min = migrad();
    std::cout << "Chi2 " << std::left << std::fixed <<
      std::setprecision(2) << min.Fval() << " for " << 
      sne.size() - nfree << " dof" << std::endl;
    std::cout << "RMS: " << std::left << std::fixed 
	      << std::setprecision(3) 
	      << fFCN.GetRMS( min.UserState().Params() ) 
	      << " mag" << std::endl;
    if (verbose) {
      fFCN.print( std::cout, min.UserState().Params() );
      if ( (includebao && includewmap) && (!useEisenstein) ) {
	auxconstraint::wmap7_bao_P09 aux(ft);
	std::cout << "Percival09+WMAP7 chi2 contribution: " 
		  << std::left << std::fixed 
		  << std::setprecision(2) << aux( min.UserState().Params() )
		  << std::endl;
      } else {
	if (includebao) {
	  if (useEisenstein) {
	    auxconstraint::sdss_lrg_bao bao( ft );
	    std::cout << "E05 BAO chi2 contribution: " 
		      << std::left << std::fixed 
		      << std::setprecision(2) << bao( min.UserState().Params() )
		      << std::endl;
	    std::cout << "Value of A at best fit: " << 
	      bao.aval(min.UserState().Params() ) << std::endl;
	  } else {
	    auxconstraint::bao_P09 bao( ft );
	    std::cout << "P09 BAO chi2 contribution: " 
		      << std::left << std::fixed 
		      << std::setprecision(2) << bao( min.UserState().Params() )
		      << std::endl;
	  }
	  if (includewmap) {
	    auxconstraint::wmap7yr_dls wmap( ft );
	    std::cout << "WMAP7 chi2 contribution: " 
		      << std::left << std::fixed 
		      << std::setprecision(2)
		      << wmap( min.UserState().Params() )
		      << std::endl;
	  }
	}
      }
    } 
    
    if ( ! min.IsValid() ) {
      std::cout << "Error -- minimization was not successful.  Parameters"
		<< " probably not accurate" << std::endl;
      std::cout << "Function value: " << min.Fval() << std::endl;
      std::cout << "EDM: " << min.Edm() << std::endl;
      std::cout << "Number of calls: " << min.NFcn() << std::endl;
      if (min.HasMadePosDefCovar() ) { 
	std::cout << "      Covar was made pos def" 
		  << std::endl;
	return -11; 
      }
      if (min.HesseFailed() ) std::cout << "      Hesse is not valid" 
					<< std::endl;
      if (min.IsAboveMaxEdm() ) std::cout << "      Edm is above max" << 
	std::endl;
      if (min.HasReachedCallLimit() ) std::cout << 
	"      Reached call limit" << std::endl;

      //Print values anyways
      for (int i = 0; i < nparams; ++i) {
	if ( ! min.UserState().Parameter(i).IsFixed() ) {
	  std::cout << upar.Name(i) << " " << min.UserState().Value(i) 
		    << " " << min.UserState().Error(i) << std::endl;
	}
      }
      return -10;  //Not successful
    }    
    
    //Get errors and print them
    MnMinos minos( fFCN, min );
    std::pair<double,double> e;
    for (int i = 0; i < nparams; ++i) {
      if ( ! min.UserState().Parameter(i).IsFixed() ) {
	e = minos(i);
	std::cout << std::left << std::setw(10) << upar.Name(i) 
		  << " " << std::right << std::setw(7) << std::fixed 
		  << std::setprecision(3) << min.UserState().Value(i) << " +" 
		  << std::right << std::fixed << std::setprecision(3)
		  << e.second << " " << std::right << std::fixed 
		  << std::setprecision(3) << e.first << " (av err: " 
		  << std::left << std::fixed << std::setprecision(3) 
		  << 0.5*( fabs(e.first)+fabs(e.second) ) << ")" << std::endl;
      } else {
	std::cout << std::left << std::setw(10) << upar.Name(i)
		  << " " << std::right << std::setw(7) << std::fixed
		  << std::setprecision(3) << min.UserState().Value(i)
		  << " (fixed)" << std::endl;
      }
    }
    std::cout << "\\scriptM_1 is equivalent to M = " <<
      min.UserState().Value("\\scriptM_1") + 5 * log10(0.7) - 42.38410 <<
      " for H_0 = 70 km/sec/Mpc" << std::endl;
    std::cout << "\\scriptM_2 is equivalent to M = " <<
      min.UserState().Value("\\scriptM_2") + 5 * log10(0.7) - 42.38410 <<
      " for H_0 = 70 km/sec/Mpc" << std::endl;
    

    //Output covariances
    if (verbose) {
      MnUserCovariance cov = min.UserCovariance();
      std::vector<double> vars;
      std::stringstream covstring;
      vars.resize( cov.Nrow() );
      for (unsigned int i = 0; i < cov.Nrow(); ++i)
	vars[i] = cov(i,i);
      std::cout << "Covariance information: " << std::endl;
      for (unsigned int i=0; i < cov.Nrow(); ++i) {
	for (unsigned int j = i; j < cov.Nrow(); ++j) {
	  covstring.str(std::string(""));
	  covstring << "cov[" << free_params[i] << ", " 
		    << free_params[j] << "]";
	  std::cout << std::left << std::setw(26) << covstring.str()
		    << " = " << std::setw(8) << std::right << std::fixed 
		    << std::setprecision(5) << cov(i,j)
		    << "  rho = " << std::setw(8) << std::fixed 
		    << std::setprecision(5) 
		    << cov(i,j)/sqrt( vars[i] * vars[j] )
		    << std::endl;
	}
      }
    }

    //Draw contours if asked
    if (docontours) {
      if (ncontours < 1) {
	std::cerr << "Error: number of contour points invalid" << std::endl;
	return 4;
      }

      MnContours contours( fFCN, min );
      std::vector<std::pair<double,double> > conts = contours(0,1,ncontours);
      MnPlot plot;
      std::cout << "Contours for: " << upar.Name(0) << " vs " 
		<< upar.Name(1) << std::endl;
      plot(min.UserState().Value(static_cast<unsigned int>(0)), 
	   min.UserState().Value(static_cast<unsigned int>(1)),
	   conts);
      
    }
    
  } catch (const CosFitterExcept& ex) {
    std::cerr << "Error performing fit" << std::endl;
    std::cerr << " To datafile: " << datafile << std::endl;
    std::cerr << ex << std::endl;
    std::cerr << "Aborting" << std::endl;
    return 2;
  }

  return 0;
}

void chifunc::printResids(std::ostream& os,
			  const std::vector<double>& par) const {
  if (isprepped != true ) throw CosFitterExcept("printResids","printResids",
                                                "Not prepped",1);
  double w, wa, om, ode, alpha, beta, scriptm1, scriptm2;
  utility::parse_param_values( par, fittype, om, ode, w, wa, alpha, beta,
                               scriptm1, scriptm2 );

  int st = lmdist.getLumDist( sne, dl, om, ode, w, wa );
  if (st != 0) return; //Failure

  os << "Summary of SN data:" << std::endl;
  double predmag;
  for (unsigned int i = 0; i < nsn; ++i) {
    predmag = dl[i] - alpha * (sne[i].widthpar-1.0) +
      + beta * sne[i].colourpar;
    if (sne[i].equiv_widthpar <= equiv_widthpar_cut)
      predmag += scriptm1; else predmag += scriptm2;
    os << std::left << std::setw(8) << sne[i].name << " z: " << std::right <<
       std::setw(6) << std::fixed << std::setprecision(4) << sne[i].zhel << 
      " s: " << std::right << std::setw(5) << std::fixed << 
      std::setprecision(3) << sne[i].widthpar << " c: " << 
      std::fixed << std::setprecision(3) << std::right << 
      std::setw(6) << sne[i].colourpar;
    os << " eq: " << std::right << std::setw(5) << std::fixed << 
      std::setprecision(3) << sne[i].equiv_widthpar;
    os << " var: " << std::setw(6) << std::fixed <<
      std::right << pre_vars[i] << " sigma: " << std::fixed << 
      std::setw(5) << std::right << sqrt(pre_vars[i]);
    os << " predmag: " << std::setw(6) << std::fixed
       << std::right << predmag << " diffmag: " 
       << std::setw(6) << std::fixed << std::right
       << sne[i].mag - predmag;
    os << std::endl;
  }
}
