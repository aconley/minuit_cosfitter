minuit_cosfitter
================

This is a supernova-centric fitter for measuring
the cosmological parameters, specifically for fitting the
dark energy equation of state.  This code minimizes the chi-square
using a variable-metric minimizer provided by
[MINUIT](http://seal.web.cern.ch/seal/MathLibs/Minuit2/html/).



###Installation
Installation is via the standard UNIX `configure` and
`make`. pofd_affine depends on the following packages:
* [MINUIT](http://seal.web.cern.ch/seal/MathLibs/Minuit2/html/)
* The [GNU scientific library](http://www.gnu.org/software/gsl/)
* It can benefit significantly from LAPACK/BLAS libraries,
   either [ATLAS](http://math-atlas.sourceforge.net/) or
   the [Intel MKL](http://software.intel.com/en-us/articles/intel-mkl/);
   on OSX there is also support for the built in accelerate libraries.
It may be necessary to tell configure where to look for these
libraries -- see `configure --help`.

### Documentation

The executable is called cosfitter -- some help can be obtained
via

	cosfitter --help

### Branches

The most interesting branch -- which will never be merged --
is a branch to handle the galaxy dependent offset in the absolute
magnitude of Type Ia SNe as described in
[Sullivan et al. (2010)](http://adsabs.harvard.edu/abs/2010MNRAS.406..782S).
This is the twoscriptm branch, and was the one actually used in the SNLS
3rd year analyses.  In this branch the executable is called
`cosfitter_twoscriptm`.

### References
* It was used in the SNLS 3rd year analysis in both
  [Conley et al. 2011](http://adsabs.harvard.edu/abs/2011ApJS..192....1C)
  and [Sullivan et al. 2011](http://adsabs.harvard.edu/abs/2011ApJ...737..102S)
