/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

// LIEGROUP_API double* GravityForces_wrapper(const double* thetalist_array) ;
CASADI_SYMBOL_EXPORT void pRNE(const double* thetalist_in, const double* dthetalist_in, const double* ddthetalist_in,
          const double* thetalist_ref_in, const double* dthetalist_ref_in, const double* ddthetalist_ref_in,
          const double* g_in, double* taulist_out) ;
#ifdef __cplusplus
}
#endif