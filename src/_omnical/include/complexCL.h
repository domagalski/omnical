#ifdef HOSTBUILD

#include <math.h>
#ifdef __APPLE__
    #include <OpenCL/opencl.h>
#else
    #include <CL/cl.h>
#endif
#ifdef AMD
    #include <CL/cl_ext.h>
#endif

typedef cl_float2 cfloat;
#else
typedef float2 cfloat;
#define I ((cfloat)(0.0, 1.0))
#endif

// Complex arithmetic functions
float creal(cfloat a);
float cimag(cfloat a);
float cmod(cfloat a);
float carg(cfloat a);
cfloat cneg(cfloat cx);
cfloat conj(cfloat cx);
cfloat cadd(cfloat a, cfloat b);
cfloat cadd(cfloat a, float b);
cfloat csub(cfloat a, cfloat b);
cfloat cmult(cfloat a, cfloat b);
cfloat cmult(cfloat a, float b);
cfloat cdiv(cfloat a, cfloat b);
cfloat cdiv(cfloat a, float b);
cfloat csqrt(cfloat a);
