#include "include/complexCL.h"

// Source: http://stackoverflow.com/a/10042415

/*
 * Return Real (Imaginary) component of complex number:
 */
float creal(cfloat a){
     return a.x;
}
float cimag(cfloat a){
     return a.y;
}

/*
 * Get the modulus of a complex number (its length):
 */
float cmod(cfloat a){
    return (sqrt(a.x*a.x + a.y*a.y));
}

/*
 * Get the argument of a complex number (its angle):
 * http://en.wikipedia.org/wiki/Complex_number#Absolute_value_and_argument
 */
float carg(cfloat a){
    if(a.x > 0){
        return atan(a.y / a.x);

    }else if(a.x < 0 && a.y >= 0){
        return atan(a.y / a.x) + M_PI;

    }else if(a.x < 0 && a.y < 0){
        return atan(a.y / a.x) - M_PI;

    }else if(a.x == 0 && a.y > 0){
        return M_PI/2;

    }else if(a.x == 0 && a.y < 0){
        return -M_PI/2;

    }else{
        return 0;
    }
}

cfloat cneg(cfloat cx){
#ifdef HOSTBUILD
    return (cfloat){-cx.x, -cx.y};
#else
    return (cfloat)(-cx.x, -cx.y);
#endif
}

// Complex conjugate a number
cfloat conj(cfloat cx){
#ifdef HOSTBUILD
    return (cfloat){cx.x, -cx.y};
#else
    return (cfloat)(cx.x, -cx.y);
#endif
}

cfloat cadd(cfloat a, cfloat b){
#ifdef HOSTBUILD
    return (cfloat){a.x + b.x, a.y + b.y};
#else
    return (cfloat)(a.x + b.x, a.y + b.y);
#endif
}

cfloat cadd(cfloat a, float b){
#ifdef HOSTBUILD
    return (cfloat){a.x + b, a.y};
#else
    return (cfloat)(a.x + b, a.y);
#endif
}

cfloat csub(cfloat a, cfloat b){
#ifdef HOSTBUILD
    return (cfloat){a.x - b.x, a.y - b.y};
#else
    return (cfloat)(a.x - b.x, a.y - b.y);
#endif
}

/*
 * Multiply two complex numbers:
 *
 *  a = (aReal + I*aImag)
 *  b = (bReal + I*bImag)
 *  a * b = (aReal + I*aImag) * (bReal + I*bImag)
 *        = aReal*bReal +I*aReal*bImag +I*aImag*bReal +I^2*aImag*bImag
 *        = (aReal*bReal - aImag*bImag) + I*(aReal*bImag + aImag*bReal)
 */
cfloat cmult(cfloat a, cfloat b){
#ifdef HOSTBUILD
    return (cfloat){ a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x};
#else
    return (cfloat)( a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
#endif
}

cfloat cmult(cfloat a, float b){
#ifdef HOSTBUILD
    return (cfloat){ a.x*b, a.y*b};
#else
    return (cfloat)( a.x*b, a.y*b);
#endif
}
/*
 * Divide two complex numbers:
 *
 *  aReal + I*aImag     (aReal + I*aImag) * (bReal - I*bImag)
 * ----------------- = ---------------------------------------
 *  bReal + I*bImag     (bReal + I*bImag) * (bReal - I*bImag)
 * 
 *        aReal*bReal - I*aReal*bImag + I*aImag*bReal - I^2*aImag*bImag
 *     = ---------------------------------------------------------------
 *            bReal^2 - I*bReal*bImag + I*bImag*bReal  -I^2*bImag^2
 * 
 *        aReal*bReal + aImag*bImag         aImag*bReal - Real*bImag 
 *     = ---------------------------- + I* --------------------------
 *            bReal^2 + bImag^2                bReal^2 + bImag^2
 * 
 */
cfloat cdiv(cfloat a, cfloat b){
#ifdef HOSTBUILD
    return (cfloat){(a.x*b.x + a.y*b.y)/(b.x*b.x + b.y*b.y), (a.y*b.x - a.x*b.y)/(b.x*b.x + b.y*b.y)};
#else
    return (cfloat)((a.x*b.x + a.y*b.y)/(b.x*b.x + b.y*b.y), (a.y*b.x - a.x*b.y)/(b.x*b.x + b.y*b.y));
#endif
}

cfloat cdiv(cfloat a, float b){
#ifdef HOSTBUILD
    return (cfloat){a.x/b, a.y/b};
#else
    return (cfloat)(a.x/b, a.y/b);
#endif
}

/*
 *  Square root of complex number.
 *  Although a complex number has two square roots, numerically we will
 *  only determine one of them -the principal square root, see wikipedia
 *  for more info: 
 *  http://en.wikipedia.org/wiki/Square_root#Principal_square_root_of_a_complex_number
 */
// cfloat csqrt(cfloat a){
//     return (cfloat)( sqrt(cmod(a)) * cos(carg(a)/2),  sqrt(cmod(a)) * sin(carg(a)/2));
// }
