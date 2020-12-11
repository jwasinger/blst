/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_FIELDS_H__
#define __BLS12_381_ASM_FIELDS_H__

#include "vect.h"
#include "consts.h"

/*
 * BLS12-381-specifc Fp shortcuts to assembly.
 */
static inline void add_fp(vec384 ret, const vec384 a, const vec384 b)
{   add_mod_384(ret, a, b, BLS12_381_P);   }

static inline void sub_fp(vec384 ret, const vec384 a, const vec384 b)
{   sub_mod_384(ret, a, b, BLS12_381_P);   }

static inline void mul_by_3_fp(vec384 ret, const vec384 a)
{   mul_by_3_mod_384(ret, a, BLS12_381_P);   }

static inline void mul_by_8_fp(vec384 ret, const vec384 a)
{   mul_by_8_mod_384(ret, a, BLS12_381_P);   }

static inline void lshift_fp(vec384 ret, const vec384 a, size_t count)
{   lshift_mod_384(ret, a, count, BLS12_381_P);   }

static inline void mul_fp(vec384 ret, const vec384 a, const vec384 b)
{   mul_mont_384(ret, a, b, BLS12_381_P, p0);   }

static inline void sqr_fp(vec384 ret, const vec384 a)
{   sqr_mont_384(ret, a, BLS12_381_P, p0);   }

static inline void cneg_fp(vec384 ret, const vec384 a, limb_t flag)
{   cneg_mod_384(ret, a, flag, BLS12_381_P);   }

#define neg_fp(r,a) cneg_fp((r),(a),1)

static inline void eucl_inverse_fp(vec384 ret, const vec384 a)
{   eucl_inverse_mod_384(ret, a, BLS12_381_P, BLS12_381_RR);   }

static inline void from_fp(vec384 ret, const vec384 a)
{   from_mont_384(ret, a, BLS12_381_P, p0);   }

/*
 * BLS12-381-specifc Fp2 shortcuts to assembly.
 */
static inline void add_fp2(vec384x ret, const vec384x a, const vec384x b)
{   add_mod_384x(ret, a, b, BLS12_381_P);   }

static inline void sub_fp2(vec384x ret, const vec384x a, const vec384x b)
{   sub_mod_384x(ret, a, b, BLS12_381_P);   }

static inline void mul_by_3_fp2(vec384x ret, const vec384x a)
{   mul_by_3_mod_384x(ret, a, BLS12_381_P);   }

static inline void mul_by_8_fp2(vec384x ret, const vec384x a)
{   mul_by_8_mod_384x(ret, a, BLS12_381_P);   }

static inline void lshift_fp2(vec384x ret, const vec384x a, size_t count)
{
    lshift_mod_384(ret[0], a[0], count, BLS12_381_P);
    lshift_mod_384(ret[1], a[1], count, BLS12_381_P);
}

static inline void mul_fp2(vec384x ret, const vec384x a, const vec384x b)
{   

// paul: flip this switch to 
#define MUL_FP2_SWITCH 2

#if MUL_FP2_SWITCH==0
  // optimized version
  mul_mont_384x(ret, a, b, BLS12_381_P, p0);   
#endif

#if MUL_FP2_SWITCH==1
  // gives same result as above mul_mont_384x()
  vec768 t0, t1, t2;
  vec384 aa, bb;
  mul_384(t0, a[0], b[0]);
  mul_384(t1, a[1], b[1]);
  add_mod_384(aa, a[0], a[1], BLS12_381_P);
  add_mod_384(bb, b[0], b[1], BLS12_381_P);
  mul_384(t2, aa, bb);
  sub_mod_384x384(t2, t2, t0, BLS12_381_P);
  sub_mod_384x384(t2, t2, t1, BLS12_381_P);
  sub_mod_384x384(t0, t0, t1, BLS12_381_P);
  redc_mont_384(ret[0], t0, BLS12_381_P, p0);
  redc_mont_384(ret[1], t2, BLS12_381_P, p0);
#endif

#if MUL_FP2_SWITCH==2
  // the unoptimized mul_fp2 shared by Kelly Olson
  vec384 aa, bb, cc;

  add_mod_384(aa, a[0], a[1], BLS12_381_P);
  add_mod_384(bb, b[0], b[1], BLS12_381_P);
  mul_mont_384(bb, bb, aa, BLS12_381_P, p0);
  mul_mont_384(aa, a[0], b[0], BLS12_381_P, p0);
  mul_mont_384(cc, a[1], b[1], BLS12_381_P, p0);
  sub_mod_384(ret[0], aa, cc, BLS12_381_P);
  sub_mod_384(ret[1], bb, aa, BLS12_381_P);
  sub_mod_384(ret[1], ret[1], cc, BLS12_381_P);
#endif

#if MUL_FP2_SWITCH==3
  // naive mul_fp2 using four muls
  vec384 t0, t1;
  mul_mont_384(t0, a[0], b[0], BLS12_381_P, p0);
  mul_mont_384(t1, a[1], b[1], BLS12_381_P, p0);
  sub_mod_384(ret[0], t0, t1, BLS12_381_P);
  mul_mont_384(t0, a[0], b[1], BLS12_381_P, p0);
  mul_mont_384(t1, a[1], b[0], BLS12_381_P, p0);
  add_mod_384(ret[1], t0, t1, BLS12_381_P);
#endif

}

static inline void sqr_fp2(vec384x ret, const vec384x a)
{   sqr_mont_384x(ret, a, BLS12_381_P, p0);   }

static inline void cneg_fp2(vec384x ret, const vec384x a, limb_t flag)
{   
    cneg_mod_384(ret[0], a[0], flag, BLS12_381_P);
    cneg_mod_384(ret[1], a[1], flag, BLS12_381_P);
}

#define neg_fp2(r,a) cneg_fp2((r),(a),1)

typedef vec384x   vec384fp2;
typedef vec384fp2 vec384fp6[3];
typedef vec384fp6 vec384fp12[2];

static void sqr_fp12(vec384fp12 ret, const vec384fp12 a);
static void cyclotomic_sqr_fp12(vec384fp12 ret, const vec384fp12 a);
static void mul_fp12(vec384fp12 ret, const vec384fp12 a, const vec384fp12 b);
static void mul_by_xy00z0_fp12(vec384fp12 ret, const vec384fp12 a,
                                               const vec384fp6 xy00z0);
static void conjugate_fp12(vec384fp12 a);
static void inverse_fp12(vec384fp12 ret, const vec384fp12 a);
/* caveat lector! |n| has to be non-zero and not more than 3! */
static void frobenius_map_fp12(vec384fp12 ret, const vec384fp12 a, size_t n);

#endif /* __BLS12_381_ASM_FIELDS_H__ */
