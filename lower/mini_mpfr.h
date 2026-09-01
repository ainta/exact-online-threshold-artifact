#ifndef MINI_MPFR_H
#define MINI_MPFR_H
#include <gmp.h>
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef long mpfr_prec_t;
typedef long mpfr_exp_t;
typedef int mpfr_sign_t;
typedef struct {
  mpfr_prec_t _mpfr_prec;
  mpfr_sign_t _mpfr_sign;
  mpfr_exp_t _mpfr_exp;
  mp_limb_t *_mpfr_d;
} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];
typedef __mpfr_struct *mpfr_ptr;
typedef const __mpfr_struct *mpfr_srcptr;
typedef enum { MPFR_RNDN=0, MPFR_RNDZ=1, MPFR_RNDU=2, MPFR_RNDD=3, MPFR_RNDA=4, MPFR_RNDF=5, MPFR_RNDNA=-1 } mpfr_rnd_t;

void mpfr_init2(mpfr_ptr, mpfr_prec_t);
void mpfr_clear(mpfr_ptr);
int mpfr_set(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_set_ui(mpfr_ptr, unsigned long, mpfr_rnd_t);
int mpfr_set_si(mpfr_ptr, long, mpfr_rnd_t);
int mpfr_set_z(mpfr_ptr, mpz_srcptr, mpfr_rnd_t);
int mpfr_set_q(mpfr_ptr, mpq_srcptr, mpfr_rnd_t);
int mpfr_set_str(mpfr_ptr, const char*, int, mpfr_rnd_t);
int mpfr_set_zero(mpfr_ptr, int);
int mpfr_add(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sub(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_mul(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_div(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_add_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_sub_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_mul_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_div_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_ui_div(mpfr_ptr, unsigned long, mpfr_srcptr, mpfr_rnd_t);
int mpfr_neg(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sqrt(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_exp(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_log(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_erfc(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_const_pi(mpfr_ptr, mpfr_rnd_t);
int mpfr_pow(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr, mpfr_srcptr);
int mpfr_cmp_ui(mpfr_srcptr, unsigned long);
unsigned long mpfr_get_ui(mpfr_srcptr, mpfr_rnd_t);
double mpfr_get_d(mpfr_srcptr, mpfr_rnd_t);
int mpfr_get_z(mpz_ptr, mpfr_srcptr, mpfr_rnd_t);
char *mpfr_get_str(char*, mpfr_exp_t*, int, size_t, mpfr_srcptr, mpfr_rnd_t);
void mpfr_free_str(char*);
const char *mpfr_get_version(void);
#ifdef __cplusplus
}
#endif
#endif
