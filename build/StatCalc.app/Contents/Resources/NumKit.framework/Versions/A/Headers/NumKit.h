//
//  NumKit.h
//
//  Created by Koichi Oshio on 2-26-2013.
//

#import <Foundation/Foundation.h>
#import <Accelerate/Accelerate.h>
#import <stdlib.h>
#import <math.h>

#define VECTOR(n)           (float *)calloc(sizeof(float), n)

// => Num_data
typedef struct {
    int     nparam;
    float   *p;
    float   *dydp;
    // add param switch later
    int     ndata;
    float   *x;     // x
    float   *y;     // measured data
    float   *f;     // model (filled after fitting)
    float   xscale;
    float   yscale;
} LSconfig;

typedef struct {
    int     n;
    float   *x;
    float   *y;
    float   xscale;
    float   yscale;
} Num_data;

typedef struct {
    int     n;
    int     *flag;
    float   *data;
} Num_param;

typedef struct {
    int     n;
    float   *min;
    float   *inc;
    int     *steps;
    int     *curr;
} Num_range;

typedef struct {
    int     nr;
    int     nc;
    float   *data;
} Num_mat;

typedef struct {
    int     n;
    float   *data;
} Num_vec;

typedef struct {
    int     nr;
    int     nc;
    float   *re;
    float   *im;
} Num_cmat;

typedef struct {
    int     n;
    float   *re;
    float   *im;
} Num_cvec;

void        Num_error(char *msg);   // error exit

// === Minimization / Least squares ===
        // non-linear least squares (uses conj gr)
int         Num_least_sq(Num_data *data, Num_param *param, float (*model)(float x, float *p, float *dy), float *mse);
int         Num_least_sq_new(Num_data *data, Num_param *param, float (^model)(float x, float *p, float *dy), float *mse);
        // conjugate gradient
int         Num_conjgr(Num_param *param, float (^cost)(float *p), void (^dcost)(float *p, float *dy), float *minval);
        // conjugate direction
int         Num_powell(Num_param *param, float(^cost)(float *p), float *minval);
        // Marquardt-Levenberg
int         Num_marquardt(Num_data *d, Num_param *p, float (*model)(float x, float *pp, float *ddy), float *mse);
        // Data structure
Num_data  * Num_alloc_data(int ndata);
Num_param * Num_alloc_param(int nparam);
Num_range * Num_alloc_range(int nparam);
void        Num_free_data(Num_data *d);
void        Num_free_param(Num_param *p);
void        Num_free_range(Num_range *r);
void        Num_normalize_data(Num_data *d);

// === private / low level ====
int         Num_search_min(Num_data *d, Num_param *p, Num_range *r, float (*model)(float x, float *p, float *dy));
float       Num_mrqmin(Num_data *d, Num_param *p, float *chisq,
                    float (*model)(float x, float *p, float *dy), float *lmd);
float       Num_mrqcof(Num_data *d, Num_param *p, float *alpha, float *beta,
                    float (*model)(float x, float *pp, float *ddy));
void        Num_mnbrak(float *ax, float *bx, float *cx, float (^f1d)(float));
float       Num_brent(float ax, float bx, float cx,
                    float (^f1d)(float), float tol, float *xmin);   // Brent minimization
void        Num_linmin(float *p, float *gr, int n, float *fret, float (^f1d)(float *));
float       Num_chk_range(float *p, float *pmin, float *pmax, float *xi, int np);
void        Num_init_brak(float *ax, float *bx, float *cx, float *gr, float (^f1d)(float));

// === Statistics ===
float       Num_nrml(float m, float v);
float       Num_mean(float *p, int n);
void        Num_avevar(float *p, int n, float *ave, float *var);
float       Num_var(float *p, int skip, int n);
void        Num_covar(Num_vec *x, Num_vec *y, Num_vec *av, Num_mat *cov);
void        Num_linreg(float *r, float *b0, float *b1, float *x, float *y, int skip, int n);
void        Num_orth(float *pr1, float *pr2, float *re1, float *im1, float *re2, float *im2, int n, float *th);
void        Num_ttest(float *p1, int skip1, int n1, float *p2, int skip2, int n2, float *tval, float *pval);
void        Num_ftest(float *p1, int skip1, int n1, float *p2, int skip2, int n2, float *fval, float *pval);
float       Num_t2p(float tval, int df);
float       Num_tlimit(float p, int df);
float       Num_f2p(float fval, int df1, int df2);
float       Num_flimit(float p, int df1, int df2, BOOL left);
float       Num_r2(float *p, float *f, int skip, int n);    // R2
float       Num_zbrent(float ax, float bx, float (^f1d)(float));    // Brent root finding
// === Complex Statistics ===
void        Num_cnrml(float *re, float *im, float mr, float mi, float v);
void        Num_cavevar(Num_cvec *p, float *avr, float *avi, float *var);
void        Num_ctmean(Num_cvec *p1, Num_cvec *p2, float *t2val, float *pval, float *th);
void        Num_ctvar(Num_cvec *p1, Num_cvec *p2, float *tval, float *pval);

// === Linear Algebra ===
// real only first
Num_mat  *  Num_new_mat(int nr, int nc);
Num_vec  *  Num_new_vec(int n);
void        Num_free_mat(Num_mat *m);
void        Num_free_vec(Num_vec *v);
void        Num_change_mat(Num_mat *m, int nr, int nc);
void        Num_change_vec(Num_vec *v, int n);
void        Num_vadd(Num_vec *v1, Num_vec *v2, Num_vec *res);             // res = v1 + v2
void        Num_vsma(Num_vec *res, Num_vec *v1, float a, Num_vec *v2);    // v1 = v1 * a + v2
float       Num_dotpr(Num_vec *v1, Num_vec *v2);            // a = v1 * v2
void        Num_mvmul(Num_vec *b, Num_mat *a, Num_vec *x); // b = Ax
void        Num_mmul(Num_mat *c, Num_mat *a, Num_mat *b);   // C = AB
void        Num_gaussj(Num_mat *a, Num_mat *b);         // A, b -> A-1, x, b can be NULL
void        Num_jacobi(Num_mat *m, Num_mat *evec, Num_vec *eval);
void        Num_trans(Num_mat *b, Num_mat *a);            // B = At, A:mxn, B:nxm
void        Num_copy_mat(Num_mat *b, Num_mat *a);         // B <- A
void        Num_clear_mat(Num_mat *a);                  // A <- 0
void        Num_unit_mat(Num_mat *a);                        // A <- I
void        Num_pinv(Num_mat *b, Num_mat *a);             // B = A+, A:mxn, B:nxm
// === complex
Num_cmat *  Num_new_cmat(int nr, int nc);
Num_cvec *  Num_new_cvec(int n);
void        Num_free_cmat(Num_cmat *m);
void        Num_free_cvec(Num_cvec *v);


// === ODE ===
void        Num_rk4(float *x, float *y, float (*deriv)(float x, float y), float step);    // not tested yet
                                                                                        // tune for speed (expand)
// === etc ===


//=== dbg
//void        dump_lsc(LSconfig *lsc);
void        dump_vec(Num_vec *v);
void        dump_mat(Num_mat *m);

