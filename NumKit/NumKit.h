//
//  NumKit.h
//
//  Created by Koichi Oshio on 2-26-2013.
//
//	=== plans (linear) ===
//	use BLAS / LAPACK
//		- complex type is interleaved
//		- use RecImage (row major, split) for dsp type tasks,
//			and convert to Num_mat (column major, interleaved)
//		- make classes
//	=== plans (Spin) ===
//		* implement scheffler version first (see how it works)
//		* oshio version (spherical PO)


#import <Foundation/Foundation.h>
#import <Accelerate/Accelerate.h>
#import <RecKit/RecKit.h>
#import <stdlib.h>
#import <math.h>

#define VECTOR(n)           (float *)calloc(sizeof(float), n)

// ======== Minimization =====
typedef struct {
    int     n;
    float   *x;
    float   *y;
    float   xscale;
    float   yscale;
} Num_data;

typedef struct {
    int     n;		// dimension
    int     *flag;	// active flag. currently not being used
    float   *data;	// param array
	float	*min;
	float	*max;
} Num_param;

typedef struct {
    int     n;
    float   *min;
    float   *inc;
    int     *steps;
    int     *curr;
} Num_range;

// ======== GMM =====
typedef struct {
    int     n;		// number of classes
    float   *mn;	// mean
    float   *sd;	// SD
    float   *ht;	// height
} Num_gmm_param;

// ======== Linear =====

// type
enum {
	NUM_REAL = 0,		// single precision float
	NUM_COMPLEX,		// single precision complex, interleaved
};

// order
enum {
	NUM_ROW_MAJOR = 0,
	NUM_COL_MAJOR,
};

// complex number is stored in interleaved format
// matrix is stored in column major order
typedef struct {
	int		type;	// real or complex
	int		order;	// make col major default (use RecImage for row major)
    int     nr;
    int     nc;
	int		ld;	// leading dim (nr if xxx, nc if xxx)
    float   *data;
} Num_mat;

typedef struct {
	int		type;
    int     n;
    float   *data;
} Num_vec;

// === Bloch simulation / phase graph ===
// == state (tmp output), real out is image
// === separate coherence (for phase graph) and vector (for image)
typedef struct {	// Woessner / Scheffler
	// order
	int		n;
	// state
	float	*Fx;	// 2n (real)
	float	*Fy;	// 2n (imag)
	float	*Zx;	// n  (real)
	float	*Zy;	// n  (imag !!!)
	// next state
	float	*Fxp;	// 2n (real)
	float	*Fyp;	// 2n (imag)
	float	*Zxp;	// n  (real)
	float	*Zyp;	// n  (imag !!!)
} Num_spin_s;

typedef struct {	// Oshio (spherical PO)
	int		n;	// order
	// state
	float	*Ipr;	// 2n (real)
	float	*Ipi;	// 2n (imag)
	float	*Imr;	// 2n (real)
	float	*Imi;	// 2n (imag)
	float	*I0r;	// 2n (real)
	float	*I0i;	// 2n (imag)
	// next state
	float	*nIpr;	// 2n (real)
	float	*nIpi;	// 2n (imag)
	float	*nImr;	// 2n (real)
	float	*nImi;	// 2n (imag)
	float	*nI0r;	// 2n (real)
	float	*nI0i;	// 2n (imag)
} Num_spin;

// == input for phase graph
typedef struct {
	int		n;
	float	*rho;
	float	*theta;
	float	tr;
	float	te;
	float	t1;
	float	t2;
} Num_rf;

// ======= Object version ====
@interface NumVector : NSObject
{
	int				type;
	int				length;
	NSMutableData	*data;
}

+ (NumVector *)vectorOfType:(int)type length:(int)n;
- (NumVector *)initWithType:(int)type length:(int)len;
- (int)length;
- (float *)data;
- (void)clear;
- (void)normal;
- (void)ramp;

@end

@interface NumMatrix : NSObject
{
	int				type;
	int				nrow;
	int				ncol;
	NSMutableData	*data;
}

+ (NumMatrix *)matrixOfType:(int)type nrow:(int)nr ncol:(int)nc;
+ (NumMatrix *)matrixWithNumMat:(Num_mat *)mat;
- (NumMatrix *)initWithType:(int)type nrow:(int)nr ncol:(int)nc;
- (float *)data;
- (int)type;
- (int)nrow;
- (int)ncol;
- (void)clear;
- (void)normal;
- (NumMatrix *)copy;
- (NumMatrix *)trans;
- (NumMatrix *)unitMatrixOfDim:(int)n;
- (NumMatrix *)multWithVector:(NumVector *)v;
- (NumMatrix *)multWithMatrix:(NumMatrix *)m;
- (NumMatrix *)multWithMatrix:(NumMatrix *)m transSelf:(BOOL)ts transMat:(BOOL)tm;
- (NumMatrix *)diagMatrixWithVector:(NumVector *)v;

- (NSDictionary *)svd;	// returns dict with entry: "U", "s", "Vt"
- (NSDictionary *)ica;	// returns dict with entry: "WX", "W"

@end
// ======= Object version ====


// === special functions ===
float	betai(float a, float b, float x);	// incomplete beta function
float	beta(float a, float b);				// beta function

// === common ===
void        Num_error(char *msg);   // error exit

// === Minimization / Least squares ===
        // non-linear least squares (uses conj gr)
int         Num_least_sq_old(Num_data *data, Num_param *param, float (*model)(float x, float *p, float *dy), float *mse);
int         Num_least_sq(Num_data *data, Num_param *param, float (^model)(float x, float *p, float *dy), float *mse);
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
void        Num_normalize_y(Num_data *d);

// === private / low level ====
int         Num_search_min(Num_data *d, Num_param *p, Num_range *r, float (^model)(float x, float *p, float *dy));
float       Num_mrqmin(Num_data *d, Num_param *p, float *chisq,
                    float (*model)(float x, float *p, float *dy), float *lmd);
float       Num_mrqcof(Num_data *d, Num_param *p, float *alpha, float *beta,
                    float (*model)(float x, float *pp, float *ddy));
void        Num_mnbrak(float *ax, float *bx, float *cx, float (^f1d)(float));
float       Num_brent(float ax, float bx, float cx,
                    float (^f1d)(float), float tol, float *xmin);   // Brent minimization
void        Num_linmin(float *p, float *pmn, float *pmx, float *gr, int n, float *fret, float (^f1d)(float *));
float       Num_chk_range(float *p, float *pmin, float *pmax, float *xi, int np);
void        Num_init_brak(float *ax, float *bx, float *cx, float max_delta, float (^f1d)(float));
float		Num_max_delta(float *p, float *pmn, float *pmx, float *xi, int np);

// === Statistics ===
float       Num_nrml(float m, float sd);
float       Num_unif(float m, float sd);
float       Num_mean(float *p, int n);
float       Num_var(float *p, int n, float m);
void        Num_avevar(float *p, int n, float *ave, float *var);
void        Num_covar(float *x, float *y, int n, float *av, float *cov);
void        Num_linreg(float *r, float *b0, float *b1, float *x, float *y, int skip, int n);
void        Num_orth(float *pr1, float *pr2, float *re1, float *im1, float *re2, float *im2, int n, float *th);
void        Num_ttest(float *p1, int n1, float *p2, int n2, float *tval, float *pval);
void        Num_ftest(float *p1, int n1, float *p2, int n2, float *fval, float *pval);
float		Num_tdist(float tval, float df);
float		Num_fdist(float fval, float df1, float df2);
float		Num_gdist(float m, float sd);
float       Num_t2p(float tval, float df);
float       Num_tlimit(float p, float df);
float       Num_f2p(float fval, float df1, float df2);
float       Num_flimit(float p, float df1, float df2, BOOL left);
float       Num_r2(float *p, float *f, int n);    // R2
float       Num_zbrent(float ax, float bx, float (^f1d)(float));    // Brent root finding
// === Complex Statistics ===
void        Num_cnrml(float *re, float *im, float mr, float mi, float v);
void		Num_cmean(float *re, float *im, int n, float *mr, float *mi);
void        Num_cavevar(float *re, float *im, int n, float *avr, float *avi, float *var);

void		Num_ctmean(float *p1, float *q1, int n1, float *p2, float *q2, int n2, float *t2val, float *pval, float *th);
void        Num_ctvar(float *p1, float *q1, int n1, float *p2, float *q2, int n2, float *tval, float *pval, float *th);

// === Gaussian Mixture Model ===
Num_gmm_param *Num_alloc_gmm(int n);
void		Num_free_gmm(Num_gmm_param * g);
float		Num_nrml_dist(float x, float m, float sd);
int			Num_gmm(Num_data *d, Num_gmm_param *g);	// *g is initial / result, return val: # iteration, -1 for error return

// === Linear Algebra ===
// conversion to/from RecImage
// === real / complex
Num_mat	*	Num_im_to_m(RecImage *im);
Num_vec *	Num_im_to_v(RecImage *im);
RecImage *	Num_m_to_im(Num_mat *m);
RecImage *	Num_v_to_im(Num_vec *v);
RecImage *	Num_alloc_im_with_m(Num_mat *m);
RecImage *	Num_alloc_im_with_v(Num_vec *v);
Num_mat *	Num_alloc_m_with_im(RecImage *im);
Num_vec *	Num_alloc_v_with_im(RecImage *im);
void		Num_copy_im_to_m(Num_mat *m, RecImage *im);
void		Num_copy_im_to_v(Num_vec *v, RecImage *im);
void		Num_copy_m_to_im(RecImage *im, Num_mat *m);
void		Num_copy_v_to_im(RecImage *im, Num_vec *v);

Num_mat *	Num_new_mat(int nr, int nc);
Num_vec *	Num_new_vec(int n);
Num_mat *	Num_new_cmat(int nr, int nc);
Num_vec *	Num_new_cvec(int n);
void		Num_make_cmat(Num_mat *m);
void		Num_make_cvec(Num_vec *v);
void		Num_copy_vec(Num_vec *b, Num_vec *a);			// b <- a
void        Num_free_mat(Num_mat *m);
void        Num_free_vec(Num_vec *v);
void        Num_clear_mat(Num_mat *a);						// A <- 0
void		Num_tocm(Num_mat *m);	// convert to column major
void		Num_torm(Num_mat *m);	// convert to row major

// vvvvvvvvvvvvvvvvvvvvvvvvvvv
// ### implement complex version
float       Num_dotpr(Num_vec *v1, Num_vec *v2);			// a = v1 * v2
void        Num_mvmul(Num_vec *b, Num_mat *a, Num_vec *x);	// b = Ax
void        Num_mmul(Num_mat *c, Num_mat *a, Num_mat *b);	// C = AB
void        Num_mtmul(Num_mat *c, Num_mat *a, BOOL ta, Num_mat *b, BOOL tb);	// C = A(t)B(t)
void		Num_col_vec(Num_vec *x, Num_mat *A, int ix);	// copy i-th col-v of A
void		Num_row_vec(Num_vec *x, Num_mat *A, int ix);	// copy i-th row-v of A
void		Num_mmul_const(Num_mat *A, float c);
// rewtite using BLAS
void        Num_copy_mat(Num_mat *b, Num_mat *a);						// B <- A
void		Num_copy_sub_mat(Num_mat *b, Num_mat *a, int rc0, int rn, int c0, int cn);
// replace following with LAPACK version
Num_mat *	Num_trans(Num_mat *a);							// B = At, A:mxn, B:nxm
void		Num_trans_ip(Num_mat *b, Num_mat *a);			// B = At
void        Num_hermit(Num_mat *b, Num_mat *a);							// B = At, A:mxn, B:nxm
void        Num_unit_mat(Num_mat *a);									// A <- I
// depricate and write new routines using LAPACK
int			Num_gaussj(Num_mat *a, Num_mat *b);							// A, b -> A-1, x, b can be NULL
//void        Num_jacobi(Num_mat *m, Num_mat *evec, Num_vec *eval);
void        Num_pinv(Num_mat *b, Num_mat *a);							// B = A+, A:mxn, B:nxm
void        Num_vadd(Num_vec *v1, Num_vec *v2, Num_vec *res);			// res = v1 + v2
void        Num_vsma(Num_vec *x, float a, Num_vec *y);					// y = a * x + y
void		Num_grmsch(Num_mat *a);										// make ortho-nomal

// complex not done yet ###
int			Num_inv(Num_mat *A, Num_mat *B);
void		Num_evd_ref(Num_mat *A, Num_mat *Evec, Num_vec *eval);
void		Num_normalize_vec(Num_vec *v);

// low level
void		Num_orth_mat(Num_mat *A);
void		Num_diag_ip(Num_mat *d, Num_vec *v);
Num_mat	*	Num_diag(Num_vec *v);
Num_mat *	Num_sort_mat(Num_vec *v);
void		Num_scale_columns(Num_mat *A, Num_vec *v);
void		Num_scale_rows(Num_mat *A, Num_vec *v);
void		Num_negate_mat(Num_mat *A);						// for SVD result

// returns struct
typedef struct {
	Num_mat	*U;
	Num_vec	*s;
	Num_mat	*Vt;
} Num_svd_result;

typedef struct {
	Num_mat	*W;		// unmixing matrix
	Num_mat	*WK;	// map
	Num_mat	*WX;	// time course
} Num_ica_result;

typedef struct {
	Num_mat	*Evec;
	Num_vec	*eval;
} Num_evd_result;

void				Num_col_center(Num_mat *A);
void				Num_row_center(Num_mat *A);
Num_evd_result *	Num_evd(Num_mat *A);
Num_svd_result *	Num_svd(Num_mat *A);
Num_svd_result *	Num_new_svd_result(Num_mat *A);
void				Num_svd_ref(Num_mat *A, Num_mat *U, Num_vec *s, Num_mat *Vt);
Num_ica_result *	Num_ica(Num_mat *U, int nvec);
void				Num_free_svd_result(Num_svd_result *r);
void				Num_free_ica_result(Num_ica_result *r);
void				Num_free_evd_result(Num_evd_result *r);

// === ODE ===
void        Num_rk4(float *x, float *y, float (^deriv)(float x, float y), float step);    // not tested yet

// === Bloch simulation / phase graph ===
Num_spin *	Num_new_spin(int order);
void		Num_free_spin(Num_spin *sp);
Num_rf *	Num_new_rf(int order);
void		Num_free_rf(Num_rf *sp);
RecImage *	Num_phase_graph(Num_rf *rf, BOOL cpmg);
RecImage *	Num_phase_graph_s(Num_rf *rf, BOOL cpmg);	// Scheffler

//=== dbg
void        dump_vec(Num_vec *v);
void        dump_mat(Num_mat *m);
void		saveAsKOImage(Num_mat *m, NSString *path);
int			diff_mat(Num_mat *a, Num_mat *b);

