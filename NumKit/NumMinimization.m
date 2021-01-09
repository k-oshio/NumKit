//
// NumMinimization.m
//
//  == plans ==
//	implement param limits (started on 5-26-2016)
//	testing git ...
//

#import "NumKit.h"

#define GOLD                1.618034
#define CGOLD               0.3819660
#define GLIMIT              100.0
#define EPS                 1.0e-10 // -10
#define BRTOL               2.0e-4  //  2.0e-4
#define TINY                EPS     // -20
#define SHFT(a, b, c, d)    (a) = (b); (b) = (c); (c) = (d)
#define SIGN(a, b)          ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define ITMAX               800     //800	// orig:100

int     NumMinimization_dbg = 0;

// private
void
Num_rewind_range(Num_range *r)
{
    int     i;
    for (i = 0; i < r->n; i++) {
        r->curr[i] = 0;
    }
}

BOOL
Num_increment_range(Num_range *r)
{
    int     i, n = r->n;
    BOOL    in_range = NO;
    for (i = n - 1; i >= 0; i--) {
        r->curr[i] += 1;
        if (r->curr[i] < r->steps[i]) {
            in_range = YES;
            break;
        }
        // else
        r->curr[i] = 0;
    }
    return in_range;
}

void
Num_current_param(Num_param *p, Num_range *r)
{
    int     i, n = p->n;
    for (i = 0; i < n; i++) {
        p->data[i] = r->min[i] + r->inc[i] * r->curr[i];
    }
}

int
Num_search_min(Num_data *d, Num_param *p, Num_range *r, float (^model)(float x, float *p, float *dy))
{
    int         i;
    int         n = p->n;
    BOOL        in_range = YES;
    float       cst, val, min_cst;
    Num_param   *ptry = Num_alloc_param(n);

    // rewind
    Num_rewind_range(r);
    for (i = 0; in_range == YES; i++) {
        int j;
        Num_current_param(ptry, r);
        cst = 0;
        for (j = 0; j < d->n; j++) {
            val = model(d->x[j], ptry->data, NULL);
            val -= d->y[j];
            cst += val * val;
        }
        if (cst < min_cst) {
            min_cst = cst;
            Num_current_param(p, r);
        }
        in_range = Num_increment_range(r);
    }
//    printf("min_cst = %f, param: %f, %f\n", min_cst, p->data[0], p->data[1]);

    Num_free_param(ptry);
    return 0;
}
// =====
// least squares using pseudo inverse
// =====
// return val : mse
// data  : in
// param : out (result)
// model : pointer to model function
RecImage *
Num_least_sq(RecImage *data, RecImage *basis, float *mse) // pseudo inverse
{
    int         i, nr, nc;
    Num_mat     *A, *B, *X;
    RecImage    *param;

    nr = [basis yDim];
    nc = [basis xDim];

// alloc
    A = Num_im_to_m(basis);
    B = Num_im_to_m(data);
    X = Num_new_mat(nr, 1);

// pseudo-inverse 
    Num_inv(X, A, B);    // solve AX = B

// set result
    param = Num_m_to_im(X);

// free
    Num_free_mat(A);
    Num_free_mat(B);
    Num_free_mat(X);

    return param;
}

// non linear least squares (conjugate gradient) (keep old one)
int
Num_nonlin_least_sq(Num_data *data, Num_param *param, float (^model)(float x, float *p, float *dy), float *mse)
{
    int     np      = param->n;
    int     nd      = data->n;
    float   *x      = data->x;
    float   *y      = data->y;
    float   *dydp;
    int     iter;

    float   (^cost)(float *);
    void    (^dcost)(float *, float *);

    dydp = (float *)malloc(sizeof(float) * np);

    cost = ^float(float *p) {
        int     i;
        float   val, cst = 0;
        for (i = 0; i < nd; i++) {
            val = model(x[i], p, NULL);
            val -= y[i];
            cst += val*val / nd;    // MSE
        }
        return cst;    
    };
    dcost = ^(float *p, float *gr) {
        int     i, j;
        float   val, e;

        for (j = 0; j < np; j++) {
            gr[j] = 0;
        }
        for (i = 0; i < nd; i++) {
            val = model(x[i], p, dydp);
            e = val - y[i];
            for (j = 0; j < np; j++) {
                gr[j] += 2 * e * dydp[j] / nd;
            }
        }
    };

    iter = Num_conjgr(param, cost, dcost, mse);
    free(dydp);
    return iter;
}

// non linear least squares (conjugate gradient)
int
Num_least_sq_old(Num_data *data, Num_param *param, float (*model)(float x, float *p, float *dy), float *mse)
{
    int     np      = param->n;
    int     nd      = data->n;
    float   *x      = data->x;
    float   *y      = data->y;
    float   *dydp;
	int     iter;

    float   (^cost)(float *);
    void    (^dcost)(float *, float *);

    dydp = (float *)malloc(sizeof(float) * np);

    cost = ^float(float *p) {
        int     i;
        float   val, cst = 0;
        for (i = 0; i < nd; i++) {
            val = model(x[i], p, NULL);
            val -= y[i];
            cst += val*val / nd;    // MSE
        }
        return cst;    
    };
    dcost = ^(float *p, float *gr) {
        int     i, j;
        float   val, e;

        for (j = 0; j < np; j++) {
            gr[j] = 0;
        }
        for (i = 0; i < nd; i++) {
            val = model(x[i], p, dydp);
            e = val - y[i];
            for (j = 0; j < np; j++) {
                gr[j] += 2 * e * dydp[j] / nd;
            }
        }
    };

    iter = Num_conjgr(param, cost, dcost, mse);
    free(dydp);
    return iter;
}

// Conjugate gradient (when gradient is available)
// param  : Num_param struct
// cost   : cost function block. input is param array p
// dcost  : cost gradient function block. input is param array p
// minval : minimum value of cost function (sq error)
int
Num_conjgr(Num_param *param, float (^cost)(float *p), void (^dcost)(float *p, float *dy), float *minval)
{
    int     np      = param->n;
    float   *p      = param->data;
	float	*pmn	= param->min;
	float	*pmx	= param->max;

	int     j,  iter;
	float   gg, gam, fp, dgg;
	float   *g, *h, *xi;
    float   err;
    float   ftol = EPS; //1.0e-10;
    float   fr;

// use conjugate gradient to minimize cost()
	g  = VECTOR(np);
	h  = VECTOR(np);
	xi = VECTOR(np);

    fp = cost(p);   // in: param, out: func value
    dcost(p, xi);
	for (j = 0; j < np; j++) {
		g[j] = -xi[j];
		xi[j] = h[j] = g[j];
	}
	for (iter = 0; iter < ITMAX; iter++) {
        if (iter < 0) {
            fr = 0.5;   // somewhat improve robustness, but not much
        } else {
            fr = 1.0;
        }
		Num_linmin(p, pmn, pmx, xi, np, &err, cost);

if (NumMinimization_dbg) {
//    printf("%d %f\n", iter, err);
    printf("%f %f\n", p[2], p[3]);
}
		if (2.0 * fabs(err - fp) <= ftol * (fabs(err) + fabs(fp) + EPS)) {
            break;
		}
		fp = cost(p);
		dcost(p, xi);
		dgg = gg = 0.0;
		for (j = 0; j < np; j++) {
			gg += g[j] * g[j];
			dgg += (xi[j] + g[j]) * xi[j];
		}
		if (gg == 0.0) {
            break;
		}
		gam = dgg / gg;
		for (j = 0; j<  np; j++) {
			g[j] = -xi[j];
			xi[j] = h[j] = g[j] + gam * h[j];
		}
	}
    if (iter >= ITMAX) {
        printf("NumKit: Too many iterations in Num_conjgr\n");
		iter = -1;
    }

    free(g); free(h); free(xi);
    *minval = err;
    return iter;
}
#define SQR(a) (a) * (a)

// conjugate direction method
//  niter(ret)   out/in          in                      out
int
Num_powell(Num_param *param, float(^cost)(float *p), float *minval)
{
    int     iter = 0;
    int     n = param->n;
    int     i, j, ibig;
    float   del, fp, fptt, t;
    float   *p = param->data;
    float   *pmn = param->min;
    float   *pmx = param->max;
    float   *pt;
    float   *ptt;
    float   *xit;
    float   ftol = EPS; //1.0e-10;
    static float   **xi;   // direction vectors (n x n)

    pt  = VECTOR(n);
    ptt = VECTOR(n);
    xit = VECTOR(n);
    xi = (float **)malloc(sizeof(float *) * n);
    for (i = 0; i < n; i++) {
        xi[i] = VECTOR(n);
    }
    // init directions
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
			xi[i][j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
		xi[i][i] = 0.01;	// not critical, but value matters
    }

    *minval = cost(p);

    for (i = 0; i < n; i++) {
        pt[i] = p[i];
    }
    for (iter = 0; iter < ITMAX; iter++) {
        fp = *minval;
        ibig = 0;
        del = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                xit[j] = xi[j][i];
            }
            fptt = *minval;
            Num_linmin(p, pmn, pmx, xit, n, minval, cost);
            if (fabs(fptt - *minval) > del) {
				del = fabs(fptt -  *minval);
				ibig = i;
            }
        }
        if (2 * fabs(fp - *minval) <= ftol * (fabs(fp) + fabs(*minval))) {
            break;
        }
        for (i = 0; i < n; i++) {
            ptt[i] = 2 * p[i] - pt[i];
            xit[i] = p[i] - pt[i];
            pt[i] = p[i];
        }
        fptt = cost(p);
        if (fptt < fp) {
			t = 2.0 * (fp - 2.0 * (*minval) + fptt) * SQR(fp - (*minval) - del) - del * SQR(fp - fptt);
			if (t < 0.0) {
                Num_linmin(p, pmn, pmx, xit, n, minval, cost);
				for (i = 0; i < n; i++) {
					xi[i][ibig] = xi[i][n - 1];
					xi[i][n - 1] = xit[i];
				}
			}
        }
    }
    free(pt);
    free(ptt);
    xit = VECTOR(n);
    for (i = 0; i < n; i++) {
        free(xi[i]);
    }
    free(xi);

    return iter;
}

#undef SQR

// debug (tmp)
void
t_dump_mat(float *m, int n)
{
    int     i, j;
    printf("===\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%8.5f ", m[i * n + j]);
        }
        printf("\n");
    }
}

void
t_dump_vec(float *v, int n)
{
    int i;
    printf("===\n");
    for (i = 0; i < n; i++) {
        printf("%8.5f ", v[i]);
    }
    printf("\n");
}

// Marquardt-Levenberg
int
Num_marquardt(Num_data *d, Num_param *p, float (*model)(float x, float *p, float *dy), float *chisq)
{
    int     iter;
    float   dchi, lmd;   // delta-chi-square, lambda

    lmd = -1;   // init
    //call Num_mrqmin until convergence is reached
    for (iter = 0; iter < 100; iter++) {
        dchi = Num_mrqmin(d, p, chisq, model, &lmd);
//printf("dchi = %e, chisq = %f, lambda = %f\n", dchi, sqrt(*chisq) * 100, lmd);
    //    if (fabs(dchi) < 1.0e-6) break;
        if (dchi >= 0) break;
    }
    // finalize
    lmd = 0;
    Num_mrqmin(d, p, chisq, model, &lmd);

    return iter;
}

// M-L single iteration
//                           work area           a           X2          model()           lambda
float
Num_mrqmin(Num_data *d, Num_param *p, float *chisq, float (*model)(float x, float *p, float *dy), float *lmd)
{
	int             j, k;
    int             np = p->n;
    static float    *beta, *da, *covar, *alpha;  // work area
    static Num_param   *atry;
	static float    ochisq;
    Num_mat         amat, bmat;
    float           dlt;

    amat.nc = amat.nr = np;
    bmat.nr = np;
    bmat.nc = 1;

    // init
	if (*lmd < 0) {
        // calloc
        covar   = VECTOR(np * np);
        alpha   = VECTOR(np * np);
        beta    = VECTOR(np);
        da      = VECTOR(np);
        atry    = Num_alloc_param(np);
    
		*lmd = 0.1; //0.001;
		for (j = 0; j < np; j++) {
            atry->data[j] = p->data[j];
        }

        // calc Hessian (alpha), gradient (beta), chisq
        //         in   out    out   out    in
		ochisq = *chisq = Num_mrqcof(d, atry, alpha, beta, model);
	}
    // final call, cleanup
	if (*lmd == 0.0) {
        free(covar);
        free(alpha);
        free(da);
        free(beta);
        Num_free_param(atry);
        Num_mrqcof(NULL, atry, alpha, beta, model);
		return -1;
	}

    // main
	for (j = 0; j < np; j++) {
        for (k = 0; k < np; k++) {
            covar[j * np + k] = alpha[j * np + k];
        }
        covar[j * np + j] = alpha[j * np + j] * (1.0 + *lmd);
	}

    // copy covar/beta to Num_mat/Num_vec first
    amat.data = covar;
    bmat.data = beta;
    Num_gaussj(&amat, &bmat);
    for (j = 0; j < np; j++) {
        da[j] = beta[j];
    }
    for (j = 0; j < np; j++) {
        atry->data[j] = p->data[j] + da[j];
    }

    // calc Hessian (covar), gradient (da), chisq
    *chisq = Num_mrqcof(d, atry, covar, da, model);
//    dlt = (*chisq - ochisq) / ochisq;
    dlt = *chisq - ochisq;

    // if chi-square decreases
	if (*chisq < ochisq) {
		*lmd *= 0.1;
		ochisq = *chisq;
		for (j = 0; j < np; j++) {
            for (k = 0; k < np; k++) {
                alpha[j * np + k] = covar[j * np + k];
            }
            beta[j] = da[j];
            p->data[j] = atry->data[j];
		}
    // if chi-square increases
	} else {
		*lmd *= 10.0;
		*chisq = ochisq;
	}

    return dlt;
}

// compute Hessian matrix, gradient and chi-square
float
Num_mrqcof(Num_data *data, Num_param *ptry, float *alpha, float *beta, float (*model)(float x, float *pp, float *ddy))
{
	int     i, j, k;
	float   ymod, wt, dy, chisq;
    int     np = ptry->n;
    static float    *dydp = NULL;

    // init work area
    if (dydp == NULL) {
        dydp = (float *)malloc(sizeof(float) * np);
    }
    // clean up
    if (data == NULL) {
        free(dydp);
        return 0;
    }
    // clear beta, lower left of alpha mat
	for (j = 0; j < np; j++) {
		beta[j] = 0;
		for (k = 0; k <= j; k++) {
            alpha[j * np + k] = 0;
        }
	}
	chisq = 0;
	for (i = 0; i < data->n; i++) {
        ymod = model(data->x[i], ptry->data, dydp);
		dy = data->y[i] - ymod; 
		for (j = 0; j < np; j++) {
            wt = dydp[j];
            beta[j] += dy * wt;
            for (k = 0; k <= j; k++) {
                alpha[j * np + k] += wt * dydp[k];
            }
		}
		chisq += dy * dy;
	}
//t_dump_mat(alpha, np);

    // fill upper right of alpha mat
	for (j = 1; j < np; j++) {
		for (k = 0; k < j; k++) {
            alpha[k * np + j] = alpha[j * np + k];
        }
    }
//t_dump_mat(alpha, np);

    return chisq;
}

// if no range is given, returns -1
// else returns max delta to hit mn/mx
float
Num_max_delta(float *p, float *pmn, float *pmx, float *xi, int np)
{
	float	mx_dlt;
	int		i;
	float	d;

	mx_dlt = -1;
	if (pmn[0] == pmx[0]) return mx_dlt;

	for (i = 0; i < np; i++) {
		if (xi[i] == 0) {
			continue;
		}
		if (xi[i] > 0) {	// chk max
			d = (pmx[i] - p[i]) / xi[i];
		} else {			// chk min
			d = (pmn[i] - p[i]) / xi[i];
		}
		if (mx_dlt < 0) {
			mx_dlt = d;
		} else
		if (mx_dlt > d) {
			mx_dlt = d;
		}
	}

	return mx_dlt;
}

// goal: ax < bx < cx, f(ax) > f(bx) < f(cx)
// start position is within range & ax = 0 -> f(bx)
//
void
Num_init_brak(float *ax, float *bx, float *cx, float max_delta, float (^f1d)(float))
{
    float   fa, fb, fc;

    fa = f1d(*ax);
    fb = f1d(*bx);
	*cx = *bx;
	fc = fb;
    if (fb > fa) { // search inside
		while (fb > fa) {
			*bx = (*ax + *bx) / 2;
			fb = f1d(*bx);
		}
    } else { // else search outside of b
		while (fc <= fb && fc < max_delta) {
			*cx = *cx + (*cx - *ax) * 2;
			fc = f1d(*cx);
		}
    }
}

// bracketing...
// make f(ax) > f(bx) < f(cx) : ax < bx < cx
// if f(cx) == NaN, this returns as OK...
// mn <= ax < bx < cx <= mx
// === old version... replaced by Num_init_brak()
void
Num_mnbrak(float *ax, float *bx, float *cx, float (^f1d)(float))
{
    float   ulim, u, r, q, fu, dum;
    float   fa, fb, fc;

    fa = f1d(*ax);
    fb = f1d(*bx);
    if (fb > fa) {
        SHFT(dum, *ax, *bx, dum);
        SHFT(dum, fb, fa, dum);
    }
    // ### should check inside first...
    // ### first outside point is out of bounds for non-scaled version
    *cx = *bx + GOLD * (*bx - *ax);
    fc = f1d(*cx);
//printf("a:%f(%f) b:%f(%f) c:%f(%f)\n", *ax, fa, *bx, fb, *cx, fc);
    while (fb > fc) {   // not found ... extend
		r = (*bx - *ax) * (fb - fc);
		q = (*bx - *cx) * (fb - fa);
		u = *bx - ((*bx - *cx) * q - (*bx - *ax) * r)
			/ (2.0 * SIGN(MAX(fabs(q - r), TINY), q - r));
		ulim = *bx + GLIMIT * (*cx - *bx);
		if ((*bx - u) * (u - *cx) > 0.0) {
			fu = f1d(u);
			if (fu < fc) {
				*ax = *bx;
				*bx = u;
				fa = fb;
				fb = fu;
				return;
			} else
            if (fu > fb) {
				*cx = u;
				fc = fu;
				return;
			}
			u = *cx + GOLD * (*cx - *bx);
			fu = f1d(u);
		} else
        if ((*cx - u)*(u - ulim) > 0.0) {
			fu = f1d(u);
			if (fu < fc) {
				SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx));
				SHFT(fb, fc, fu, f1d(u));
			}
		} else
        if ((u - ulim) * (ulim - *cx) >= 0.0) {
			u = ulim;
			fu = f1d(u);
		} else {
			u = *cx + GOLD * (*cx - *bx);
			fu = f1d(u);
		}
		SHFT(*ax, *bx, *cx, u);
		SHFT(fa, fb, fc, fu);
    }
}

// debug ???
float
Num_chk_range(float *p, float *pmin, float *pmax, float *gr, int np)
{
    float   bx = 1.0;
    float   d;
    int     i;

printf("========\n");
    for (i = 0; i < np; i++) {
        if (gr[i] < 0) {
            d = (pmin[i] - p[i]) / gr[i];
        } else {
            d = (pmax[i] - p[i]) / gr[i];
        }
        if (i == 0) {
            bx = d;
        } else {
            bx = MIN(bx, d);
        }
        printf("gr[%d] = %e, range: %e : %e : %e [%e]\n", i, gr[i], pmin[i], p[i], pmax[i], d);
    }
    bx *= 0.5;
    printf("=== bx = %e ===\n", bx);
    return bx;
}

//
// Brent minimization (NR)
// input: ax < bx < cx, f(ax) > f(bx) < f(cx)
// output: xmin, returns f(xmin)
//
float
Num_brent(float ax, float bx, float cx, float (^f1d)(float), float tol, float *xmin)
{
	int     iter;
	float   a, b, d, etemp;
    float   fu, fv, fw, fx, p, q, r;
    float   tol1, tol2, u, v, w,x , xm;
	float   e = 0.0;

//	a = MIN(ax, cx);
//	b = MAX(ax, cx);
	a = ax;
	b = cx;
	x = w = v = bx;
    // a - x - b
	fw = fv = fx = f1d(x);
	for (iter = 0; iter < ITMAX; iter++) {
		xm = (a + b) / 2;
		tol2 = 2.0 * (tol1 = tol * fabs(x) + EPS);
		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0) p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			} else {
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2) {
					d = SIGN(tol1, xm - x);
                }
			}
		} else {
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = f1d(u);
		if (fu <= fx) {
			if (u >= x) a = x; else b = x;
			SHFT(v, w, x, u);
			SHFT(fv, fw, fx, fu);
		} else {
			if (u < x) a = u; else b = u;
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			} else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}
//	printf("NumKit: Too many iterations in brent");
	*xmin = x;

	return fx;
}

// p[]: param, in/out, xi[]: gradient, np: nparam, fret:val return, func():in
// pmn, pmx: in, param range (added on 5-27-2016)
void
Num_linmin(float *p, float *pmn, float *pmx, float *xi, int np, float *fret, float (^func)(float *))
{
	int             i;
	float           xmin, ax, bx, cx;
	float			max_delta;
    __block float   *xtmp;
    float           (^f1d)(float);

    xtmp = VECTOR(np);
    f1d = ^float(float x) { // var = ^returntype (argtype arg)
        int     j;
        for (j = 0; j < np; j++) {
            xtmp[j] = p[j] + x * xi[j];
        }
        return func(xtmp);
    };

	/**/
	max_delta = Num_max_delta(p, pmn, pmx, xi, np);;
//printf("mx_delta = %f\n", max_delta);

	ax = 0.0;
	bx = 1.0;	// initial range is hardcoded (Num Recipe)

	Num_mnbrak(&ax, &bx, &cx, f1d);
//	Num_init_brak(&ax, &bx, &cx, max_delta, f1d);	// ok for conjgr, NG for powel -> chk

//printf("mnbrak:  a(%4.2f) %4.3e \t- b(%4.2f) %4.3e \t- c(%4.2f) %4.3e\n",
//	ax, f1d(ax), bx, f1d(bx), cx, f1d(cx));

	*fret = Num_brent(ax, bx, cx, f1d, BRTOL, &xmin);
//printf("brent:x = %g, y = %g\n", xmin, *fret);

	for (i = 0; i < np; i++) {
		p[i] += xi[i] * xmin;
	}
    free(xtmp);
}

Num_data *
Num_alloc_data(int n)
{
    Num_data    *data = (Num_data *)malloc(sizeof(Num_data));

    data->n     = n;
    data->x     = VECTOR(n);
    data->y     = VECTOR(n);
    data->xscale = data->yscale = 1.0;
    return data;
}

Num_param *
Num_alloc_param(int n)
{
    Num_param   *param = (Num_param *)malloc(sizeof(Num_param));
    int         i;
    param->n = n;
    param->data = VECTOR(n);
    param->min = VECTOR(n);
    param->max = VECTOR(n);
    for (i = 0; i < n; i++) {
        param->min[i] = 0.0;
        param->max[i] = 1.0;
    }

//    param->flag = (int *)malloc(sizeof(int) * n);

    return param;
}

Num_range *
Num_alloc_range(int n)
{
    Num_range   *range = (Num_range *)malloc(sizeof(Num_range));
    range->n = n;
    range->min = VECTOR(n);
    range->inc = VECTOR(n);
    range->steps = (int *)malloc(sizeof(int) * n);
    range->curr  = (int *)malloc(sizeof(int) * n);
    return range;
}

void
Num_free_data(Num_data *d)
{
    if (d) {
        free(d->x);
        free(d->y);
        free(d);
    }
}

void
Num_free_param(Num_param *p)
{
    if (p) {
        free(p->data);
		free(p->min);
		free(p->max);
    //    free(p->flag);
        free(p);
    }
}

void
Num_free_range(Num_range *r)
{
    if (r) {
        free(r->min);
        free(r->inc);
        free(r->steps);
        free(r->curr);
        free(r);
    }
}

void
Num_normalize_data(Num_data *d)
{
    int     i;
    float   val, xs, ys;

    xs = ys = 0;
    for (i = 0; i < d->n; i++) {
        val = fabs(d->x[i]);
        if (xs < val) xs = val;
        val = fabs(d->y[i]);
        if (ys < val) ys = val;
    }
    for (i = 0; i < d->n; i++) {
        d->x[i] /= xs;
        d->y[i] /= ys;
    }
    d->xscale = xs;
    d->yscale = ys;
}

void
Num_normalize_y(Num_data *d)
{
    int     i;
    float   val, ys;

    ys = 0;
    for (i = 0; i < d->n; i++) {
        val = fabs(d->y[i]);
        if (ys < val) ys = val;
    }
    for (i = 0; i < d->n; i++) {
        d->y[i] /= ys;
    }
    d->xscale = 1;
    d->yscale = ys;
}

void
Num_dump_data(Num_data *d)
{
    int     i;
    printf("Num_data (n = %d, xscale:%f, yscale:%f)\n", d->n, d->xscale, d->yscale);
    for (i = 0; i < d->n; i++) {
        printf("%f %f\n", d->x[i], d->y[i]);
    }
}

Num_data  *
Num_copy_data(Num_data *d)
{
    Num_data    *d2;
    int         i;

    d2 = Num_alloc_data(d->n);
    for (i = 0; i < d->n; i++) {
        d2->x[i] = d->x[i];
        d2->y[i] = d->y[i];
    }
    d2->xscale = d->xscale;
    d2->yscale = d->yscale;

    return d2;
}

void
Num_dump_param(Num_param *p)
{
    int     i;
    printf("Num_param (min/max) = %f / %f\n", p->min, p->max);
    for (i = 0; i < p->n; i++) {
        printf("%d %f\n", i, p->data[i]);
    }
}

// ========== gaussian amoeba ========

// make this general purpose (dim != 4)
typedef struct {
    float   err;    // square error
    int     dim;
    float   *x;     // param
} gs_pt;

void
//next_gs_pt(gs_pt *pt, Num_param *param, Num_data *data, NumMatrix *cosd, float (^model_b)(Num_param *prm, float x))
next_gs_pt(gs_pt *pt, Num_param *param, NumMatrix *cosd, float (^model_b)(Num_param *prm))
{
    NumMatrix   *dp;
    Num_param   *prm = Num_alloc_param(param->n);
    float       *p, t, val, err;
    int         i;

// generate step vector with cov
    dp = [cosd matNrml];
    p = [dp data];
    for (i = 0; i < param->n; i++) {
        val = param->data[i] + p[i];
        if (prm->min[i] > val) {
            prm->data[i] = prm->min[i];
        } else
        if (prm->max[i] < val) {
            prm->data[i] = prm->max[i];
        } else {
            prm->data[i] = val;
        }
    }
    pt->err = model_b(prm);
    for (i = 0; i < param->n; i++) {
        pt->x[i] = prm->data[i];
    }
    Num_free_param(prm);
}

void
copy_gs_pt(gs_pt *pt1, gs_pt *pt2)
{
    int     i;
    pt2->err = pt1->err;
    for (i = 0; i < pt1->dim; i++) {
        pt2->x[i] = pt1->x[i];
    }
}

// ### mn[4] -> *mn
NumMatrix *
calcCov(gs_pt *buf, int bufLen)
{
    int         i, j, k, dim;
    float       var, *mn;   // var, mean
    NumMatrix   *cov;
    float       *p;

    dim = buf[0].dim;
    cov = [NumMatrix matrixOfType:NUM_REAL nRow:dim nCol:dim];
    p = [cov data];
    mn = (float *)malloc(sizeof(float) * dim);
    for (j = 0; j < dim; j++) {
        // calc mean
        mn[j] = 0;
        for (i = 0; i < bufLen; i++) {
            mn[j] += buf[i].x[j];
        }
        mn[j] /= bufLen;
    }
// calc covar
    var = 0;
    for (j = 0; j < dim; j++) {   // row
        for (k = j; k < dim; k++) {   // col
            var = 0;
            for (i = 0; i < bufLen; i++) {
                var += (buf[i].x[j] - mn[j]) * (buf[i].x[k] - mn[j]);
            }
            var /= bufLen;
            p[k*dim + j] = var;
            if (i != k) {
                p[j*dim + k] = var;
            }
        }
    }
    free(mn);
    return cov;
}

void
init_gs_pt(gs_pt *pt, int dim)
{
    int     i;
        
    pt->err = -1;
    pt->dim = dim;
    pt->x = (float *)malloc(sizeof(float) * dim);
    for (i = 0; i < dim; i++) {
        pt->x[i] = 0;
    }
}

// expects scaled inputs
int
Num_gauss_amoeba(Num_param *param, float (^model_b)(Num_param *prm), float *mse)
{
    int         i, j;
    int         dim;    // number of parameters
    int         iter, max_iter, lap;
    gs_pt       *gbuf, pt;
    int         buf_len;
    int         min_ix, max_ix;
    float       min_err, max_err;
    NumMatrix   *cov, *cosd;   // current covariance / coSD of amoeba
    Num_param   *next_prm;

    buf_len = 30;
    dim = param->n;
    max_iter = 100000; //100000;
    cov = [NumMatrix unitMatrixOfDim:dim];
    cov = [cov multByConst:0.001];     // initial cov
    cosd = [cov matSqrt];
    next_prm = Num_alloc_param(dim);

    // alloc buffer
    gbuf = (gs_pt *)malloc(sizeof(gs_pt) * buf_len);
    init_gs_pt(&pt, dim);
    
    // fill buffer (initialize)
    for (i = 0; i < buf_len; i++) {
        // alloc x
        init_gs_pt(gbuf + i, dim);
        
        // gen rand pt, eval error, and fill buffer
        next_gs_pt(&pt, param, cosd, model_b);
        
        copy_gs_pt(&pt, gbuf + i);
//    printf("%d %f\n", i, pt.err);
    }

    // find min/max pt ###
    min_ix = max_ix = 0;
    min_err = max_err = gbuf[0].err;
    for (i = 1; i < buf_len; i++) {
//        printf("%d %f\n", i, gbuf[i].err);
        if (min_err > gbuf[i].err) {
            min_err = gbuf[i].err;
            min_ix = i;
        }
        if (max_err < gbuf[i].err) {
            max_err = gbuf[i].err;
            max_ix = i;
        }
    }

    // move center to current min (cov is I)
    for (j = 0; j < dim; j++) {
        param->data[j] = gbuf[min_ix].x[j];
    }            

    // === main amoeba ====
    lap = 0;
    for (iter = 0; iter < max_iter; iter++) {
        if (gbuf[max_ix].err - gbuf[min_ix].err < 1.0e-8) break; // 1.0e-8
        next_gs_pt(&pt, param, cosd, model_b);
        if (pt.err > max_err) {
            continue;
        }
        for (j = 0; j < dim; j++) {
            if (pt.x[j] <= param->min[j] || pt.x[j] >= param->max[j]) {
                continue;
            }
        }
        if (pt.err >= min_err) { // within range
            // replace max
            copy_gs_pt(&pt, gbuf + max_ix);
            // find max
            max_ix = 0;
            max_err = gbuf[0].err;
            for (i = 0; i < buf_len; i++) {
                if (max_err < gbuf[i].err) {
                    max_err = gbuf[i].err;
                    max_ix = i;
                }
            }
        } else {    // new min
            copy_gs_pt(&pt, gbuf + max_ix);  // replace max with new min
            min_ix = max_ix;    //
            min_err = pt.err;
            // find max
            max_ix = 0;
            max_err = gbuf[0].err;
            for (i = 0; i < buf_len; i++) {
               if (max_err < gbuf[i].err) {
                   max_err = gbuf[i].err;
                   max_ix = i;
               }
            }
            // update current pt
            // move center to current min
            for (j = 0; j < dim; j++) {
                param->data[j] = gbuf[min_ix].x[j];
            }            
            // update cosd
            cov = calcCov(gbuf, buf_len);
            cosd = [cov matSqrt];

            if (NumMinimization_dbg) {
                printf("%f %f\n", param->data[2], param->data[3]);   // tc1, tc2
            //    printf("%d %e %e %e\n", iter, min_err, max_err, max_err - min_err);
            }
        }
    }
    if (gbuf) {
        for (i = 0; i < buf_len; i++) {
            if (gbuf[i].x) {
                free(gbuf[i].x);
            }
        }
        free(gbuf);
    }
    *mse = pt.err;

    return iter;
}

// bi-exponential curve fitting
// input is not altered (copy is scaled), param is output (in orig scale)
int
Num_biex_fit(Num_param *param, Num_data *data, float *mse)
{
    Num_data    *p_dat;     // log
    Num_param   *p_prm;     // polynomial coeff
    Num_data    *e_dat;     // linear
    Num_param   *e_prm;     // bi-ex coeff
    float       a, t;
    int         i, j, niter;
    int         order = 3;  // 5
    float       r1, r2, c;
    float       (^model_b)(Num_param *prm);
    float       mn, mx;

    mn = mx = data->y[0];
    for (i = 0; i < data->n; i++) {
        if (mn > data->y[i]) {
            mn = data->y[i];
        }
        if (mx < data->y[i]) {
            mx = data->y[i];
        }
    }
    if (mx - mn < 1.0e-6) {
    //    printf("range\n");
        return 0;
    }
    p_dat = Num_copy_data(data);    // alloc and copy
    Num_normalize_data(p_dat);
    e_dat = Num_copy_data(p_dat);

    // take log
    for (i = 0; i < p_dat->n; i++) {
        p_dat->y[i] = log(p_dat->y[i]);
    }

// ======== initial guess using polynomial expantion / least squares =======
//  polynomial fitting
    p_prm = Num_alloc_param(order);
    Num_poly_fit(p_prm, p_dat);

// expand -> continuous
    for (i = 0; i < 100; i++) {
        a = 0;
        t = (float)i / 100;
        for (j = 0; j < order; j++) {
            a += p_prm->data[j] * pow(t, j);
        }
        a = exp(a);
    }

    e_prm = Num_alloc_param(4);

// first estimate tc assuming frac = 0.5
// expand
//    a_to_tau(p_prm, e_prm);
    c = 2.0 * p_prm->data[2]; 
    if (c < 0) c = 0; // single
    c = sqrt(c);

    r1 = -p_prm->data[1] + c;
    r2 = -p_prm->data[1] - c;
    if (r2 <= 0) {
        e_prm->data[3] = 100;
    } else {
        e_prm->data[3] = 1.0 / r2;
    }
    e_prm->data[2] = 1.0 / r1;
    e_prm->data[0] = e_prm->data[1] = exp(p_prm->data[0])/2;

// range
    e_prm->min[0] = 0; e_prm->max[0] = 5.0;    // pd1
    e_prm->min[1] = 0; e_prm->max[1] = 5.0;    // pd2
    e_prm->min[2] = 0; e_prm->max[2] = 2.0;    // tc1
    e_prm->min[3] = 0; e_prm->max[3] = 2.0;    // tc2

// === adaptive gaussian search (general purpose) ===
    model_b = ^float(Num_param *prm) {
                    float   a1, a2, t1, t2, x, y, err;
                    int     i;
                    a1 = prm->data[0];
                    a2 = prm->data[1];
                    t1 = prm->data[2];
                    t2 = prm->data[3];

                    err = 0;
                    for (i = 0; i < e_dat->n; i++) {
                        if (a1 <= 0 || a2 <= 0 || t1 <= 0 || t2 <= 0) {
                            return 100.0;    // error return
                        } else {
                            x = e_dat->x[i];
                            y = a1 * exp(-x / t1) + a2 * exp(-x / t2);
                        }
                        err += (y - e_dat->y[i]) * (y - e_dat->y[i]);
                    }
                    return sqrt(err);
                };
    // === main func call ===
    niter = Num_gauss_amoeba(e_prm, model_b, mse);

    if (e_prm->data[2] > e_prm->data[3]) {
        float   tmp;
        tmp = e_prm->data[2];
        e_prm->data[2] = e_prm->data[3];
        e_prm->data[3] = tmp;
        tmp = e_prm->data[0];
        e_prm->data[0] = e_prm->data[1];
        e_prm->data[1] = tmp;
    }

    param->data[0] = e_prm->data[0] * p_dat->yscale;
    param->data[1] = e_prm->data[1] * p_dat->yscale;
    param->data[2] = e_prm->data[2] * p_dat->xscale;
    param->data[3] = e_prm->data[3] * p_dat->xscale;

    Num_free_param(p_prm);
    Num_free_param(e_prm);
    Num_free_data(p_dat);
    Num_free_data(e_dat);
    
    return niter;
}

// single exponential curve fitting
int
Num_exp_fit(Num_param *param, Num_data *data, float *mse)
{
    Num_data    *e_dat;     // linear
    Num_param   *e_prm;     // bi-ex coeff
    int         niter;
    float       (^model_b)(Num_param *prm);

    e_dat = Num_copy_data(data);    // alloc and copy
    Num_normalize_data(e_dat);

    e_prm = Num_alloc_param(2);
    e_prm->data[0] = 1.0;
    e_prm->data[1] = 0.5;

// gaussian amoeba (general purpose) ### -> prm, x -> y
    model_b = ^float(Num_param *prm) {
                    float   pd, tc, x, y, err;
                    int     i;
                    pd  = prm->data[0];
                    tc = prm->data[1];

                    err = 0;
                    for (i = 0; i < e_dat->n; i++) {
                        if (pd <= 0 || tc <= 0) {
                            return 100.0;    // error return
                        } else {
                            x = e_dat->x[i];
                            y = pd * exp(-x / tc);
                        }
                        err += (y - e_dat->y[i]) * (y - e_dat->y[i]);
                    }
                    return sqrt(err);
                };

Num_dump_param(e_prm);
    niter = Num_gauss_amoeba(e_prm, model_b, mse);
Num_dump_param(e_prm);

    param->data[0] = e_prm->data[0] * e_dat->yscale;
    param->data[1] = e_prm->data[1] * e_dat->xscale;

    Num_free_param(e_prm);
    Num_free_data(e_dat);
    
    return niter;
}


// move to NumKit when done
float
Num_poly_fit(Num_param *param, Num_data *data)
{
    float       err;
    RecImage    *X, *B, *A;
    float       *px, *pb, *pa;
    int         i, j, ndata, order;

    ndata = data->n;
    order = param->n;

    A = [RecImage imageOfType:RECIMAGE_REAL xDim:order yDim:ndata];
    B = [RecImage imageOfType:RECIMAGE_REAL xDim:1 yDim:ndata];
    pa = [A data];
    pb = [B data];
    for (i = 0; i < ndata; i++) {
        pb[i] = data->y[i];
        pa[i*order] = 1.0;
        for (j = 1; j < order; j++) {
            pa[i*order + j] = pa[i*order + j - 1] * data->x[i];
        }
    }

    // solve
    X = Num_least_sq(B, A, &err); // pseudo inverse

    // copy output to param
    px = [X data];
    for (i = 0; i < order; i++) {
        param->data[i] = px[i];
    }

    return err;
}


#undef  GOLD
#undef  GLIMIT
#undef  TINY
#undef  EPS
#undef  BRTOL
#undef  SHFT
#undef  SIGN

#undef  ITMAX
#undef  CGOLD
#undef  VECTOR
