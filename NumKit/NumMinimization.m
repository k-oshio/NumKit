//
// NumMinimization.m
//
//  == plans ==
//	implement param limits (started on 5-26-2016)
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
#define ITMAX               100     //800	// orig:100

// todo: ... global search

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

// data  : in
// param : in (initial est) / out (result)
// model : pointer to model function. input is indep var x, param array pp, and dy/dp array.
// minval: minimum value of cost function (sq error).
int
Num_least_sq(Num_data *data, Num_param *param, float (^model)(float x, float *p, float *dy), float *mse)
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
		//printf("iter = %d, cost = %f, p0:%f p1:%f err:%f\n", iter, fp, p[0], p[1], err);
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
//    if (iter >= ITMAX) {
//        printf("NumKit: Too many iterations in Num_conjgr\n");
//    }

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
    Num_param    *param = (Num_param *)malloc(sizeof(Num_param));
    param->n = n;
    param->data = VECTOR(n);
    param->min = VECTOR(n);
    param->max = VECTOR(n);
    param->flag = (int *)malloc(sizeof(int) * n);
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
        free(p->flag);
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
