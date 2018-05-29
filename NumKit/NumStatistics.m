//
//  NumStatistics
//

#import "NumKit.h"


#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define setmin(x); {if (fabs(x) < FPMIN) x = FPMIN;}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void
Num_error(char *msg)   // error exit
{
    printf("%s", msg);
//   exit(-1);
}

float
Num_mean(float *p, int n)
{
    float   mn;
    vDSP_sve(p, 1, &mn, n);
    return mn / n;
}

float
Num_var(float *p, int n, float m)
{
    float   v;
    float   s, ep;
    int     i;

// variance
    v = ep = 0;
    for (i = 0; i < n; i++) {
        s = p[i] - m;
        ep += s;
        v += s * s;
    }
    v = (v - ep * ep/n) / (n - 1);  // ep*ep/n is correction term for rounding error
    return v;
}

// mean / var estimator
// -> use Accelerate framework
// p[] is not modified (no need to copy)
void
Num_avevar(float *p, int n, float *ave, float *var)
{
    float   m, v;

    m = Num_mean(p, n);		// sum(p) / n;
    v = Num_var(p, n, m);	// sum(p - m) / (n - 1)

    *ave = m;
    *var = v;
}

// x:   real,       vec(n_sample)
// y:   imag,       vec(n_sample)
// av:  mean,       vec(2)
// cov: covariance, vec(3) : xx, yy, xy
void
Num_covar(float *x, float *y, int n, float *av, float *cov)
{
    float   mnx, mny, vx, vy, vxy;
    float   xx, yy;
    int     i;

    if (n < 2) return;

    mnx = Num_mean(x, n);
    mny = Num_mean(y, n);
    vx = vy = vxy = 0;
    for (i = 0; i < n; i++) {
        xx = x[i] - mnx;
        yy = y[i] - mny;
        vx += xx * xx;
        vy += yy * yy;
        vxy += xx * yy;
    }
    vx  /= n;
    vy  /= n;
    vxy /= n;
    
    av[0] = mnx;
    av[1] = mny;
    cov[0] = vx;
    cov[1] = vy;
    cov[2] = vxy;
}

// not tested yet
void
Num_linreg(float *r, float *b0, float *b1, float *x, float *y, int skip, int n)
{
    float   mx, my, sx, sy;
    float   sxy, rr, a, b;
    int     i;

//    Num_avevar(x, n, &mx, &sx);
//    Num_avevar(y, n, &my, &sy);
    mx = Num_mean(x, n);
    sx = Num_var(x, n, mx);
    my = Num_mean(y, n);
    sy = Num_var(y, n, mx);
    sxy = 0;
    for (i = 0; i < n; i++) {
        sxy += x[i] * y[i];
    }
    rr = (sxy - n * mx * my) / ((n - 1) * sx * sy);
    a = sy * rr / sx;
    b = my - a * mx;

    *r = rr;
    *b0 = b;
    *b1 = a;
}

// Gaussian, N(0.0, 1.0)
float
Num_nrml(float m, float sd)
{
	int		i, n = 10;
	float	r = 0;
    float   mx = pow(2.0, 31) - 1;

	for (i = 0; i < n; i++) {
		r += (float)random() / mx;
	}
	r = (r - (float)n/2) * sqrt(12.0/n);

	return r * sd + m;
}// uniform dist, (0.0, 1.0)
float
Num_unif(float m, float sd)
{
	float	r;
    float   mx = pow(2.0, 31) - 1;

	r = (float)random() / mx;	// [0..1]

	return (r - 0.5) * 2 * sqrt(3) * sd + m;
}


// find variation direction of complex random number,
// and rotate it to real axis
void
Num_orth(float *pr1, float *pr2, float *re1, float *im1, float *re2, float *im2, int n, float *th)
{
    float   cs, sn, r;
    float   mr1, mr2;
    float   mi1, mi2;
    float   re, im;
    int     i;

    mr1 = Num_mean(re1, n); mi1 = Num_mean(im1, n);
    mr2 = Num_mean(re2, n); mi2 = Num_mean(im2, n);
    cs = mr2 - mr1;
    sn = mi2 - mi1;
    r = sqrt(cs*cs + sn*sn);
    if (r != 0) {
        cs /= r;
        sn /= r;
    } else {
        cs = 0.707;
        sn = 0.707;
    }

    for (i = 0; i < n; i++) {
        re = re1[i]; im = im1[i];
        pr1[i] = re * cs + im * sn;
    //    im1[i] = re * sn - im * cs;
        re = re2[i]; im = im2[i];
        pr2[i] = re * cs + im * sn;
    //    im2[i] = re * sn - im * cs;
    }
    *th = atan2(sn, cs);
}

// continued beta fraction (used by betai)
float
betacf(float a, float b, float x)
{
    int     m, m2;
    float   aa, c, d, del, h;

    c = 1.0;
    d = 1.0 - (a + b) * x / (a + 1.0); setmin(d);
    d = 1.0 / d;
    h = d;

    for (m = 1; m <= MAXIT; m++) {
        m2 = 2 * m;
        aa = m * (b - m) * x / ((a + m2 - 1.0) * (a + m2));
        d = 1.0 + aa * d; setmin(d);
        c = 1.0 + aa / c; setmin(c);
        d = 1.0 / d;
        h *= d * c;

        aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0));
        d = 1.0 + aa * d; setmin(d);
        c = 1.0 + aa / c; setmin(c);
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }
    if (m > MAXIT) Num_error("a or b too big, or MAXIT too small");

    return h;
}

// regularized incomplete beta function = cumulative beta distribution
float
betai(float a, float b, float x)
{
	float bt;

	if (x < 0.0 || x > 1.0) {
        Num_error("x has to be [0..1]");
        return 0;
    }
	if (x == 0.0 || x == 1.0) {
        bt = 0.0;
	} else {
		bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1.0 - x));
    }
	if (x < (a + 1.0) / (a + b + 2.0)) {
		return bt * betacf(a, b, x) / a;
	} else {
		return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
    }
}

float
beta(float a, float b)				// beta function
{
	if (a * b == 0) {
		return 0;
	} else {
		return tgamma((double)a) * tgamma((double)b) / tgamma((double)a + b);
	}
}

// t distribution
float
Num_tdist(float tval, float df)
{
    float   n12, n2;
    float   a;

    n2 = df / 2.0;
    n12 = (df + 1) / 2.0;
    a = tgamma(n12) / tgamma(df/2.0) / sqrt(df * M_PI);

    return a * pow(1.0 + tval * tval/df, -n12);
}

// F distribution
float
Num_fdist(float fval, float df1, float df2)
{
	double df12 = df1 / df2;

	return 1.0 / beta(df1 / 2, df2 / 2)
			* pow(df12, df1 / 2)
			* pow(fval, df1 / 2 - 1)
			* pow(1 + fval * df12, -(df1 + df2) / 2);
}

// gaussian distribution
float
Num_gdist(float x, float sd)
{
	return 1.0 / sqrt(2 * M_PI) / sd * exp(-(x * x / sd / sd / 2));
}

float
Num_t2p(float tval, float df)
{
    float p;
    p = betai(0.5 * df, 0.5, df / (tval * tval + df));
    return p * 0.5; // one-sided
}

float
Num_tlimit(float p, float df)
{
    float   (^f1d)(float) = ^float(float x) { return Num_t2p(x, df) - p; };

    return Num_zbrent(0, 100, f1d);
}

float
Num_f2p(float fval, float df1, float df2)  // F: var1/var2
{
    float   p;
	p =  2.0 * betai(0.5 * df2, 0.5 * df1, df2 / (df2 + df1 * fval));
	if (p > 1.0) p = 2.0 - p;
    return p * 0.5;   // one-sided
}

float
Num_flimit(float p, float df1, float df2, BOOL left)
{
    float   (^f1d)(float) = ^float(float x) { return Num_f2p(x, df1, df2) - p; };

    if (left) {
        if (df1 < 2 || df2 < 2) return 0.0;
        return Num_zbrent(1, 0, f1d);
    } else {
        if (df1 < 2 || df2 < 2) return 100.0;
        return Num_zbrent(1, 100, f1d);
    }
}

// 1: sample, 2:ref
void
Num_ttest(float *p1, int n1, float *p2, int n2, float *tval, float *pval)
{
    float   m1, m2, v1, v2;
    float   var, serr, df;
    float   t;

    Num_avevar(p1, n1, &m1, &v1);
    Num_avevar(p2, n2, &m2, &v2);
    df = n1 + n2 - 2;   // degree of freedom
    var = ((n1 - 1) * v1 + (n2 - 1) * v2) / df; // variance (estimate)
    serr = sqrt(var * (1.0 / n1 + 1.0 / n2));    // std error
    t = (m1 - m2) / serr;

    *tval = t;
    *pval = Num_t2p(t, df);
}

void
Num_ftest(float *p1, int n1, float *p2, int n2, float *fval, float *pval)
{
	float   m1, m2, v1, v2, df1, df2;
    float   f, p;

    Num_avevar(p1, n1, &m1, &v1);
    Num_avevar(p2, n2, &m2, &v2);

    f = v1 / v2;
    df1 = n1 - 1;
    df2 = n2 - 1;

	p = Num_f2p(f, df1, df2);

    *fval = f;
    *pval = p;
}

float
Num_r2(float *p, float *f, int n)
{
    float   m, v;
    float   sum_tot, sum_err, val;
    int     i;

//    Num_avevar(p, n, &m, &v);
    m = Num_mean(p, n);
    v = Num_var(p, n, m);
    sum_tot = 0;
    for (i = 0; i < n; i++) {
        val = p[i] - m;
        sum_tot += val * val;
    }
    sum_err = 0;
    for (i = 0; i < n; i++) {
        val = p[i] - f[i];
        sum_err += val * val;
    }
    return 1.0 - sum_err / sum_tot;
}

float
Num_zbrent(float ax, float bx, float (^f1d)(float))
{
    float   tol = EPS;
	int     iter;
	float   a = ax, b = bx, c = bx, d, e, min1, min2;
	float   fa = f1d(a), fb = f1d(b), fc, p, q, r, s, tol1, xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
		Num_error("Root must be bracketed in zbrent");
        return 0;
    }

	fc = fb;
	for (iter = 0;iter < MAXIT; iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c  = a;
			fc = fa;
			e = d = b - a;
		}
		if (fabs(fc) < fabs(fb)) {
			a  = b;
			b  = c;
			c  = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
		xm = 0.5 * (c - b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb / fa;
			if (a == c) {
				p = 2.0 * xm * s;
				q = 1.0 - s;
			} else {
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if (p > 0.0) q = -q;
			p = fabs(p);
			min1 = 3.0 * xm * q - fabs(tol1 * q);
			min2 = fabs(e * q);
			if (2.0 * p < (min1 < min2 ? min1 : min2)) {
				e = d;
				d = p / q;
			} else {
				d = xm;
				e = d;
			}
		} else {
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (fabs(d) > tol1) {
			b += d;
		} else {
			b += SIGN(tol1, xm);
        }
		fb = f1d(b);
	}

    Num_error("a or b too big, or MAXIT too small");
	return 0;
}

// === Complex Statistics ===
void
Num_cnrml(float *re, float *im, float mr, float mi, float sd)
{
    float   sd2 = sd / 2;
    *re = Num_nrml(mr, sd2);
    *im = Num_nrml(mi, sd2);
}

void
Num_cmean(float *re, float *im, int n, float *mr, float *mi)
{
	vDSP_meanv(re, 1, mr, n);
	vDSP_meanv(im, 1, mi, n);
}

void
Num_cavevar(float *p, float *q, int n, float *avr, float *avi, float *var)
{
    int     i;
    float   mr, mi, re, im, v;

// mean
	Num_cmean(p, q, n, &mr, &mi);

// variance
    v = 0;
    for (i = 0; i < n; i++) {
        re = p[i]- mr;
        im = q[i]- mi;
        v += re * re + im * im;
    }

    *avr = mr;
    *avi = mi;
    *var = v / (n - 1);
}

// complex mean difference test using F
// 1: sample, 2:ref
void
Num_ctmean(float *p1, float *q1, int n1, float *p2, float *q2, int n2, float *t2val, float *pval, float *th)
{
    float   var, serr2, t2, mr1, mi1, mr2, mi2, v1, v2;
    int     df;

    Num_cavevar(p1, q1, n1, &mr1, &mi1, &v1);
    Num_cavevar(p2, q2, n2, &mr2, &mi2, &v2);

    df = n1 + n2 - 2;
    var = (v1 * (n1 - 1) + v2 * (n2 - 1)) / df;
    serr2 = var * (1.0 / n1 + 1.0 / n2);

    mr1 -= mr2; mi1 -= mi2;
    t2 = mr1 * mr1 + mi1 * mi1;
    t2 /= serr2;

    *t2val = t2;
    *pval = Num_f2p(t2, 2, 2 * n1 - 1) * 2;
 //   *pval = Num_f2p(t2, 2, n1 - 1);
    *th = atan2(mi1, mr1);
}

// 1: sample, 2:ref
// ver2: F-test
void
Num_ctvar(float *p1, float *q1, int n1, float *p2, float *q2, int n2, float *fval, float *pval, float *th)
{
    float   f, mr1, mi1, mr2, mi2, v1, v2;
    int     df1, df2;

    Num_cavevar(p1, q1, n1, &mr1, &mi1, &v1);
    Num_cavevar(p2, q2, n2, &mr2, &mi2, &v2);

    df1 = n1 - 1;
    df2 = n2 - 1;
    f = v1 / v2;

    *fval = f;
    *pval = Num_f2p(f, df1, df2) * 0.5;   // two-sided
//    *pval = Num_f2p(f, df1, df2) * 0.3535;  // one-sided ????
    if (v1 > v2) {
        *th = M_PI * 0.1;   // pos (red)
    } else {
        *th = M_PI * 0.9;   // neg (blue)
    }
}

// === Gaussian Mixture Model ===
float
Num_nrml_dist(float x, float m, float sd)
{
	return exp(-((x - m) * (x - m)) / (2 * sd * sd));
}

Num_gmm_param *
Num_alloc_gmm(int n)
{
	Num_gmm_param *g;
	g = (Num_gmm_param *)malloc(sizeof(Num_gmm_param));
	g->n = n;
	g->mn = (float *)malloc(sizeof(float) * n);
	g->sd = (float *)malloc(sizeof(float) * n);
	g->ht = (float *)malloc(sizeof(float) * n);
	return g;
}

void
Num_free_gmm(Num_gmm_param * g)
{
	if (g) {
		if (g->mn) free(g->mn);
		if (g->sd) free(g->sd);
		if (g->ht) free(g->ht);
	}
	free(g);
}

int
Num_gmm(Num_data *d, Num_gmm_param *g)	// *g is initial / result, return val: # iteration, -1 for error return
{
	int			i, j;
	int			iter, MAX_ITER = 500;
	int			N, K;
	float		*yest;
	float		**resp;
	float		sum, den, area, sd;
	float		err, prev;

// init dim
	N = d->n;
	K = g->n;

// alloc work area etc
	yest = (float *)malloc(sizeof(float) * N);
	resp = (float **)malloc(sizeof(float *) * K);
	for (j = 0; j < K; j++) {
		resp[j] = (float *)malloc(sizeof(float) * N);
	}

// ===== EM loop
	err = prev = 0;
	for (iter = 0; iter < MAX_ITER; iter++) {
	// 2. E step
		// total
		for (i = 0; i < N; i++) {
			yest[i] = 0;
			for (j = 0; j < K; j++) {
				yest[i] += g->ht[j] * Num_nrml_dist(d->x[i], g->mn[j], g->sd[j]);
			}
		}
		// responsibilities
		for (j = 0; j < K; j++) {
			for (i = 0; i < N; i++) {
				if (yest[i] != 0) {
					resp[j][i] = g->ht[j] * Num_nrml_dist(d->x[i], g->mn[j], g->sd[j]) / yest[i];
				} else {	// this won't happen for gaussian p
					resp[j][i] = 0;
				}
			}
		}

		// 3. M step
		for (j = 0; j < K; j++) {
			// area of each peak
			den = 0;
			sum = 0;
			for (i = 0; i < N; i++) {
				den += resp[j][i] * d->y[i];
				sum += resp[j][i] * yest[i];
			}
			g->ht[j] *= den / sum;

			// mean
			sum = 0;
			for (i = 0; i < N; i++) {
				sum += resp[j][i] * d->x[i] * d->y[i];
			}
			g->mn[j] = sum / den;

			// sd
			sum = 0;
			for (i = 0; i < N; i++) {
				sd = fabs(d->x[i] - g->mn[j]) / g->sd[j];
				if (sd > 3.0) continue;
				sum += resp[j][i] * (d->x[i] - g->mn[j]) * (d->x[i] - g->mn[j]) * d->y[i];
			}
			g->sd[j] = sqrt(sum / den);
		}

		// 5. termination condition
		area = sd = 0;
		for (i = 0; i < N; i++) {
			sd += (d->y[i] - yest[i]) * (d->y[i] - yest[i]);
			if (area < yest[i]) {
				area = yest[i];
			}
		}
		sd = sqrt(sd / N);
		err = sd / area;
	//	printf("%d %f\n", iter, sd / area);
		if (fabs(err - prev) < 1e-8) break;
		prev = err;

	}
// ===== EM loop
	
	free(yest);
	for (j = 0; j < K; j++) {
		free(resp[j]);
	}
	free(resp);

	if (iter == MAX_ITER) {
		return -1;
	} else {
		return iter;
	}
}

#undef  MAXIT
#undef  EPS
#undef  FPMIN
#undef  SIGN